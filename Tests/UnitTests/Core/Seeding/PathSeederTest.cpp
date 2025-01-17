// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/detail/CuboidalDetectorHelper.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Seeding/PathSeeder.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/MeasurementsCreator.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <numbers>

BOOST_AUTO_TEST_SUITE(PathSeeder)

using namespace Acts;
using namespace Acts::UnitLiterals;

using Axis = Acts::Axis<AxisType::Equidistant, AxisBoundaryType::Open>;
using Grid = Acts::Grid<std::vector<SourceLink>, Axis, Axis>;

using TrackParameters = CurvilinearTrackParameters;

GeometryContext gctx;

// Parameters for the geometry
const double halfY = 10.;
const double halfZ = 10.;
const double deltaX = 10.;
const double deltaYZ = 1.;

const Vector4 trueVertex(-5., 0., 0., 0);
const std::vector<double> truePhis = {-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15};
const double trueTheta = std::numbers::pi / 2.;
const double trueQOverP = 1. / 1._GeV;

// Intersection finding to get the
// region of interest for seeding
class NoFieldIntersectionFinder {
 public:
  double m_tol = 1e-4;

  std::vector<const Surface*> m_surfaces;

  // Find the intersections along the path
  // and return them in the order of the path
  // length
  std::vector<std::pair<GeometryIdentifier, Vector2>> operator()(
      const GeometryContext& geoCtx,
      const TrackParameters& trackParameters) const {
    Vector3 position = trackParameters.position();
    Vector3 direction = trackParameters.direction();

    std::vector<std::pair<GeometryIdentifier, Vector2>> sIntersections;
    // Intersect the surfaces
    for (auto& surface : m_surfaces) {
      // Get the intersection
      auto sMultiIntersection = surface->intersect(
          geoCtx, position, direction,
          BoundaryTolerance::AbsoluteCartesian(m_tol, m_tol));

      // Take the closest
      auto closestForward = sMultiIntersection.closestForward();

      // Store if the intersection is reachable
      if (closestForward.status() == IntersectionStatus::reachable &&
          closestForward.pathLength() > 0.0) {
        sIntersections.push_back(
            {closestForward.object()->geometryId(),
             surface
                 ->globalToLocal(geoCtx, closestForward.position(),
                                 Vector3{0, 1, 0})
                 .value()});
        continue;
      }
    }
    return sIntersections;
  }
};

// A simple path width provider to set
// the grid lookup boundaries around the
// intersection point
class PathWidthProvider {
 public:
  std::pair<double, double> width;

  std::pair<double, double> operator()(
      const GeometryContext& /*gctx*/,
      const GeometryIdentifier& /*geoId*/) const {
    return width;
  }
};

// Estimator of the particle's energy,
// vertex, momentum direction at the IP
// and the direction at the first hit
class TrackEstimator {
 public:
  Vector3 m_ip;
  SourceLinkSurfaceAccessor m_surfaceAccessor;

  std::pair<TrackParameters, TrackParameters> operator()(
      const GeometryContext& geoCtx, const SourceLink& pivot) const {
    auto testSourceLink = pivot.get<detail::Test::TestSourceLink>();
    Vector3 pivot3 = m_surfaceAccessor(pivot)->localToGlobal(
        geoCtx, testSourceLink.parameters, Vector3{0, 1, 0});

    Vector3 direction = (pivot3 - m_ip).normalized();

    Vector4 ip = {m_ip.x(), m_ip.y(), m_ip.z(), 0};
    double qOverP = 1_e / 1._GeV;
    double phi = Acts::VectorHelpers::phi(direction);
    double theta = Acts::VectorHelpers::theta(direction);
    ParticleHypothesis particle = ParticleHypothesis::electron();

    TrackParameters ipParams(ip, phi, theta, qOverP, std::nullopt, particle);
    TrackParameters firstLayerParams(ip, phi, theta, qOverP, std::nullopt,
                                     particle);

    return {ipParams, firstLayerParams};
  }
};

// Construct grid with the source links
struct ConstructSourceLinkGrid {
  SourceLinkSurfaceAccessor m_surfaceAccessor;

  std::unordered_map<GeometryIdentifier, Grid> construct(
      std::vector<SourceLink> sourceLinks) {
    // Lookup table for each layer
    std::unordered_map<GeometryIdentifier, Grid> lookupTable;

    // Construct a binned grid for each layer
    for (int i : {14, 15, 16, 17}) {
      Axis xAxis(-halfY, halfY, 100);
      Axis yAxis(-halfZ, halfZ, 100);

      Grid grid(std::make_tuple(xAxis, yAxis));

      GeometryIdentifier geoId;

      geoId.setSensitive(i);

      lookupTable.insert({geoId, grid});
    }
    // Fill the grid with source links
    for (auto& sl : sourceLinks) {
      auto tsl = sl.get<detail::Test::TestSourceLink>();
      auto id = tsl.m_geometryId;

      auto bin = lookupTable.at(id).localBinsFromPosition(tsl.parameters);
      lookupTable.at(id).atLocalBins(bin).push_back(sl);
    }

    return lookupTable;
  }
};

// Construct a simple telescope detector
std::shared_ptr<Experimental::Detector> constructTelescopeDetector() {
  RotationMatrix3 rotation;
  double angle = 90_degree;
  Vector3 xPos(cos(angle), 0., sin(angle));
  Vector3 yPos(0., 1., 0.);
  Vector3 zPos(-sin(angle), 0., cos(angle));
  rotation.col(0) = xPos;
  rotation.col(1) = yPos;
  rotation.col(2) = zPos;

  // Create a bunch of surface bounds
  auto surf0bounds = std::make_unique<RectangleBounds>(halfY, halfZ);

  auto surf1bounds = std::make_unique<RectangleBounds>(halfY, halfZ);

  auto surf2bounds = std::make_unique<RectangleBounds>(halfY, halfZ);

  auto surf3bounds = std::make_unique<RectangleBounds>(halfY, halfZ);

  // Create a bunch of surfaces
  auto transform0 = Transform3::Identity();
  auto surf0 = Surface::makeShared<PlaneSurface>(
      transform0 * Transform3(rotation), std::move(surf0bounds));
  auto geoId0 = GeometryIdentifier();
  geoId0.setSensitive(1);

  auto transform1 =
      Transform3::Identity() * Translation3(Vector3(2 * deltaX, 0, 0));
  auto surf1 = Surface::makeShared<PlaneSurface>(
      transform1 * Transform3(rotation), std::move(surf1bounds));
  auto geoId1 = GeometryIdentifier();
  geoId1.setSensitive(2);

  auto transform2 =
      Transform3::Identity() * Translation3(Vector3(4 * deltaX, 0, 0));
  auto surf2 = Surface::makeShared<PlaneSurface>(
      transform2 * Transform3(rotation), std::move(surf2bounds));
  auto geoId2 = GeometryIdentifier();
  geoId2.setSensitive(3);

  auto transform3 =
      Transform3::Identity() * Translation3(Vector3(6 * deltaX, 0, 0));
  auto surf3 = Surface::makeShared<PlaneSurface>(
      transform3 * Transform3(rotation), std::move(surf3bounds));
  auto geoId3 = GeometryIdentifier();
  geoId3.setSensitive(4);

  // Create a bunch of volume bounds
  auto vol0bounds = std::make_unique<CuboidVolumeBounds>(
      deltaX, halfY + deltaYZ, halfZ + deltaYZ);

  auto vol1bounds = std::make_unique<CuboidVolumeBounds>(
      deltaX, halfY + deltaYZ, halfZ + deltaYZ);

  auto vol2bounds = std::make_unique<CuboidVolumeBounds>(
      deltaX, halfY + deltaYZ, halfZ + deltaYZ);

  auto vol3bounds = std::make_unique<CuboidVolumeBounds>(
      deltaX, halfY + deltaYZ, halfZ + deltaYZ);

  // Create a bunch of volumes
  auto vol0 = Experimental::DetectorVolumeFactory::construct(
      Experimental::defaultPortalAndSubPortalGenerator(), gctx, "vol0",
      transform0, std::move(vol0bounds), {surf0}, {},
      Experimental::tryNoVolumes(), Experimental::tryAllPortalsAndSurfaces());

  auto vol1 = Experimental::DetectorVolumeFactory::construct(
      Experimental::defaultPortalAndSubPortalGenerator(), gctx, "vol1",
      transform1, std::move(vol1bounds), {surf1}, {},
      Experimental::tryNoVolumes(), Experimental::tryAllPortalsAndSurfaces());

  auto vol2 = Experimental::DetectorVolumeFactory::construct(
      Experimental::defaultPortalAndSubPortalGenerator(), gctx, "vol2",
      transform2, std::move(vol2bounds), {surf2}, {},
      Experimental::tryNoVolumes(), Experimental::tryAllPortalsAndSurfaces());

  auto vol3 = Experimental::DetectorVolumeFactory::construct(
      Experimental::defaultPortalAndSubPortalGenerator(), gctx, "vol3",
      transform3, std::move(vol3bounds), {surf3}, {},
      Experimental::tryNoVolumes(), Experimental::tryAllPortalsAndSurfaces());

  std::vector<std::shared_ptr<Experimental::DetectorVolume>> volumes = {
      vol0, vol1, vol2, vol3};

  // Connect the volumes
  auto portalContainer = Experimental::detail::CuboidalDetectorHelper::connect(
      gctx, volumes, AxisDirection::AxisX, {}, Logging::INFO);

  // Make sure that the geometry ids are
  // independent of the potential Id generation
  // changes
  int id = 1;

  // Volume ids: 1-3
  for (auto& volume : volumes) {
    volume->assignGeometryId(Acts::GeometryIdentifier(id));
    id++;
  }
  // Intervolume portal ids: 6,7,10,11
  for (auto& volume : volumes) {
    for (auto& port : volume->portalPtrs()) {
      if (port->surface().geometryId() == Acts::GeometryIdentifier(0)) {
        port->surface().assignGeometryId(Acts::GeometryIdentifier(id));
        id++;
      }
    }
  }
  // Surface ids
  for (auto& surf : {surf0, surf1, surf2, surf3}) {
    auto geoId = GeometryIdentifier();
    geoId.setSensitive(id);
    surf->assignGeometryId(geoId);
    id++;
  }

  auto detector = Experimental::Detector::makeShared(
      "TelescopeDetector", volumes, Experimental::tryRootVolumes());

  return detector;
}

std::vector<SourceLink> createSourceLinks(
    const GeometryContext& geoCtx, const Experimental::Detector& detector) {
  NoFieldIntersectionFinder intersectionFinder;

  for (auto volume : detector.volumes()) {
    for (auto surface : volume->surfaces()) {
      intersectionFinder.m_surfaces.push_back(surface);
    }
  }

  std::vector<SourceLink> sourceLinks;
  for (double phi : truePhis) {
    TrackParameters trackParameters(trueVertex, phi, trueTheta, trueQOverP,
                                    std::nullopt,
                                    ParticleHypothesis::electron());

    auto intersections = intersectionFinder(geoCtx, trackParameters);

    SquareMatrix2 cov = SquareMatrix2::Identity();

    for (auto& [id, refPoint] : intersections) {
      detail::Test::TestSourceLink sourceLink(eBoundLoc0, eBoundLoc1, refPoint,
                                              cov, id, id.value());

      SourceLink sl{sourceLink};
      sourceLinks.push_back(sl);
    }
  }

  return sourceLinks;
}

BOOST_AUTO_TEST_CASE(PathSeederZeroField) {
  using SurfaceAccessor =
      detail::Test::Experimental::TestSourceLinkSurfaceAccessor;

  // Create detector
  auto detector = constructTelescopeDetector();

  // Create source links
  auto sourceLinks = createSourceLinks(gctx, *detector);

  // Prepare the PathSeeder
  auto pathSeederCfg = Acts::PathSeeder::Config();

  // Grid to bin the source links
  SurfaceAccessor surfaceAccessor{*detector};
  auto sourceLinkGridConstructor = ConstructSourceLinkGrid();
  sourceLinkGridConstructor.m_surfaceAccessor
      .connect<&SurfaceAccessor::operator()>(&surfaceAccessor);

  // Create the grid
  std::unordered_map<GeometryIdentifier, Grid> sourceLinkGrid =
      sourceLinkGridConstructor.construct(sourceLinks);

  // Estimator of the IP and first hit
  // parameters of the track
  TrackEstimator trackEstimator;
  trackEstimator.m_ip = Vector3(-5., 0., 0.);
  trackEstimator.m_surfaceAccessor.connect<&SurfaceAccessor::operator()>(
      &surfaceAccessor);
  pathSeederCfg.trackEstimator.connect<&TrackEstimator::operator()>(
      &trackEstimator);

  // Intersection finder
  NoFieldIntersectionFinder intersectionFinder;
  for (auto volume : detector->volumes()) {
    for (auto surface : volume->surfaces()) {
      intersectionFinder.m_surfaces.push_back(surface);
    }
  }
  pathSeederCfg.intersectionFinder
      .connect<&NoFieldIntersectionFinder::operator()>(&intersectionFinder);

  // Path width provider
  PathWidthProvider pathWidthProvider;

  pathSeederCfg.pathWidthProvider.connect<&PathWidthProvider::operator()>(
      &pathWidthProvider);

  GeometryIdentifier geoId;
  geoId.setSensitive(14);
  pathSeederCfg.refLayerIds.push_back(geoId);

  // Create the PathSeeder
  Acts::PathSeeder pathSeeder(pathSeederCfg);

  // Get the seeds
  pathWidthProvider.width = {1e-3, 1e-3};

  std::vector<Acts::PathSeeder::PathSeed> seeds;

  // SeedTreeContainer seeds;
  pathSeeder.findSeeds(gctx, sourceLinkGrid, seeds);

  // Check the seeds
  BOOST_CHECK_EQUAL(seeds.size(), 7);

  for (std::size_t i = 0; i < seeds.size(); i++) {
    auto seed = seeds.at(i);
    BOOST_CHECK_EQUAL(seed.second.size(), 4);
    BOOST_CHECK_EQUAL(seed.first.phi(), truePhis.at(i));
    BOOST_CHECK_EQUAL(seed.first.theta(), trueTheta);
    BOOST_CHECK_EQUAL(seed.first.qOverP(), trueQOverP);
  }
}

BOOST_AUTO_TEST_SUITE_END()
