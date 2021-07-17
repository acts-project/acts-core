// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/MeasurementsCreator.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "Acts/Tests/CommonHelpers/TestSourceLink.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/KalmanGlobalCovariance.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsAlignment/Kernel/detail/AlignmentEngine.hpp"

#include <cmath>
#include <random>
#include <string>

namespace {
using namespace Acts;
using namespace ActsAlignment;
using namespace Acts::Test;
using namespace Acts::UnitLiterals;

using StraightPropagator =
    Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
using ConstantFieldStepper = Acts::EigenStepper<>;
using ConstantFieldPropagator =
    Acts::Propagator<ConstantFieldStepper, Acts::Navigator>;

using KalmanUpdater = Acts::GainMatrixUpdater;
using KalmanSmoother = Acts::GainMatrixSmoother;
using KalmanFitter =
    Acts::KalmanFitter<ConstantFieldPropagator, KalmanUpdater, KalmanSmoother>;

// Create a test context
const GeometryContext geoCtx;
const MagneticFieldContext magCtx;
const CalibrationContext calCtx;

std::normal_distribution<double> normalDist(0., 1.);
std::default_random_engine rng(42);

///
/// @brief Contruct a telescope-like detector
///
struct TelescopeTrackingGeometry {
  /// Default constructor for the Cubit tracking geometry
  ///
  /// @param gctx the geometry context for this geometry at building time
  TelescopeTrackingGeometry(std::reference_wrapper<const GeometryContext> gctx)
      : geoContext(gctx) {
    using namespace UnitLiterals;

    // Construct the rotation
    double rotationAngle = 90_degree;
    Vector3 xPos(cos(rotationAngle), 0., sin(rotationAngle));
    Vector3 yPos(0., 1., 0.);
    Vector3 zPos(-sin(rotationAngle), 0., cos(rotationAngle));
    rotation.col(0) = xPos;
    rotation.col(1) = yPos;
    rotation.col(2) = zPos;

    // Boundaries of the surfaces
    rBounds =
        std::make_shared<const RectangleBounds>(RectangleBounds(0.1_m, 0.1_m));

    // Material of the surfaces
    MaterialSlab matProp(Acts::Test::makeSilicon(), 0.5_mm);
    surfaceMaterial = std::make_shared<HomogeneousSurfaceMaterial>(matProp);
  }

  ///
  /// Call operator to build the standard cubic tracking geometry
  ///
  std::shared_ptr<const TrackingGeometry> operator()() {
    using namespace UnitLiterals;

    // Set translation vectors
    std::vector<Vector3> translations;
    translations.reserve(6);
    translations.push_back({-500_mm, 0., 0.});
    translations.push_back({-300_mm, 0., 0.});
    translations.push_back({-100_mm, 0., 0.});
    translations.push_back({100_mm, 0., 0.});
    translations.push_back({300_mm, 0., 0.});
    translations.push_back({500_mm, 0., 0.});

    // Construct layer configs
    std::vector<CuboidVolumeBuilder::LayerConfig> lConfs;
    lConfs.reserve(6);
    unsigned int i;
    for (i = 0; i < translations.size(); i++) {
      CuboidVolumeBuilder::SurfaceConfig sConf;
      sConf.position = translations[i];
      sConf.rotation = rotation;
      sConf.rBounds = rBounds;
      sConf.surMat = surfaceMaterial;
      // The thickness to construct the associated detector element
      sConf.thickness = 1._um;
      sConf.detElementConstructor =
          [](const Transform3& trans,
             std::shared_ptr<const RectangleBounds> bounds, double thickness) {
            return new Acts::Test::DetectorElementStub(trans, bounds,
                                                       thickness);
          };
      CuboidVolumeBuilder::LayerConfig lConf;
      lConf.surfaceCfg = sConf;
      lConfs.push_back(lConf);
    }

    // Construct volume config
    CuboidVolumeBuilder::VolumeConfig vConf;
    vConf.position = {0., 0., 0.};
    vConf.length = {1.2_m, 1._m, 1._m};
    vConf.layerCfg = lConfs;
    vConf.name = "Tracker";

    // Construct volume builder config
    CuboidVolumeBuilder::Config conf;
    conf.position = {0., 0., 0.};
    conf.length = {1.2_m, 1._m, 1._m};
    conf.volumeCfg = {vConf};  // one volume

    // Build detector
    CuboidVolumeBuilder cvb;
    cvb.setConfig(conf);
    TrackingGeometryBuilder::Config tgbCfg;
    tgbCfg.trackingVolumeBuilders.push_back(
        [=](const auto& context, const auto& inner, const auto& vb) {
          return cvb.trackingVolume(context, inner, vb);
        });
    TrackingGeometryBuilder tgb(tgbCfg);
    std::shared_ptr<const TrackingGeometry> detector =
        tgb.trackingGeometry(geoCtx);

    // Build and return tracking geometry
    return detector;
  }

  RotationMatrix3 rotation = RotationMatrix3::Identity();
  std::shared_ptr<const RectangleBounds> rBounds = nullptr;
  std::shared_ptr<const ISurfaceMaterial> surfaceMaterial = nullptr;

  std::reference_wrapper<const GeometryContext> geoContext;
};

// Build detector
TelescopeTrackingGeometry geometryStore(geoCtx);
const auto geometry = geometryStore();

// Construct a straight-line propagator.
StraightPropagator makeStraightPropagator(
    std::shared_ptr<const Acts::TrackingGeometry> geo) {
  Acts::Navigator::Config cfg{geo};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);
  Acts::StraightLineStepper stepper;
  return StraightPropagator(std::move(stepper), std::move(navigator));
}

// Construct a propagator using a constant magnetic field along z.
ConstantFieldPropagator makeConstantFieldPropagator(
    std::shared_ptr<const Acts::TrackingGeometry> geo, double bz) {
  Acts::Navigator::Config cfg{geo};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);
  auto field =
      std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, bz));
  ConstantFieldStepper stepper(std::move(field));
  return ConstantFieldPropagator(std::move(stepper), std::move(navigator));
}

// Construct initial track parameters.
Acts::CurvilinearTrackParameters makeParameters() {
  // create covariance matrix from reasonable standard deviations
  Acts::BoundVector stddev;
  stddev[Acts::eBoundLoc0] = 100_um;
  stddev[Acts::eBoundLoc1] = 100_um;
  stddev[Acts::eBoundTime] = 25_ns;
  stddev[Acts::eBoundPhi] = 0.5_degree;
  stddev[Acts::eBoundTheta] = 0.5_degree;
  stddev[Acts::eBoundQOverP] = 1 / 100_GeV;
  Acts::BoundSymMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();

  auto loc0 = 0. + stddev[Acts::eBoundLoc0] * normalDist(rng);
  auto loc1 = 0. + stddev[Acts::eBoundLoc1] * normalDist(rng);
  auto t = 42_ns + stddev[Acts::eBoundTime] * normalDist(rng);
  auto phi = 0_degree + stddev[Acts::eBoundPhi] * normalDist(rng);
  auto theta = 90_degree + stddev[Acts::eBoundTheta] * normalDist(rng);
  auto qOverP = 1_e / 1_GeV + stddev[Acts::eBoundQOverP] * normalDist(rng);

  // define a track in the transverse plane along x
  Acts::Vector4 mPos4(-1_m, loc0, loc1, t);

  return Acts::CurvilinearTrackParameters(mPos4, phi, theta, 1_e / qOverP, 1_e,
                                          cov);
}

// detector resolutions
const MeasurementResolution resPixel = {MeasurementType::eLoc01,
                                        {30_um, 50_um}};
const MeasurementResolutionMap resolutions = {
    {GeometryIdentifier().setVolume(1), resPixel},
};

// simulation propagator
const auto simPropagator = makeStraightPropagator(geometry);

// reconstruction propagator and fitter
const auto kfLogger = getDefaultLogger("KalmanFilter", Logging::INFO);
const auto kfZeroPropagator = makeConstantFieldPropagator(geometry, 0_T);
const auto kfZero = KalmanFitter(kfZeroPropagator);

// alignment
const auto alignLogger = getDefaultLogger("Alignment", Logging::VERBOSE);
const auto alignZero = Alignment(kfZero);

struct KalmanFitterInputTrajectory {
  // The source links
  std::vector<TestSourceLink> sourcelinks;
  // The start parameters
  std::optional<CurvilinearTrackParameters> startParameters;
};

///
/// Function to create trajectories for kalman fitter
///
std::vector<KalmanFitterInputTrajectory> createTrajectories(
    //    const std::shared_ptr<const TrackingGeometry>& geo,
    size_t nTrajectories,
    const std::array<double, 2>& localSigma = {1000_um, 1000_um},
    const double& pSigma = 0.025_GeV) {
  std::vector<KalmanFitterInputTrajectory> trajectories;
  trajectories.reserve(nTrajectories);

  for (unsigned int iTrack = 0; iTrack < nTrajectories; iTrack++) {
    if (iTrack % 10 == 0) {
      std::cout << "Processing track: " << iTrack << "..." << std::endl;
    }
    // Set initial parameters for the particle track
    //    Vector4D mPos4(-1_m, 100_um * gauss(rng), 100_um * gauss(rng),
    //                   42_ns);
    //    Vector3 mDir(1_GeV, 0.01_GeV * gauss(rng),
    //                  0.01_GeV * gauss(rng));
    //    CurvilinearTrackParameters mStart(mPos4, mDir, 1_e / 1_GeV);

    auto start = makeParameters();
    // Launch and collect - the measurements
    auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                           resolutions, rng);

    // Extract measurements from result of propagation.
    KalmanFitterInputTrajectory traj;
    traj.startParameters = start;
    traj.sourcelinks = measurements.sourceLinks;

    // Smear the start parameters to be used as input of KF
    //    Covariance cov;
    //    cov << std::pow(localSigma[0], 2), 0., 0., 0., 0., 0., 0.,
    //        std::pow(localSigma[1], 2), 0., 0., 0., 0., 0., 0., pSigma, 0.,
    //        0., 0., 0., 0., 0., pSigma, 0., 0., 0., 0., 0., 0., 0.01, 0., 0.,
    //        0., 0., 0., 0., 1.;
    //    Vector4D rPos4(mPos4.x(), mPos4.y() + localSigma[0] * gauss(rng),
    //                   mPos4.z() + localSigma[1] * gauss(rng), 42_ns);
    //    Vector3 rDir(mDir.x(), mDir.y() + pSigma * gauss(rng),
    //                  mDir.z() + pSigma * gauss(rng));
    //    CurvilinearTrackParameters rStart(rPos4, rDir, 1_e / 1_GeV, cov);

    trajectories.push_back(std::move(traj));
  }
  return trajectories;
}
}  // namespace

///
/// @brief Unit test for KF-based alignment algorithm
///
BOOST_AUTO_TEST_CASE(ZeroFieldKalmanAlignment) {
  // Create the trajectories
  const auto& trajectories = createTrajectories(100);

  // Construct the KalmanFitter options
  KalmanFitterOptions<TestSourceLinkCalibrator, VoidOutlierFinder> kfOptions(
      geoCtx, magCtx, calCtx, TestSourceLinkCalibrator(), VoidOutlierFinder(),
      LoggerWrapper{*kfLogger}, PropagatorPlainOptions());
  // this is the default option. set anyways for consistency
  kfOptions.referenceSurface = nullptr;

  // Construct an non-updating alignment updater
  AlignedTransformUpdater voidAlignUpdater =
      [](Acts::DetectorElementBase* /*unused*/,
         const Acts::GeometryContext& /*unused*/,
         const Acts::Transform3& /*unused*/) { return true; };

  // Construct the alignment options
  AlignmentOptions<
      KalmanFitterOptions<TestSourceLinkCalibrator, VoidOutlierFinder>>
      alignOptions(kfOptions, voidAlignUpdater, LoggerWrapper{*alignLogger});
  alignOptions.maxIterations = 1;

  // Set the surfaces to be aligned
  unsigned int iSurface = 0;
  std::unordered_map<const Surface*, size_t> idxedAlignSurfaces;
  geometry->visitSurfaces([&](const Surface* surface) {
    // Missing out the forth layer
    if (surface and surface->associatedDetectorElement() and
        surface->geometryId().layer() != 8) {
      alignOptions.alignedDetElements.push_back(
          const_cast<DetectorElementBase*>(
              surface->associatedDetectorElement()));
      idxedAlignSurfaces.emplace(surface, iSurface);
      iSurface++;
    }
  });

  // The alignment mask
  const auto alignMask =
      std::bitset<Acts::eAlignmentSize>(std::string("111111"));

  // Test the method to evaluate alignment state for a single track
  const auto& inputTraj = trajectories.front();
  kfOptions.referenceSurface = &(*inputTraj.startParameters).referenceSurface();
  auto evaluateRes = alignZero.evaluateTrackAlignmentState(
      kfOptions.geoContext, inputTraj.sourcelinks, *inputTraj.startParameters,
      kfOptions, idxedAlignSurfaces, alignMask, alignOptions.logger);
  BOOST_CHECK(evaluateRes.ok());
  const auto& alignState = evaluateRes.value();
  std::cout << "Chi2/dof = " << alignState.chi2 / alignState.alignmentDof
            << std::endl;

  // Check the dimensions
  BOOST_CHECK_EQUAL(alignState.measurementDim, 12);
  BOOST_CHECK_EQUAL(alignState.trackParametersDim, 36);
  // Check the alignment dof
  BOOST_CHECK_EQUAL(alignState.alignmentDof, 30);
  BOOST_CHECK_EQUAL(alignState.alignedSurfaces.size(), 5);
  // Check the measurements covariance
  BOOST_CHECK_EQUAL(alignState.measurementCovariance.rows(), 12);
  const SymMatrix2 measCov = alignState.measurementCovariance.block<2, 2>(2, 2);
  SymMatrix2 cov2D;
  cov2D << 30_um * 30_um, 0, 0, 50_um * 50_um;
  CHECK_CLOSE_ABS(measCov, cov2D, 1e-10);
  // Check the track parameters covariance matrix. Its rows/columns scales
  // with the number of measurement states
  BOOST_CHECK_EQUAL(alignState.trackParametersCovariance.rows(), 36);
  // Check the projection matrix
  BOOST_CHECK_EQUAL(alignState.projectionMatrix.rows(), 12);
  BOOST_CHECK_EQUAL(alignState.projectionMatrix.cols(), 36);
  const ActsMatrix<2, 6> proj = alignState.projectionMatrix.block<2, 6>(0, 0);
  const ActsMatrix<2, 6> refProj = ActsMatrix<2, 6>::Identity();
  CHECK_CLOSE_ABS(proj, refProj, 1e-10);
  // Check the residual
  BOOST_CHECK_EQUAL(alignState.residual.size(), 12);
  // Check the residual covariance
  BOOST_CHECK_EQUAL(alignState.residualCovariance.rows(), 12);
  // Check the alignment to residual derivative
  BOOST_CHECK_EQUAL(alignState.alignmentToResidualDerivative.rows(), 12);
  BOOST_CHECK_EQUAL(alignState.alignmentToResidualDerivative.cols(), 30);
  // Check the chi2 derivative
  BOOST_CHECK_EQUAL(alignState.alignmentToChi2Derivative.size(), 30);
  BOOST_CHECK_EQUAL(alignState.alignmentToChi2SecondDerivative.rows(), 30);

  // Test the align method
  std::vector<std::vector<TestSourceLink>> trajCollection;
  trajCollection.reserve(100);
  std::vector<CurvilinearTrackParameters> sParametersCollection;
  sParametersCollection.reserve(100);
  for (const auto& traj : trajectories) {
    trajCollection.push_back(traj.sourcelinks);
    sParametersCollection.push_back(*traj.startParameters);
  }
  auto alignRes =
      alignZero.align(trajCollection, sParametersCollection, alignOptions);

  // BOOST_CHECK(alignRes.ok());
}
