// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Detray/DetrayGeometryConverter.hpp"

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Navigation/PortalNavigation.hpp"
#include "Acts/Plugins/Json/DetrayJsonHelper.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <algorithm>

#include <detray/io/frontend/detector_writer.hpp>

using namespace detray;

detray::io::transform_payload Acts::DetrayGeometryConverter::convertTransform(
    const Transform3& t) {
  detray::io::transform_payload tfPayload;
  Vector3 translation = t.translation();
  tfPayload.tr = {translation.x(), translation.y(), translation.z()};
  RotationMatrix3 rotation = t.rotation().transpose();
  tfPayload.rot = {rotation(0, 0), rotation(0, 1), rotation(0, 2),
                   rotation(1, 0), rotation(1, 1), rotation(1, 2),
                   rotation(2, 0), rotation(2, 1), rotation(2, 2)};
  return tfPayload;
}

detray::io::mask_payload Acts::DetrayGeometryConverter::convertMask(
    const Acts::SurfaceBounds& bounds, bool portal) {
  detray::io::mask_payload maskPayload;
  auto [shape, boundaries] = DetrayJsonHelper::maskFromBounds(bounds, portal);
  maskPayload.shape = static_cast<io::mask_payload::mask_shape>(shape);
  maskPayload.boundaries = static_cast<std::vector<real_io>>(boundaries);
  // default maskPayload.volume_link

  return maskPayload;
}

detray::io::surface_payload Acts::DetrayGeometryConverter::convertSurface(
    const GeometryContext& gctx, const Surface& surface, bool portal) {
  detray::io::surface_payload surfacePayload;

  surfacePayload.transform = convertTransform(surface.transform(gctx));
  surfacePayload.source = surface.geometryId().value();
  surfacePayload.barcode = std::nullopt;
  surfacePayload.type = static_cast<detray::surface_id>(
      portal ? surface_id::e_portal
             : (surface.associatedDetectorElement() != nullptr
                    ? detray::surface_id::e_sensitive
                    : detray::surface_id::e_passive));
  surfacePayload.mask = convertMask(surface.bounds());
  return surfacePayload;
}

std::vector<detray::io::surface_payload>
Acts::DetrayGeometryConverter::convertPortal(
    DetrayConversionUtils::Cache& cCache, const GeometryContext& gctx,
    const Experimental::Portal& portal, std::size_t ip,
    const Experimental::DetectorVolume& volume,
    const std::vector<Acts::OrientedSurface>& orientedSurfaces) {
  std::vector<detray::io::surface_payload> portals{};

  const RegularSurface& surface = portal.surface();
  const auto& volumeLinks = portal.portalNavigation();

  // First assumption for outside link (along direction)
  std::size_t outside = 1u;

  // Find out if you need to take the outside or inside volume
  // for planar surfaces that's easy
  if (surface.type() != Acts::Surface::SurfaceType::Cylinder) {
    // Get the two volume center
    const auto volumeCenter = volume.transform(gctx).translation();
    const auto surfaceCenter = surface.center(gctx);
    const auto surfaceNormal = surface.normal(gctx, surfaceCenter);
    // Get the direction from the volume to the surface, correct link
    const auto volumeToSurface = surfaceCenter - volumeCenter;
    if (volumeToSurface.dot(surfaceNormal) < 0.) {
      outside = 0u;
    }
  } else {
    // This is a cylinder portal, inner cover reverses the normal
    if (ip == 3u) {
      outside = 0u;
    }
  }

  const auto& outsideLink = volumeLinks[outside];
  // Grab the corresponding volume link
  // If it is a single link, we are done
  const auto* instance = outsideLink.instance();
  // Single link cast
  auto singleLink =
      dynamic_cast<const Acts::Experimental::SingleDetectorVolumeNavigation*>(
          instance);

  auto [surfaceAdjusted, insidePointer] = orientedSurfaces[ip];
  // Assign the geometry id to the surface
  surfaceAdjusted->assignGeometryId(surface.geometryId());

  // Single link detected - just write it out, we use the oriented surface
  // in order to make sure the size is adjusted
  if (singleLink != nullptr) {
    // Single link can be written out
    std::size_t vLink = cCache.volumeIndex(singleLink->object());
    auto portalPayload = convertSurface(gctx, *surfaceAdjusted, true);
    portalPayload.mask.volume_link.link = vLink;
    portals.push_back(portalPayload);
  } else {
    // Multi link detected - 1D
    auto multiLink1D =
        dynamic_cast<const Experimental::BoundVolumesGrid1Navigation*>(
            instance);
    if (multiLink1D != nullptr) {
      // Resolve the multi link 1D
      auto boundaries =
          multiLink1D->indexedUpdater.grid.axes()[0u]->getBinEdges();
      const auto& cast = multiLink1D->indexedUpdater.casts[0u];
      const auto& transform = multiLink1D->indexedUpdater.transform;
      const auto& volumes = multiLink1D->indexedUpdater.extractor.dVolumes;

      // Apply the correction from local to global boundaries
      double gCorr = VectorHelpers::cast(transform.translation(), cast);
      std::for_each(boundaries.begin(), boundaries.end(),
                    [&gCorr](double& b) { b -= gCorr; });

      // Get the volume indices
      auto surfaceType = surfaceAdjusted->type();
      std::vector<unsigned int> vIndices = {};
      for (const auto& v : volumes) {
        vIndices.push_back(cCache.volumeIndex(v));
      }

      // Pick the surface dimension
      std::array<double, 2u> clipRange = {0., 0.};
      std::vector<double> boundValues = surfaceAdjusted->bounds().values();
      if (surfaceType == Surface::SurfaceType::Cylinder &&
          cast == AxisDirection::AxisZ) {
        double zPosition = surfaceAdjusted->center(gctx).z();
        clipRange = {
            zPosition - boundValues[CylinderBounds::BoundValues::eHalfLengthZ],
            zPosition + boundValues[CylinderBounds::BoundValues::eHalfLengthZ]};
      } else if (surfaceType == Surface::SurfaceType::Disc &&
                 cast == AxisDirection::AxisR) {
        clipRange = {boundValues[RadialBounds::BoundValues::eMinR],
                     boundValues[RadialBounds::BoundValues::eMaxR]};
      } else {
        throw std::runtime_error(
            "PortalDetrayGeometryConverter: surface type not (yet) supported "
            "for "
            "detray "
            "conversion, only cylinder and disc are currently supported.");
      }

      // Need to clip the parameter space to the surface dimension
      std::vector<unsigned int> clippedIndices = {};
      std::vector<double> clippedBoundaries = {};
      clippedBoundaries.push_back(clipRange[0u]);
      for (const auto [ib, b] : enumerate(boundaries)) {
        if (ib > 0) {
          unsigned int vI = vIndices[ib - 1u];
          double highEdge = boundaries[ib];
          if (boundaries[ib - 1] >= clipRange[1u]) {
            break;
          }
          if (highEdge <= clipRange[0u] ||
              std::abs(highEdge - clipRange[0u]) < 1e-5) {
            continue;
          }
          if (highEdge > clipRange[1u]) {
            highEdge = clipRange[1u];
          }
          clippedIndices.push_back(vI);
          clippedBoundaries.push_back(highEdge);
        }
      }
      // Interpret the clipped information
      //
      // Clipped cylinder case
      if (surfaceType == Surface::SurfaceType::Cylinder) {
        for (auto [ib, b] : enumerate(clippedBoundaries)) {
          if (ib > 0) {
            // Create sub surfaces
            std::array<double, CylinderBounds::BoundValues::eSize>
                subBoundValues = {};
            for (auto [ibv, bv] : enumerate(boundValues)) {
              subBoundValues[ibv] = bv;
            }
            subBoundValues[CylinderBounds::BoundValues::eHalfLengthZ] =
                0.5 * (b - clippedBoundaries[ib - 1u]);
            auto subBounds = std::make_shared<CylinderBounds>(subBoundValues);
            auto subTransform = Transform3::Identity();
            subTransform.pretranslate(Vector3(
                0., 0.,
                clippedBoundaries[ib - 1u] +
                    subBoundValues[CylinderBounds::BoundValues::eHalfLengthZ]));

            auto subSurface =
                Surface::makeShared<CylinderSurface>(subTransform, subBounds);
            subSurface->assignGeometryId(surface.geometryId());

            auto portalPayload = convertSurface(gctx, *subSurface, true);
            portalPayload.mask.volume_link.link = clippedIndices[ib - 1u];
            portals.push_back(portalPayload);
          }
        }
      } else {
        for (auto [ib, b] : enumerate(clippedBoundaries)) {
          if (ib > 0) {
            // Create sub surfaces
            std::array<double, RadialBounds::BoundValues::eSize>
                subBoundValues = {};
            for (auto [ibv, bv] : enumerate(boundValues)) {
              subBoundValues[ibv] = bv;
            }
            subBoundValues[RadialBounds::BoundValues::eMinR] =
                clippedBoundaries[ib - 1u];
            subBoundValues[RadialBounds::BoundValues::eMaxR] = b;
            auto subBounds = std::make_shared<RadialBounds>(subBoundValues);
            auto subSurface = Surface::makeShared<DiscSurface>(
                portal.surface().transform(gctx), subBounds);

            subSurface->assignGeometryId(surface.geometryId());
            auto portalPayload = convertSurface(gctx, *subSurface, true);
            portalPayload.mask.volume_link.link = clippedIndices[ib - 1u];
            portals.push_back(portalPayload);
          }
        }
      }

    } else {
      // Write surface with invalid link
      auto portalPayload = convertSurface(gctx, *surfaceAdjusted, true);
      using NavigationLink =
          typename DetrayHostDetector::surface_type::navigation_link;
      portalPayload.mask.volume_link.link =
          std::numeric_limits<NavigationLink>::max();

      portals.push_back(portalPayload);
    }
  }
  return portals;
}

detray::io::volume_payload Acts::DetrayGeometryConverter::convertVolume(
    DetrayConversionUtils::Cache& cCache, const GeometryContext& gctx,
    const Acts::Experimental::DetectorVolume& volume,
    const Acts::Logger& logger) {
  ACTS_DEBUG("DetrayGeometryConverter: converting volume "
             << volume.name() << " with " << volume.surfaces().size()
             << " surfaces and " << volume.portals().size() << " portals");

  detray::io::volume_payload volumePayload;
  std::size_t volumeIndex = cCache.volumeIndex(&volume);
  volumePayload.name = volume.name();
  volumePayload.index.link = volumeIndex;
  volumePayload.transform = convertTransform(volume.transform(gctx));

  // Remember the link
  cCache.volumeLinks[volume.geometryId()] = volumePayload.index.link;

  std::multimap<GeometryIdentifier, unsigned long> localSurfaceLinks;

  // iterate over surfaces and portals keeping the same surf_pd.index_in_coll
  std::size_t sIndex = 0;
  for (const auto surface : volume.surfaces()) {
    io::surface_payload surfacePayload = convertSurface(gctx, *surface, false);
    // Set the index in the collection & remember it in the cache
    surfacePayload.index_in_coll = sIndex++;
    localSurfaceLinks.insert(
        {surface->geometryId(), surfacePayload.index_in_coll.value()});
    // Set mask to volume link
    surfacePayload.mask.volume_link.link =
        volumePayload.index.link;  // link surface' mask to volume
    volumePayload.surfaces.push_back(surfacePayload);
  }

  auto orientedSurfaces =
      volume.volumeBounds().orientedSurfaces(volume.transform(gctx));

  // Iterate over portals
  int portalCounter = 0;
  for (const auto& [ip, p] : enumerate(volume.portals())) {
    auto portals =
        convertPortal(cCache, gctx, *p, ip, volume, orientedSurfaces);
    ACTS_VERBOSE(" > portal " << ip << " split into " << portals.size()
                              << " surfaces");
    GeometryIdentifier geoID = p->surface().geometryId();
    std::for_each(portals.begin(), portals.end(), [&](auto& portalPayload) {
      // Set the index in the collection & remember it in the cache
      portalPayload.index_in_coll = sIndex++;
      localSurfaceLinks.insert({geoID, portalPayload.index_in_coll.value()});
      // Add it to the surfaces
      volumePayload.surfaces.push_back(portalPayload);
      portalCounter++;
    });
  }
  cCache.localSurfaceLinks[volumeIndex] = localSurfaceLinks;
  ACTS_DEBUG(" > " << volume.portals().size()
                   << " initial ACTS portals split into final " << portalCounter
                   << " detray portals");
  ACTS_VERBOSE(" > Local surface link cache has " << localSurfaceLinks.size()
                                                  << " entries");

  return volumePayload;
}

detray::io::detector_payload Acts::DetrayGeometryConverter::convertDetector(
    DetrayConversionUtils::Cache& cCache, const GeometryContext& gctx,
    const Acts::Experimental::Detector& detector, const Acts::Logger& logger) {
  ACTS_DEBUG("DetrayGeometryConverter: converting detector"
             << detector.name() << " with " << detector.volumes().size()
             << " volumes.");

  // The detector payload which will be handed back
  detray::io::detector_payload detectorPayload;

  for (const auto volume : detector.volumes()) {
    detectorPayload.volumes.push_back(
        convertVolume(cCache, gctx, *volume, logger));
  }

  return detectorPayload;
}
