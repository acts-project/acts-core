// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/VolumeMaterialMapper.hpp"

#include "Acts/EventData/NeutralTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/ApproachDescriptor.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Material/AccumulatedVolumeMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/InterpolatedMaterialMap.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialGridHelper.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/PropagatorError.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/Propagator/VolumeCollector.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/BinAdjustmentVolume.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <iosfwd>
#include <stdexcept>
#include <tuple>

namespace {
using EAxis = Acts::detail::EquidistantAxis;
using Grid2D =
    Acts::detail::Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis>;
using Grid3D =
    Acts::detail::Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis, EAxis>;
using MaterialGrid2D =
    Acts::detail::Grid<Acts::Material::ParametersVector, EAxis, EAxis>;
using MaterialGrid3D =
    Acts::detail::Grid<Acts::Material::ParametersVector, EAxis, EAxis, EAxis>;

}  // namespace

Acts::VolumeMaterialMapper::VolumeMaterialMapper(
    const Config& cfg, StraightLinePropagator propagator,
    std::unique_ptr<const Logger> slogger)
    : m_cfg(cfg),
      m_propagator(std::move(propagator)),
      m_logger(std::move(slogger)) {}

Acts::VolumeMaterialMapper::State Acts::VolumeMaterialMapper::createState(
    const GeometryContext& gctx, const MagneticFieldContext& mctx,
    const TrackingGeometry& tGeometry) const {
  // Parse the geometry and find all surfaces with material proxies
  auto world = tGeometry.highestTrackingVolume();

  // The Surface material mapping state
  State mState(gctx, mctx);
  resolveMaterialVolume(mState, *world);
  collectMaterialSurfaces(mState, *world);
  return mState;
}

void Acts::VolumeMaterialMapper::resolveMaterialVolume(
    State& mState, const TrackingVolume& tVolume) const {
  ACTS_VERBOSE("Checking volume '" << tVolume.volumeName()
                                   << "' for material surfaces.")

  ACTS_VERBOSE("- Insert Volume ...");
  checkAndInsert(mState, tVolume);

  // Step down into the sub volume
  if (tVolume.confinedVolumes()) {
    ACTS_VERBOSE("- Check children volume ...");
    for (auto& sVolume : tVolume.confinedVolumes()->arrayObjects()) {
      // Recursive call
      resolveMaterialVolume(mState, *sVolume);
    }
  }
  if (!tVolume.denseVolumes().empty()) {
    for (auto& sVolume : tVolume.denseVolumes()) {
      // Recursive call
      resolveMaterialVolume(mState, *sVolume);
    }
  }
}

void Acts::VolumeMaterialMapper::checkAndInsert(
    State& mState, const TrackingVolume& volume) const {
  auto volumeMaterial = volume.volumeMaterial();
  // Check if the volume has a proxy
  if (volumeMaterial != nullptr) {
    auto geoID = volume.geometryId();
    size_t volumeID = geoID.volume();
    ACTS_DEBUG("Material volume found with volumeID " << volumeID);
    ACTS_DEBUG("       - ID is " << geoID);

    RecordedMaterialVolumePoint mat;
    mState.recordedMaterial[geoID] = mat;

    // We need a dynamic_cast to either a volume material proxy or
    // proper surface material
    auto psm = dynamic_cast<const ProtoVolumeMaterial*>(volumeMaterial);
    // Get the bin utility: try proxy material first
    const BinUtility* bu = (psm != nullptr) ? (&psm->binUtility()) : nullptr;
    if (bu != nullptr) {
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - (proto) binning is " << *bu);
      // Now update
      BinUtility buAdjusted = adjustBinUtility(*bu, volume);
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - adjusted binning is " << buAdjusted);
      mState.materialBin[geoID] = buAdjusted;
      return;
    }
    // Second attempt: binned material
    auto bmp = dynamic_cast<
        const InterpolatedMaterialMap<MaterialMapper<MaterialGrid3D>>*>(
        volumeMaterial);
    if (bmp != nullptr) {
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - binning is " << *bu);
      mState.materialBin[geoID] = *bu;
      return;
    } else {
      // Create a homogeneous type of material
      ACTS_DEBUG("       - this is homogeneous material.");
      BinUtility buHomogeneous;
      mState.materialBin[geoID] = buHomogeneous;
      return;
    }
  }
}

void Acts::VolumeMaterialMapper::collectMaterialSurfaces(
    State& mState, const TrackingVolume& tVolume) const {
  ACTS_VERBOSE("Checking volume '" << tVolume.volumeName()
                                   << "' for material surfaces.")

  ACTS_VERBOSE("- boundary surfaces ...");
  // Check the boundary surfaces
  for (auto& bSurface : tVolume.boundarySurfaces()) {
    if (bSurface->surfaceRepresentation().surfaceMaterial() != nullptr) {
      mState.surfaceMaterial[bSurface->surfaceRepresentation().geometryId()] =
          bSurface->surfaceRepresentation().surfaceMaterialSharedPtr();
    }
  }

  ACTS_VERBOSE("- confined layers ...");
  // Check the confined layers
  if (tVolume.confinedLayers() != nullptr) {
    for (auto& cLayer : tVolume.confinedLayers()->arrayObjects()) {
      // Take only layers that are not navigation layers
      if (cLayer->layerType() != navigation) {
        // Check the representing surface
        if (cLayer->surfaceRepresentation().surfaceMaterial() != nullptr) {
          mState.surfaceMaterial[cLayer->surfaceRepresentation().geometryId()] =
              cLayer->surfaceRepresentation().surfaceMaterialSharedPtr();
        }
        // Get the approach surfaces if present
        if (cLayer->approachDescriptor() != nullptr) {
          for (auto& aSurface :
               cLayer->approachDescriptor()->containedSurfaces()) {
            if (aSurface != nullptr) {
              if (aSurface->surfaceMaterial() != nullptr) {
                mState.surfaceMaterial[aSurface->geometryId()] =
                    aSurface->surfaceMaterialSharedPtr();
              }
            }
          }
        }
        // Get the sensitive surface is present
        if (cLayer->surfaceArray() != nullptr) {
          // Sensitive surface loop
          for (auto& sSurface : cLayer->surfaceArray()->surfaces()) {
            if (sSurface != nullptr) {
              if (sSurface->surfaceMaterial() != nullptr) {
                mState.surfaceMaterial[sSurface->geometryId()] =
                    sSurface->surfaceMaterialSharedPtr();
              }
            }
          }
        }
      }
    }
  }
  // Step down into the sub volume
  if (tVolume.confinedVolumes()) {
    for (auto& sVolume : tVolume.confinedVolumes()->arrayObjects()) {
      // Recursive call
      collectMaterialSurfaces(mState, *sVolume);
    }
  }
}

void Acts::VolumeMaterialMapper::createExtraHits(
    RecordedMaterialVolumePoint& matPoint, Acts::MaterialSlab properties,
    Vector3 position, Vector3 direction) const {
  std::vector<Acts::Vector3> extraPosition;
  std::vector<Acts::Vector3> extraRemainderPositions;

  int volumeStep = floor(properties.thickness() / m_cfg.mappingStep);
  float remainder = properties.thickness() - m_cfg.mappingStep * volumeStep;
  properties.scaleThickness(m_cfg.mappingStep / properties.thickness());
  direction = direction * (m_cfg.mappingStep / direction.norm());

  for (int extraStep = 0; extraStep < volumeStep; extraStep++) {
    // Create additional extrapolated points for the grid mapping
    extraPosition.push_back(position + extraStep * direction);
  }
  matPoint.push_back(std::pair(properties, extraPosition));

  if (remainder > 0) {
    // adjust the thickness of the last extrapolated step
    properties.scaleThickness(remainder / properties.thickness());
    extraRemainderPositions.push_back(position + volumeStep * direction);
    matPoint.push_back(std::pair(properties, extraRemainderPositions));
  }
}

void Acts::VolumeMaterialMapper::finalizeMaps(State& mState) const {
  // iterate over the volumes
  for (auto& recMaterial : mState.recordedMaterial) {
    ACTS_DEBUG("Create the material for volume  " << recMaterial.first);
    if (mState.materialBin[recMaterial.first].dimensions() == 0) {
      // Accumulate all the recorded material onto a signle point
      ACTS_DEBUG("Homogeneous material volume");
      Acts::AccumulatedVolumeMaterial homogeneousAccumulation;
      for (const auto& rm : recMaterial.second) {
        homogeneousAccumulation.accumulate(rm.first);
      }
      Acts::Material mat = homogeneousAccumulation.average();
      mState.volumeMaterial[recMaterial.first] =
          std::make_unique<HomogeneousVolumeMaterial>(std::move(mat));
    } else if (mState.materialBin[recMaterial.first].dimensions() == 2) {
      // Accumulate all the recorded material onto a grid
      ACTS_DEBUG("Grid material volume");
      std::function<Acts::Vector2(Acts::Vector3)> transfoGlobalToLocal;
      Grid2D Grid = createGrid2D(mState.materialBin[recMaterial.first],
                                 transfoGlobalToLocal);
      MaterialGrid2D matGrid =
          mapMaterialPoints(Grid, recMaterial.second, transfoGlobalToLocal);
      MaterialMapper<MaterialGrid2D> matMap(transfoGlobalToLocal, matGrid);
      mState.volumeMaterial[recMaterial.first] = std::make_unique<
          InterpolatedMaterialMap<MaterialMapper<MaterialGrid2D>>>(
          std::move(matMap), mState.materialBin[recMaterial.first]);
    } else if (mState.materialBin[recMaterial.first].dimensions() == 3) {
      // Accumulate all the recorded material onto a grid
      ACTS_DEBUG("Grid material volume");
      std::function<Acts::Vector3(Acts::Vector3)> transfoGlobalToLocal;
      Grid3D Grid = createGrid3D(mState.materialBin[recMaterial.first],
                                 transfoGlobalToLocal);
      MaterialGrid3D matGrid =
          mapMaterialPoints(Grid, recMaterial.second, transfoGlobalToLocal);
      MaterialMapper<MaterialGrid3D> matMap(transfoGlobalToLocal, matGrid);
      mState.volumeMaterial[recMaterial.first] = std::make_unique<
          InterpolatedMaterialMap<MaterialMapper<MaterialGrid3D>>>(
          std::move(matMap), mState.materialBin[recMaterial.first]);
    } else {
      throw std::invalid_argument(
          "Incorrect bin dimension, only 0, 2 and 3 are accepted");
    }
  }
}

void Acts::VolumeMaterialMapper::mapMaterialTrack(
    State& mState, RecordedMaterialTrack& mTrack) const {
  using VectorHelpers::makeVector4;

  // Neutral curvilinear parameters
  NeutralCurvilinearTrackParameters start(makeVector4(mTrack.first.first, 0),
                                          mTrack.first.second,
                                          1 / mTrack.first.second.norm());

  // Prepare Action list and abort list
  using BoundSurfaceCollector = SurfaceCollector<BoundSurfaceSelector>;
  using MaterialVolumeCollector = VolumeCollector<MaterialVolumeSelector>;
  using ActionList = ActionList<BoundSurfaceCollector, MaterialVolumeCollector>;
  using AbortList = AbortList<EndOfWorldReached>;

  auto propLogger = getDefaultLogger("Propagator", Logging::INFO);
  PropagatorOptions<ActionList, AbortList> options(
      mState.geoContext, mState.magFieldContext, LoggerWrapper{*propLogger});

  // Now collect the material volume by using the straight line propagator
  const auto& result = m_propagator.propagate(start, options).value();
  auto mcResult = result.get<BoundSurfaceCollector::result_type>();
  auto mvcResult = result.get<MaterialVolumeCollector::result_type>();

  auto mappingSurfaces = mcResult.collected;
  auto mappingVolumes = mvcResult.collected;

  // Retrieve the recorded material from the recorded material track
  auto& rMaterial = mTrack.second.materialInteractions;
  ACTS_VERBOSE("Retrieved " << rMaterial.size()
                            << " recorded material steps to map.")

  // These should be mapped onto the mapping surfaces found
  ACTS_VERBOSE("Found     " << mappingVolumes.size()
                            << " mapping volumes for this track.");
  ACTS_VERBOSE("Mapping volumes are :")
  for (auto& mVolumes : mappingVolumes) {
    ACTS_VERBOSE(" - Volume : " << mVolumes.volume->geometryId()
                                << " at position = (" << mVolumes.position.x()
                                << ", " << mVolumes.position.y() << ", "
                                << mVolumes.position.z() << ")");

    mappingVolumes.push_back(mVolumes);
  }
  // Run the mapping process, i.e. take the recorded material and map it
  // onto the mapping volume:
  auto rmIter = rMaterial.begin();
  auto sfIter = mappingSurfaces.begin();
  auto volIter = mappingVolumes.begin();

  // Use those to minimize the lookup
  GeometryIdentifier lastID = GeometryIdentifier();
  GeometryIdentifier currentID = GeometryIdentifier();
  auto currentRecMaterial = mState.recordedMaterial.end();

  // store end position of the last material slab
  Acts::Vector3 lastPositionEnd = {0, 0, 0};
  Acts::Vector3 direction = {0, 0, 0};

  if (volIter != mappingVolumes.end()) {
    lastPositionEnd = volIter->position;
  }

  // loop over all the material hit in the track or until there no more volume
  // to map onto
  while (rmIter != rMaterial.end() && volIter != mappingVolumes.end()) {
    if (volIter != mappingVolumes.end() &&
        !volIter->volume->inside(rmIter->position)) {
      // Check if the material point is past the entry point to the current
      // volume
      double distVol = (volIter->position - mTrack.first.first).norm();
      double distMat = (rmIter->position - mTrack.first.first).norm();
      if (distMat - distVol > s_epsilon) {
        // Switch to next material volume
        ++volIter;
      }
    }
    if (volIter != mappingVolumes.end() &&
        volIter->volume->inside(rmIter->position)) {
      currentID = volIter->volume->geometryId();
      if (not(currentID == lastID)) {
        // Let's (re-)assess the information
        lastID = currentID;
        lastPositionEnd = volIter->position;
        currentRecMaterial = mState.recordedMaterial.find(currentID);
      }
      // If the curent volume has a ProtoVolumeMaterial
      // and the material hit has a non 0 thickness
      if (currentRecMaterial != mState.recordedMaterial.end() &&
          rmIter->materialSlab.thickness() > 0) {
        // check if there is vacuum between this material point and the last one
        float vacuumThickness = (rmIter->position - lastPositionEnd).norm();
        if (vacuumThickness > s_epsilon) {
          auto properties = Acts::MaterialSlab(vacuumThickness);
          // creat vacuum hits
          createExtraHits(currentRecMaterial->second, properties,
                          lastPositionEnd, direction);
        }
        // determine the position of the last material slab using the track
        // direction
        direction = rmIter->direction;
        direction =
            direction * (rmIter->materialSlab.thickness() / direction.norm());
        lastPositionEnd = rmIter->position + direction;
        // create additional material point
        createExtraHits(currentRecMaterial->second, rmIter->materialSlab,
                        rmIter->position, direction);
      }

      // check if we have reached the end of the volume or the last hit of the
      // track.
      if ((rmIter + 1) == rMaterial.end() ||
          !volIter->volume->inside((rmIter + 1)->position)) {
        // find the boundary surface corresponding to the end of the volume
        while (sfIter != mappingSurfaces.end()) {
          if (sfIter->surface->geometryId().volume() == lastID.volume() ||
              ((volIter + 1) != mappingVolumes.end() &&
               sfIter->surface->geometryId().volume() ==
                   (volIter + 1)->volume->geometryId().volume())) {
            double distVol = (volIter->position - mTrack.first.first).norm();
            double distSur = (sfIter->position - mTrack.first.first).norm();
            if (distSur - distVol > s_epsilon) {
              float vacuumThickness =
                  (sfIter->position - lastPositionEnd).norm();
              // if the last material slab stop before the boundary surface
              // create vacuum hits
              if (vacuumThickness > s_epsilon) {
                auto properties = Acts::MaterialSlab(vacuumThickness);
                createExtraHits(currentRecMaterial->second, properties,
                                lastPositionEnd, direction);
                lastPositionEnd = sfIter->position;
              }
              break;
            }
          }
          sfIter++;
        }
      }
      rmIter->volume = volIter->volume;
    }
    ++rmIter;
  }
}
