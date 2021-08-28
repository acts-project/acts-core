// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <functional>
#include <unordered_map>

namespace Acts {

namespace MultiTrajectoryHelpers {

/// @brief Struct for brief trajectory summary info
///
struct TrajectoryState {
  size_t nStates = 0;
  size_t nMeasurements = 0;
  size_t nOutliers = 0;
  size_t nHoles = 0;
  double chi2Sum = 0;
  std::vector<double> measurementChi2 = {};
  std::vector<double> outlierChi2 = {};
  size_t NDF = 0;
  std::vector<unsigned int> measurementVolume = {};
  std::vector<unsigned int> measurementLayer = {};
  std::vector<unsigned int> outlierVolume = {};
  std::vector<unsigned int> outlierLayer = {};
  size_t nSharedHits = 0;
};

// Container for trajectory summary info at a specific volume
using VolumeTrajectoryStateContainer =
    std::unordered_map<GeometryIdentifier::Value, TrajectoryState>;

/// @brief Getter for global trajectory info
///
/// @tparam source_link_t Type of source link
///
/// @param multiTraj The MultiTrajectory object
/// @param entryIndex The entry index of trajectory to investigate
///
/// @return The trajectory summary info
template <typename source_link_t>
TrajectoryState trajectoryState(
    const Acts::MultiTrajectory<source_link_t>& multiTraj,
    const size_t& entryIndex) {
  TrajectoryState trajState;
  multiTraj.visitBackwards(entryIndex, [&](const auto& state) {
    // Get the volume Id of this surface
    const auto& geoID = state.referenceSurface().geometryId();
    const auto& volume = geoID.volume();
    const auto& layer = geoID.layer();
    trajState.nStates++;
    trajState.chi2Sum += state.chi2();
    trajState.NDF += state.calibratedSize();
    auto typeFlags = state.typeFlags();
    if (typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
      if (typeFlags.test(Acts::TrackStateFlag::SharedHitFlag))
        trajState.nSharedHits++;
      trajState.nMeasurements++;
      trajState.measurementChi2.push_back(state.chi2());
      trajState.measurementVolume.push_back(volume);
      trajState.measurementLayer.push_back(layer);
    } else if (typeFlags.test(Acts::TrackStateFlag::OutlierFlag)) {
      trajState.nOutliers++;
      trajState.outlierChi2.push_back(state.chi2());
      trajState.outlierVolume.push_back(volume);
      trajState.outlierLayer.push_back(layer);
    } else if (typeFlags.test(Acts::TrackStateFlag::HoleFlag)) {
      trajState.nHoles++;
    }
  });
  return trajState;
}

/// @brief Getter for trajectory info for different sub-detectors
///
/// @tparam source_link_t Type of source link
///
/// @param multiTraj The MultiTrajectory object
/// @param entryIndex The entry index of trajectory to investigate
/// track states at different sub-detectors.
/// @param volumeIds The container for sub-detector Ids
///
/// @return The trajectory summary info at different sub-detectors (i.e.
/// different volumes)
template <typename source_link_t>
VolumeTrajectoryStateContainer trajectoryState(
    const Acts::MultiTrajectory<source_link_t>& multiTraj,
    const size_t& entryIndex,
    const std::vector<GeometryIdentifier::Value>& volumeIds) {
  VolumeTrajectoryStateContainer trajStateContainer;
  multiTraj.visitBackwards(entryIndex, [&](const auto& state) {
    // Get the volume Id of this surface
    const auto& geoID = state.referenceSurface().geometryId();
    const auto& volume = geoID.volume();
    const auto& layer = geoID.layer();
    // Check if the track info for this sub-detector is requested
    auto it = std::find(volumeIds.begin(), volumeIds.end(), volume);
    if (it == volumeIds.end()) {
      return true;
    }
    // The trajectory state for this volume
    auto& trajState = trajStateContainer[volume];
    trajState.nStates++;
    trajState.chi2Sum += state.chi2();
    trajState.NDF += state.calibratedSize();
    auto typeFlags = state.typeFlags();
    if (typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
      if (typeFlags.test(Acts::TrackStateFlag::SharedHitFlag))
        trajState.nSharedHits++;
      trajState.nMeasurements++;
      trajState.measurementChi2.push_back(state.chi2());
      trajState.measurementVolume.push_back(volume);
      trajState.measurementLayer.push_back(layer);
    } else if (typeFlags.test(Acts::TrackStateFlag::OutlierFlag)) {
      trajState.nOutliers++;
      trajState.outlierChi2.push_back(state.chi2());
      trajState.outlierVolume.push_back(volume);
      trajState.outlierLayer.push_back(layer);
    } else if (typeFlags.test(Acts::TrackStateFlag::HoleFlag)) {
      trajState.nHoles++;
    }
    return true;
  });
  return trajStateContainer;
}

/// @brief Transforms the filtered parameters from a @c TrackStateProxy to free
/// parameters
///
/// @tparam track_state_proxy_t Type of the @c TrackStateProxy
/// @param [in] gctx Geometry context
/// @param [in] trackStateProxy TrackStateProxy
///
/// @return Free parameters representation of the filtered parameters
template <typename track_state_proxy_t>
FreeVector freeFiltered(const GeometryContext& gctx,
                        const track_state_proxy_t& trackStateProxy) {
  return detail::transformBoundToFreeParameters(
      trackStateProxy.referenceSurface(), gctx, trackStateProxy.filtered());
}

/// @brief Transforms the smoothed parameters from a @c TrackStateProxy to free
/// parameters
///
/// @tparam track_state_proxy_t Type of the @c TrackStateProxy
/// @param [in] gctx Geometry context
/// @param [in] trackStateProxy TrackStateProxy
///
/// @return Free parameters representation of the smoothed parameters
template <typename track_state_proxy_t>
FreeVector freeSmoothed(const GeometryContext& gctx,
                        const track_state_proxy_t& trackStateProxy) {
  return detail::transformBoundToFreeParameters(
      trackStateProxy.referenceSurface(), gctx, trackStateProxy.smoothed());
}
}  // namespace MultiTrajectoryHelpers

}  // namespace Acts
