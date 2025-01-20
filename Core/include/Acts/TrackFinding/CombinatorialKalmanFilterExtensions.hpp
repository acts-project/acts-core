// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/Result.hpp"

#include <vector>

namespace Acts {

namespace CkfTypes {

/// expected max number of track states that are expected to be added by
/// stateCandidateCreator
/// @note if the number of states exceeds this number dynamic memory allocation will occur.
///       the number is chosen to yield a container size of 64 bytes.

static constexpr std::size_t s_maxBranchesPerSurface = 10;
template <typename T>
using BranchVector = boost::container::small_vector<T, s_maxBranchesPerSurface>;

using BoundState = std::tuple<BoundTrackParameters, BoundMatrix, double>;

}  // namespace CkfTypes

/// Return type of the `BranchStopper` delegate for the
/// CombinatorialKalmanFilter
enum class CombinatorialKalmanFilterBranchStopperResult {
  Continue,
  StopAndDrop,
  StopAndKeep,
};

/// Extension struct which holds the delegates to customize the CKF behavior
template <typename track_container_t>
struct CombinatorialKalmanFilterExtensions {
  using traj_t = typename track_container_t::TrackStateContainerBackend;
  using TrackProxy = typename track_container_t::TrackProxy;
  using TrackStateProxy = typename track_container_t::TrackStateProxy;

  using BranchStopperResult = CombinatorialKalmanFilterBranchStopperResult;

  using Updater = typename KalmanFitterExtensions<traj_t>::Updater;
  using BranchStopper =
      Delegate<BranchStopperResult(const TrackProxy&, const TrackStateProxy&)>;

  /// The updater incorporates measurement information into the track parameters
  Updater updater{DelegateFuncTag<detail::voidFitterUpdater<traj_t>>{}};

  /// The branch stopper is called during the filtering by the Actor.
  BranchStopper branchStopper{DelegateFuncTag<voidBranchStopper>{}};

  /// @brief Delegate the extension of the trajectory onto the given surface to
  ///        an external unit.
  ///
  /// @note Expected to create track states for measurements associated to the
  ///       given surface which match the given bound state. Moreover the
  ///       The "filtered" data is not expected to be set, but the outlier
  ///       flag should be set for states that are considered to be outlier.
  ///
  /// @param geoContext The current geometry context
  /// @param calibrationContext pointer to the current calibration context
  /// @param surface the surface at which new track states are to be created
  /// @param boundState the current bound state of the trajectory
  /// @param prevTip Index pointing at previous trajectory state (i.e. tip)
  /// @param trackStateCandidates a temporary buffer that can be used to collect track states
  /// @param trajectory the trajectory to which the new states are to be added
  /// @param logger a logger for messages
  /// @return indices of new track states which extend the trajectory given by prevTip

  using TrackStateCreator =
      Delegate<Result<CkfTypes::BranchVector<TrackIndexType>>(
          const GeometryContext& geoContext,
          const CalibrationContext& calibrationContext, const Surface& surface,
          const CkfTypes::BoundState& boundState, TrackIndexType prevTip,
          std::vector<TrackStateProxy>& trackStateCandidates,
          traj_t& trajectory, const Logger& logger)>;

  /// The delegate to create new track states.
  /// @note a reference implementation can be found in @ref TrackStateCreator
  ///   which makes uses of @ref MeasurementSelector and SourceLinkAccessor
  TrackStateCreator createTrackStates;

 private:
  /// Default branch stopper which will never stop
  /// @return false
  static BranchStopperResult voidBranchStopper(
      const TrackProxy& /*track*/, const TrackStateProxy& /*trackState*/) {
    return BranchStopperResult::Continue;
  }
};
}  // namespace Acts
