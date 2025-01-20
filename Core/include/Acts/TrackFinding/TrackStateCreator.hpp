// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// for definitions of Calibrator, MeasurementSelector
#include "Acts/TrackFinding/CombinatorialKalmanFilterExtensions.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"

namespace Acts {

/// @brief Create track states for selected measurements associated to a surface.
///
/// - First get a source link range covering relevant measurements associated to
///   the given surface. This task is delegated to a SourceLinkAccessor.
/// - Then create temporary track states for all measurements defined
///   by a source link range, calibrate the measurements and fill the
///   the calibrated data of these track states using a dedicated calibrator
/// - The measurement selection is delegated to a dedicated measurement
///   selector.
/// - Finally add branches to the given trajectory for the selected, temporary
///   track states. The track states of these branches still lack the filtered
///    data which is to be filled by the next stage e.g. the
///    CombinatorialKalmanFilter.
/// All track states, the temporary track states and track states for selected
/// measurements, are created in the given trajectory. The resulting container
/// may become big. Thus, it is advisable to copy selected tracks and their
/// track states to a separate container after each track finding step.
///
template <typename source_link_iterator_t, typename track_container_t>
struct TrackStateCreator {
  using TrackStatesResult =
      Acts::Result<CkfTypes::BranchVector<TrackIndexType>>;
  using TrackStateContainerBackend =
      typename track_container_t::TrackStateContainerBackend;
  using TrackProxy = typename track_container_t::TrackProxy;
  using TrackStateProxy = typename track_container_t::TrackStateProxy;
  using BoundState = std::tuple<BoundTrackParameters, BoundMatrix, double>;
  using candidate_container_t =
      typename std::vector<typename track_container_t::TrackStateProxy>;

  // delegate definition to get source link ranges for a surface
  using SourceLinkAccessor =
      Delegate<std::pair<source_link_iterator_t, source_link_iterator_t>(
          const Surface&)>;

  // delegate to get calibrted measurements from a source link iterator
  using Calibrator =
      typename KalmanFitterExtensions<TrackStateContainerBackend>::Calibrator;

  // delegate to select measurements from a track state range
  using MeasurementSelector =
      Delegate<Result<std::pair<typename candidate_container_t::iterator,
                                typename candidate_container_t::iterator>>(
          candidate_container_t& trackStates, bool&, const Logger&)>;

  /// The source link accessor will return an source link range for a surface
  /// which link to the associated measurements.
  SourceLinkAccessor sourceLinkAccessor;

  /// The Calibrator is a dedicated calibration algorithm that allows to
  /// calibrate measurements using track information, this could be e.g. sagging
  /// for wires, module deformations, etc.
  Calibrator calibrator{DelegateFuncTag<
      detail::voidFitterCalibrator<TrackStateContainerBackend>>{}};

  MeasurementSelector measurementSelector{
      DelegateFuncTag<voidMeasurementSelector>{}};

 public:
  /// @brief extend the trajectory onto the given surface.
  ///
  /// @param gctx The geometry context to be used for this task
  /// @param calibrationContext The calibration context used to fill the calibrated data
  /// @param surface The surface onto which the trajectory is extended
  /// @param boundState the predicted bound state on the given surface
  /// @param prevTip the tip of the trajectory which is to be extended
  /// @param trackStateCandidates a temporary buffer which can be used to
  ///        to keep track of newly created temporary track states.
  /// @param trajectory the trajectory to be extended.
  /// @param logger a logger for messages.
  ///
  /// @return a list of indices of newly created track states which extend the
  ///    trajectory onto the given surface and match the bound state, or an
  ///    error.
  ///
  /// Extend or branch the trajectory onto the given surface. This may create
  /// new track states using measurements which match the predicted bound state.
  /// This may create multiple branches. The new track states still miss the
  /// "filtered" data.
  Result<CkfTypes::BranchVector<TrackIndexType>> createTrackStates(
      const GeometryContext& gctx, const CalibrationContext& calibrationContext,
      [[maybe_unused]] const Surface& surface, const BoundState& boundState,
      TrackIndexType prevTip,
      std::vector<TrackStateProxy>& trackStateCandidates,
      TrackStateContainerBackend& trajectory, const Logger& logger) const {
    TrackStatesResult tsRes = TrackStatesResult::success({});
    using SourceLinkRange = decltype(sourceLinkAccessor(surface));
    SourceLinkRange slRange = sourceLinkAccessor(surface);
    if (slRange.first != slRange.second) {
      auto [slBegin, slEnd] = slRange;
      tsRes = createSourceLinkTrackStates(
          gctx, calibrationContext, surface, boundState, slBegin, slEnd,
          prevTip, trackStateCandidates, trajectory, logger);
    }
    return tsRes;
  }

  /// Create track states for selected measurements given by the source links
  ///
  /// @param gctx The current geometry context
  /// @param calibrationContext pointer to the current calibration context
  /// @param surface the surface the sourceLinks are associated to
  /// @param boundState Bound state from the propagation on this surface
  /// @param slBegin Begin iterator for sourceLinks
  /// @param slEnd End iterator for sourceLinks
  /// @param prevTip Index pointing at previous trajectory state (i.e. tip)
  /// @param trackStateCandidates a temporary buffer which can be used to
  ///        to keep track of newly created temporary track states.
  /// @param trajectory the trajectory to which new track states for selected measurements will be added
  /// @param logger the logger for messages.
  Result<CkfTypes::BranchVector<TrackIndexType>> createSourceLinkTrackStates(
      const GeometryContext& gctx, const CalibrationContext& calibrationContext,
      [[maybe_unused]] const Surface& surface, const BoundState& boundState,
      source_link_iterator_t slBegin, source_link_iterator_t slEnd,
      TrackIndexType prevTip,
      std::vector<TrackStateProxy>& trackStateCandidates,
      TrackStateContainerBackend& trajectory, const Logger& logger) const {
    using PM = TrackStatePropMask;

    using ResultTrackStateList =
        Acts::Result<CkfTypes::BranchVector<TrackIndexType>>;
    ResultTrackStateList resultTrackStateList{
        CkfTypes::BranchVector<TrackIndexType>()};
    const auto& [boundParams, jacobian, pathLength] = boundState;

    trackStateCandidates.clear();
    if constexpr (std::ranges::random_access_range<source_link_iterator_t>) {
      trackStateCandidates.reserve(std::distance(slBegin, slEnd));
    }

    // Calibrate all the source links on the surface since the selection has
    // to be done based on calibrated measurement
    for (auto it = slBegin; it != slEnd; ++it) {
      // get the source link
      const auto sourceLink = *it;

      // prepare the track state
      PM mask = PM::Predicted | PM::Jacobian | PM::Calibrated;
      if (it != slBegin) {
        // not the first TrackState, only need uncalibrated and calibrated
        mask = PM::Calibrated;
      }

      ACTS_VERBOSE("Create temp track state with mask: " << mask);
      // Temporary and final track states are created in the same
      // trajectory, which could lead to very large containers.

      // CAREFUL! This trackstate has a previous index that is not in this
      // MultiTrajectory Visiting backwards from this track state will
      // fail!
      auto ts = trajectory.makeTrackState(mask, prevTip);

      if (it == slBegin) {
        // only set these for first
        ts.predicted() = boundParams.parameters();
        if (boundParams.covariance()) {
          ts.predictedCovariance() = *boundParams.covariance();
        }
        ts.jacobian() = jacobian;
      } else {
        // subsequent track states can reuse
        auto& first = trackStateCandidates.front();
        ts.shareFrom(first, PM::Predicted);
        ts.shareFrom(first, PM::Jacobian);
      }

      ts.pathLength() = pathLength;
      ts.setReferenceSurface(boundParams.referenceSurface().getSharedPtr());

      // now calibrate the track state
      calibrator(gctx, calibrationContext, sourceLink, ts);

      trackStateCandidates.push_back(ts);
    }

    bool isOutlier = false;
    Result<std::pair<typename std::vector<TrackStateProxy>::iterator,
                     typename std::vector<TrackStateProxy>::iterator>>
        selectorResult =
            measurementSelector(trackStateCandidates, isOutlier, logger);
    if (!selectorResult.ok()) {
      ACTS_ERROR("Selection of calibrated measurements failed: "
                 << selectorResult.error().message());
      resultTrackStateList =
          ResultTrackStateList::failure(selectorResult.error());
    } else {
      auto selectedTrackStateRange = *selectorResult;
      resultTrackStateList = processSelectedTrackStates(
          selectedTrackStateRange.first, selectedTrackStateRange.second,
          trajectory, isOutlier, logger);
    }

    return resultTrackStateList;
  }

  /// Create track states for the given trajectory from candidate track states
  ///
  /// @param begin begin iterator of the list of candidate track states
  /// @param end end iterator of the list of candidate track states
  /// @param trackStates the trajectory to which the new track states are added
  /// @param isOutlier true if the candidate(s) is(are) an outlier(s).
  /// @param logger the logger for messages
  Result<CkfTypes::BranchVector<TrackIndexType>> processSelectedTrackStates(
      typename std::vector<TrackStateProxy>::const_iterator begin,
      typename std::vector<TrackStateProxy>::const_iterator end,
      TrackStateContainerBackend& trackStates, bool isOutlier,
      const Logger& logger) const {
    using PM = TrackStatePropMask;

    using ResultTrackStateList =
        Acts::Result<CkfTypes::BranchVector<TrackIndexType>>;
    ResultTrackStateList resultTrackStateList{
        CkfTypes::BranchVector<TrackIndexType>()};
    CkfTypes::BranchVector<TrackIndexType>& trackStateList =
        *resultTrackStateList;
    trackStateList.reserve(end - begin);

    std::optional<TrackStateProxy> firstTrackState{std::nullopt};
    for (auto it = begin; it != end; ++it) {
      auto& candidateTrackState = *it;

      PM mask = PM::Predicted | PM::Filtered | PM::Jacobian | PM::Calibrated;
      if (it != begin) {
        // subsequent track states don't need storage for these as they will
        // be shared
        mask &= ~PM::Predicted & ~PM::Jacobian;
      }
      if (isOutlier) {
        // outlier won't have separate filtered parameters
        mask &= ~PM::Filtered;
      }

      // copy this trackstate into fitted states MultiTrajectory
      auto trackState =
          trackStates.makeTrackState(mask, candidateTrackState.previous());
      ACTS_VERBOSE("Create SourceLink output track state #"
                   << trackState.index() << " with mask: " << mask);

      if (it != begin) {
        // assign indices pointing to first track state
        trackState.shareFrom(*firstTrackState, PM::Predicted);
        trackState.shareFrom(*firstTrackState, PM::Jacobian);
      } else {
        firstTrackState = trackState;
      }

      // either copy ALL or everything except for predicted and jacobian
      trackState.copyFrom(candidateTrackState, mask, false);

      auto typeFlags = trackState.typeFlags();
      typeFlags.set(TrackStateFlag::ParameterFlag);
      typeFlags.set(TrackStateFlag::MeasurementFlag);
      if (trackState.referenceSurface().surfaceMaterial() != nullptr) {
        typeFlags.set(TrackStateFlag::MaterialFlag);
      }
      if (isOutlier) {
        // propagate information that this is an outlier state
        ACTS_VERBOSE(
            "Creating outlier track state with tip = " << trackState.index());
        typeFlags.set(TrackStateFlag::OutlierFlag);
      }

      trackStateList.push_back(trackState.index());
    }
    return resultTrackStateList;
  }

  /// Default measurement selector which will return all measurements
  /// @param candidates Measurement track state candidates
  static Result<std::pair<typename std::vector<TrackStateProxy>::iterator,
                          typename std::vector<TrackStateProxy>::iterator>>
  voidMeasurementSelector(typename std::vector<TrackStateProxy>& candidates,
                          bool& /*isOutlier*/, const Logger& /*logger*/) {
    return std::pair{candidates.begin(), candidates.end()};
  };
};

}  // namespace Acts
