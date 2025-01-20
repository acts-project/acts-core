// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterError.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/VoidFitterComponents.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <type_traits>

namespace Acts {

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
  using candidate_container_t =
      typename std::vector<typename track_container_t::TrackStateProxy>;
  using TrackProxy = typename track_container_t::TrackProxy;
  using TrackStateProxy = typename track_container_t::TrackStateProxy;

  using BranchStopperResult = CombinatorialKalmanFilterBranchStopperResult;

  using Calibrator = typename KalmanFitterExtensions<traj_t>::Calibrator;
  using Updater = typename KalmanFitterExtensions<traj_t>::Updater;
  using MeasurementSelector =
      Delegate<Result<std::pair<typename candidate_container_t::iterator,
                                typename candidate_container_t::iterator>>(
          candidate_container_t& trackStates, bool&, const Logger&)>;
  using BranchStopper =
      Delegate<BranchStopperResult(const TrackProxy&, const TrackStateProxy&)>;

  /// The Calibrator is a dedicated calibration algorithm that allows to
  /// calibrate measurements using track information, this could be e.g. sagging
  /// for wires, module deformations, etc.
  Calibrator calibrator{
      DelegateFuncTag<detail::voidFitterCalibrator<traj_t>>{}};

  /// The updater incorporates measurement information into the track parameters
  Updater updater{DelegateFuncTag<detail::voidFitterUpdater<traj_t>>{}};

  /// The measurement selector is called during the filtering by the Actor.
  MeasurementSelector measurementSelector{
      DelegateFuncTag<voidMeasurementSelector>{}};

  /// The branch stopper is called during the filtering by the Actor.
  BranchStopper branchStopper{DelegateFuncTag<voidBranchStopper>{}};

 private:
  /// Default measurement selector which will return all measurements
  /// @param candidates Measurement track state candidates
  static Result<std::pair<typename std::vector<TrackStateProxy>::iterator,
                          typename std::vector<TrackStateProxy>::iterator>>
  voidMeasurementSelector(typename std::vector<TrackStateProxy>& candidates,
                          bool& /*isOutlier*/, const Logger& /*logger*/) {
    return std::pair{candidates.begin(), candidates.end()};
  };

  /// Default branch stopper which will never stop
  /// @return false
  static BranchStopperResult voidBranchStopper(
      const TrackProxy& /*track*/, const TrackStateProxy& /*trackState*/) {
    return BranchStopperResult::Continue;
  }
};

/// Delegate type that retrieves a range of source links to for a given surface
/// to be processed by the CKF
template <typename source_link_iterator_t>
using SourceLinkAccessorDelegate =
    Delegate<std::pair<source_link_iterator_t, source_link_iterator_t>(
        const Surface&)>;

/// expected max number of track states that are expected to be added by
/// stateCandidateCreator
/// @note if the number of states exceeds this number dynamic memory allocation will occur.
///       the number is chosen to yield a container size of 64 bytes.
static constexpr std::size_t s_maxBranchesPerSurface = 10;

namespace CkfTypes {

template <typename T>
using BranchVector = boost::container::small_vector<T, s_maxBranchesPerSurface>;

}  // namespace CkfTypes

/// Combined options for the combinatorial Kalman filter.
///
/// @tparam source_link_iterator_t Type of the source link iterator
/// @tparam track_container_t Type of the track container
template <typename source_link_iterator_t, typename track_container_t>
struct CombinatorialKalmanFilterOptions {
  using SourceLinkIterator = source_link_iterator_t;
  using SourceLinkAccessor = SourceLinkAccessorDelegate<source_link_iterator_t>;

  using TrackStateContainerBackend =
      typename track_container_t::TrackStateContainerBackend;
  using TrackStateProxy = typename track_container_t::TrackStateProxy;

  /// PropagatorOptions with context
  ///
  /// @param gctx The geometry context for this track finding/fitting
  /// @param mctx The magnetic context for this track finding/fitting
  /// @param cctx The calibration context for this track finding/fitting
  /// @param accessor_ The source link accessor
  /// @param extensions_ The extension struct
  /// @param pOptions The plain propagator options
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  CombinatorialKalmanFilterOptions(
      const GeometryContext& gctx, const MagneticFieldContext& mctx,
      std::reference_wrapper<const CalibrationContext> cctx,
      SourceLinkAccessor accessor_,
      CombinatorialKalmanFilterExtensions<track_container_t> extensions_,
      const PropagatorPlainOptions& pOptions, bool mScattering = true,
      bool eLoss = true)
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        sourceLinkAccessor(std::move(accessor_)),
        extensions(extensions_),
        propagatorPlainOptions(pOptions),
        multipleScattering(mScattering),
        energyLoss(eLoss) {}

  /// Contexts are required and the options must not be default-constructible.
  CombinatorialKalmanFilterOptions() = delete;

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  /// The source link accessor
  SourceLinkAccessor sourceLinkAccessor;

  /// The filter extensions
  CombinatorialKalmanFilterExtensions<track_container_t> extensions;

  /// The trivial propagator options
  PropagatorPlainOptions propagatorPlainOptions;

  /// The target surface
  /// @note This is useful if the filtering should be terminated at a
  ///       certain surface
  const Surface* targetSurface = nullptr;

  using BoundState = std::tuple<BoundTrackParameters, BoundMatrix, double>;

  /// Delegate definition to create track states for selected measurements
  ///
  /// @note expected to iterator over the given sourceLink range,
  ///       select measurements, and create track states for
  ///       which new tips are to be created, more over the outlier
  ///       flag should be set for states that are outlier.
  ///
  /// @param geoContext The current geometry context
  /// @param calibrationContext pointer to the current calibration context
  /// @param surface the surface at which new track states are to be created
  /// @param boundState the current bound state of the trajectory
  /// @param slBegin Begin iterator for sourceLinks
  /// @param slEnd End iterator for sourceLinks
  /// @param prevTip Index pointing at previous trajectory state (i.e. tip)
  /// @param bufferTrajectory a temporary trajectory which can be used to create temporary track states
  /// @param trackStateCandidates a temporary buffer that can be used to collect track states
  /// @param trajectory the trajectory to which the new states are to be added
  /// @param logger a logger for messages
  using TrackStateCandidateCreator =
      Delegate<Result<CkfTypes::BranchVector<TrackIndexType>>(
          const GeometryContext& geoContext,
          const CalibrationContext& calibrationContext, const Surface& surface,
          const BoundState& boundState, source_link_iterator_t slBegin,
          source_link_iterator_t slEnd, TrackIndexType prevTip,
          TrackStateContainerBackend& bufferTrajectory,
          std::vector<TrackStateProxy>& trackStateCandidates,
          TrackStateContainerBackend& trajectory, const Logger& logger)>;

  /// The delegate to create new track states.
  TrackStateCandidateCreator trackStateCandidateCreator;

  /// Whether to consider multiple scattering.
  bool multipleScattering = true;

  /// Whether to consider energy loss.
  bool energyLoss = true;

  /// Skip the pre propagation call. This effectively skips the first surface
  /// @note This is useful if the first surface should not be considered in a second reverse pass
  bool skipPrePropagationUpdate = false;
};

template <typename track_container_t>
struct CombinatorialKalmanFilterResult {
  using TrackStateContainerBackend =
      typename track_container_t::TrackStateContainerBackend;
  using TrackProxy = typename track_container_t::TrackProxy;
  using TrackStateProxy = typename track_container_t::TrackStateProxy;

  /// The track container to store the found tracks
  track_container_t* tracks{nullptr};

  /// Fitted states that the actor has handled.
  TrackStateContainerBackend* trackStates{nullptr};

  /// Indices into `tracks` which mark active branches
  std::vector<TrackProxy> activeBranches;

  /// Indices into `tracks` which mark active branches
  std::vector<TrackProxy> collectedTracks;

  /// Track state candidates buffer
  std::vector<TrackStateProxy> trackStateCandidates;

  /// Indicator if track finding has been done
  bool finished = false;

  /// Last encountered error
  Result<void> lastError{Result<void>::success()};

  /// Path limit aborter
  PathLimitReached pathLimitReached;
};

/// Combinatorial Kalman filter to find tracks.
///
/// @tparam propagator_t Type of the propagator
///
/// The CombinatorialKalmanFilter contains an Actor and a Sequencer sub-class.
/// The Sequencer has to be part of the Navigator of the Propagator in order to
/// initialize and provide the measurement surfaces.
///
/// The Actor is part of the Propagation call and does the Kalman update.
/// Updater and Calibrator are given to the Actor for further use:
/// - The Updater is the implemented kalman updater formalism, it
///   runs via a visitor pattern through the measurements.
///
/// Measurements are not required to be ordered for the
/// CombinatorialKalmanFilter, measurement ordering needs to be figured out by
/// the navigation of the propagator.
///
/// The void components are provided mainly for unit testing.
///
template <typename propagator_t, typename track_container_t>
class CombinatorialKalmanFilter {
 public:
  /// Default constructor is deleted
  CombinatorialKalmanFilter() = delete;
  /// Constructor from arguments
  CombinatorialKalmanFilter(propagator_t pPropagator,
                            std::unique_ptr<const Logger> _logger =
                                getDefaultLogger("CKF", Logging::INFO))
      : m_propagator(std::move(pPropagator)),
        m_logger(std::move(_logger)),
        m_actorLogger{m_logger->cloneWithSuffix("Actor")},
        m_updaterLogger{m_logger->cloneWithSuffix("Updater")} {}

 private:
  using BoundState = std::tuple<BoundTrackParameters, BoundMatrix, double>;
  using TrackStateContainerBackend =
      typename track_container_t::TrackStateContainerBackend;
  using TrackProxy = typename track_container_t::TrackProxy;
  using TrackStateProxy = typename track_container_t::TrackStateProxy;

  /// The propagator for the transport and material update
  propagator_t m_propagator;

  std::unique_ptr<const Logger> m_logger;
  std::shared_ptr<const Logger> m_actorLogger;
  std::shared_ptr<const Logger> m_updaterLogger;

  const Logger& logger() const { return *m_logger; }

  struct DefaultTrackStateCreator {
    typename CombinatorialKalmanFilterExtensions<track_container_t>::Calibrator
        calibrator;
    typename CombinatorialKalmanFilterExtensions<
        track_container_t>::MeasurementSelector measurementSelector;

    /// Create track states for selected measurements given by the source links
    ///
    /// @param gctx The current geometry context
    /// @param calibrationContext pointer to the current calibration context
    /// @param surface the surface the sourceLinks are associated to
    /// @param boundState Bound state from the propagation on this surface
    /// @param slBegin Begin iterator for sourceLinks
    /// @param slEnd End iterator for sourceLinks
    /// @param prevTip Index pointing at previous trajectory state (i.e. tip)
    /// @param bufferTrajectory a buffer for temporary candidate track states
    /// @param trackStateCandidates a buffer for temporary track state proxies for candidates
    /// @param trajectory the trajectory to which new track states for selected measurements will be added
    /// @param logger the logger for messages.
    template <typename source_link_iterator_t>
    Result<CkfTypes::BranchVector<TrackIndexType>> createSourceLinkTrackStates(
        const GeometryContext& gctx,
        const CalibrationContext& calibrationContext,
        [[maybe_unused]] const Surface& surface, const BoundState& boundState,
        source_link_iterator_t slBegin, source_link_iterator_t slEnd,
        TrackIndexType prevTip, TrackStateContainerBackend& bufferTrajectory,
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
        // CAREFUL! This trackstate has a previous index that is not in this
        // MultiTrajectory Visiting backwards from this track state will
        // fail!
        auto ts = bufferTrajectory.makeTrackState(mask, prevTip);

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
  };

  /// @brief Propagator Actor plugin for the CombinatorialKalmanFilter
  ///
  /// @tparam source_link_accessor_t The type of source link accessor
  /// @tparam parameters_t The type of parameters used for "local" parameters.
  ///
  /// The CombinatorialKalmanFilter Actor does not rely on the measurements to
  /// be sorted along the track.
  template <typename source_link_accessor_t, typename parameters_t>
  class Actor {
   public:
    using BoundState = std::tuple<parameters_t, BoundMatrix, double>;
    using CurvilinearState =
        std::tuple<CurvilinearTrackParameters, BoundMatrix, double>;
    /// Broadcast the result_type
    using result_type = CombinatorialKalmanFilterResult<track_container_t>;

    using BranchStopperResult = CombinatorialKalmanFilterBranchStopperResult;

    /// The target surface aborter
    SurfaceReached targetReached{std::numeric_limits<double>::lowest()};

    /// Whether to consider multiple scattering.
    bool multipleScattering = true;

    /// Whether to consider energy loss.
    bool energyLoss = true;

    /// Skip the pre propagation call. This effectively skips the first surface
    bool skipPrePropagationUpdate = false;

    /// Calibration context for the finding run
    const CalibrationContext* calibrationContextPtr{nullptr};

    /// @brief CombinatorialKalmanFilter actor operation
    ///
    /// @tparam propagator_state_t Type of the Propagator state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper is the stepper in use
    /// @param navigator is the navigator in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void act(propagator_state_t& state, const stepper_t& stepper,
             const navigator_t& navigator, result_type& result,
             const Logger& /*logger*/) const {
      assert(result.trackStates && "No MultiTrajectory set");

      if (result.finished) {
        return;
      }

      if (state.stage == PropagatorStage::prePropagation &&
          skipPrePropagationUpdate) {
        ACTS_VERBOSE("Skip pre-propagation update (first surface)");
        return;
      }

      ACTS_VERBOSE("CombinatorialKalmanFilter step");

      assert(!result.activeBranches.empty() && "No active branches");

      // Initialize path limit reached aborter
      if (result.pathLimitReached.internalLimit ==
          std::numeric_limits<double>::max()) {
        detail::setupLoopProtection(state, stepper, result.pathLimitReached,
                                    true, logger());
      }

      // Update:
      // - Waiting for a current surface
      if (auto surface = navigator.currentSurface(state.navigation);
          surface != nullptr) {
        // There are three scenarios:
        // 1) The surface is in the measurement map
        // -> Select source links
        // -> Perform the kalman update for selected non-outlier source links
        // -> Add track states in multitrajectory. Multiple states mean branch
        // splitting.
        // -> Call branch stopper to justify each branch
        // -> If there is non-outlier state, update stepper information
        // 2) The surface is not in the measurement map but with material or is
        // an active surface
        // -> Add a hole or passive material state in multitrajectory
        // -> Call branch stopper to justify the branch
        // 3) The surface is neither in the measurement map nor with material
        // -> Do nothing
        ACTS_VERBOSE("Perform filter step");
        auto res = filter(surface, state, stepper, navigator, result);
        if (!res.ok()) {
          ACTS_ERROR("Error in filter: " << res.error());
          result.lastError = res.error();
        }

        if (result.finished) {
          return;
        }
      }

      assert(!result.activeBranches.empty() && "No active branches");

      const bool isEndOfWorldReached =
          endOfWorldReached.checkAbort(state, stepper, navigator, logger());
      const bool isVolumeConstraintReached = volumeConstraintAborter.checkAbort(
          state, stepper, navigator, logger());
      const bool isPathLimitReached = result.pathLimitReached.checkAbort(
          state, stepper, navigator, logger());
      const bool isTargetReached =
          targetReached.checkAbort(state, stepper, navigator, logger());
      if (isEndOfWorldReached || isVolumeConstraintReached ||
          isPathLimitReached || isTargetReached) {
        if (isEndOfWorldReached) {
          ACTS_VERBOSE("End of world reached");
        } else if (isVolumeConstraintReached) {
          ACTS_VERBOSE("Volume constraint reached");
        } else if (isPathLimitReached) {
          ACTS_VERBOSE("Path limit reached");
        } else if (isTargetReached) {
          ACTS_VERBOSE("Target surface reached");

          // Bind the parameter to the target surface
          auto res = stepper.boundState(state.stepping, *targetReached.surface);
          if (!res.ok()) {
            ACTS_ERROR("Error while acquiring bound state for target surface: "
                       << res.error() << " " << res.error().message());
            result.lastError = res.error();
          } else {
            const auto& [boundParams, jacobian, pathLength] = *res;
            auto currentBranch = result.activeBranches.back();
            // Assign the fitted parameters
            currentBranch.parameters() = boundParams.parameters();
            currentBranch.covariance() = *boundParams.covariance();
            currentBranch.setReferenceSurface(
                boundParams.referenceSurface().getSharedPtr());
          }

          stepper.releaseStepSize(state.stepping,
                                  ConstrainedStep::Type::Navigator);
        }

        // Record the active branch and remove it from the list
        storeLastActiveBranch(result);
        result.activeBranches.pop_back();

        // Reset propagation state to track state at next active branch
        reset(state, stepper, navigator, result);
      }
    }

    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    bool checkAbort(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/, const result_type& result,
                    const Logger& /*logger*/) const {
      return !result.lastError.ok() || result.finished;
    }

    /// @brief CombinatorialKalmanFilter actor operation: reset propagation
    ///
    /// @tparam propagator_state_t Type of Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper is the stepper in use
    /// @param navigator is the navigator in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void reset(propagator_state_t& state, const stepper_t& stepper,
               const navigator_t& navigator, result_type& result) const {
      if (result.activeBranches.empty()) {
        ACTS_VERBOSE("Stop CKF with " << result.collectedTracks.size()
                                      << " found tracks");
        result.finished = true;

        return;
      }

      auto currentBranch = result.activeBranches.back();
      auto currentState = currentBranch.outermostTrackState();

      ACTS_VERBOSE("Propagation jumps to branch with tip = "
                   << currentBranch.tipIndex());

      // Reset the stepping state
      stepper.resetState(state.stepping, currentState.filtered(),
                         currentState.filteredCovariance(),
                         currentState.referenceSurface(),
                         state.options.stepping.maxStepSize);

      // Reset the navigation state
      // Set targetSurface to nullptr for forward filtering
      auto navigationOptions = state.navigation.options;
      navigationOptions.startSurface = &currentState.referenceSurface();
      navigationOptions.targetSurface = nullptr;
      state.navigation = navigator.makeState(navigationOptions);
      navigator.initialize(state.navigation, stepper.position(state.stepping),
                           stepper.direction(state.stepping),
                           state.options.direction);

      // No Kalman filtering for the starting surface, but still need
      // to consider the material effects here
      materialInteractor(navigator.currentSurface(state.navigation), state,
                         stepper, navigator, MaterialUpdateStage::PostUpdate);

      // Set path limit based on loop protection
      detail::setupLoopProtection(state, stepper, result.pathLimitReached, true,
                                  logger());

      // Set path limit based on target surface
      targetReached.checkAbort(state, stepper, navigator, logger());
    }

    /// @brief CombinatorialKalmanFilter actor operation:
    /// - filtering for all measurement(s) on surface
    /// - store selected track states in multiTrajectory
    /// - update propagator state to the (last) selected track state
    ///
    /// @tparam propagator_state_t Type of the Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param surface The surface where the update happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param navigator The navigator in use
    /// @param result The mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    Result<void> filter(const Surface* surface, propagator_state_t& state,
                        const stepper_t& stepper, const navigator_t& navigator,
                        result_type& result) const {
      using PM = TrackStatePropMask;

      bool isSensitive = surface->associatedDetectorElement() != nullptr;
      bool hasMaterial = surface->surfaceMaterial() != nullptr;
      bool isMaterialOnly = hasMaterial && !isSensitive;
      bool expectMeasurements = isSensitive;

      if (isSensitive) {
        ACTS_VERBOSE("Measurement surface " << surface->geometryId()
                                            << " detected.");
      } else if (isMaterialOnly) {
        ACTS_VERBOSE("Material surface " << surface->geometryId()
                                         << " detected.");
      } else {
        ACTS_VERBOSE("Passive surface " << surface->geometryId()
                                        << " detected.");
        return Result<void>::success();
      }

      using SourceLinkRange = decltype(m_sourceLinkAccessor(*surface));
      std::optional<SourceLinkRange> slRange = std::nullopt;
      bool hasMeasurements = false;
      if (isSensitive) {
        slRange = m_sourceLinkAccessor(*surface);
        hasMeasurements = slRange->first != slRange->second;
      }
      bool isHole = isSensitive && expectMeasurements && !hasMeasurements;

      if (isHole) {
        ACTS_VERBOSE("Detected hole before measurement selection on surface "
                     << surface->geometryId());
      }

      // Transport the covariance to the surface
      if (isMaterialOnly) {
        stepper.transportCovarianceToCurvilinear(state.stepping);
      } else {
        stepper.transportCovarianceToBound(state.stepping, *surface);
      }

      // Update state and stepper with pre material effects
      materialInteractor(surface, state, stepper, navigator,
                         MaterialUpdateStage::PreUpdate);

      // Bind the transported state to the current surface
      auto boundStateRes = stepper.boundState(state.stepping, *surface, false);
      if (!boundStateRes.ok()) {
        return boundStateRes.error();
      }
      auto& boundState = *boundStateRes;
      auto& [boundParams, jacobian, pathLength] = boundState;
      boundParams.covariance() = state.stepping.cov;

      auto currentBranch = result.activeBranches.back();
      TrackIndexType prevTip = currentBranch.tipIndex();

      // Create trackstates for all source links (will be filtered later)
      using TrackStatesResult =
          Acts::Result<CkfTypes::BranchVector<TrackIndexType>>;
      TrackStatesResult tsRes = TrackStatesResult::success({});
      if (hasMeasurements) {
        auto [slBegin, slEnd] = *slRange;

        tsRes = trackStateCandidateCreator(
            state.geoContext, *calibrationContextPtr, *surface, boundState,
            slBegin, slEnd, prevTip, *result.trackStates,
            result.trackStateCandidates, *result.trackStates, logger());
        if (!tsRes.ok()) {
          ACTS_ERROR(
              "Processing of selected track states failed: " << tsRes.error());
          return tsRes.error();
        }
      }
      const CkfTypes::BranchVector<TrackIndexType>& newTrackStateList = *tsRes;

      if (!newTrackStateList.empty()) {
        Result<unsigned int> procRes =
            processNewTrackStates(state.geoContext, newTrackStateList, result);
        if (!procRes.ok()) {
          ACTS_ERROR("Processing of selected track states failed: "
                     << procRes.error());
          return procRes.error();
        }
        unsigned int nBranchesOnSurface = *procRes;

        if (nBranchesOnSurface == 0) {
          ACTS_VERBOSE("All branches on surface " << surface->geometryId()
                                                  << " have been stopped");

          reset(state, stepper, navigator, result);

          return Result<void>::success();
        }

        // `currentBranch` is invalidated after `processNewTrackStates`
        currentBranch = result.activeBranches.back();
        prevTip = currentBranch.tipIndex();
      } else {
        if (expectMeasurements) {
          ACTS_VERBOSE("Detected hole after measurement selection on surface "
                       << surface->geometryId());
        }

        auto stateMask = PM::Predicted | PM::Jacobian;

        // Add a hole or material track state to the multitrajectory
        TrackIndexType currentTip = addNonSourcelinkState(
            stateMask, boundState, result, expectMeasurements, prevTip);
        currentBranch.tipIndex() = currentTip;
        auto currentState = currentBranch.outermostTrackState();
        if (expectMeasurements) {
          currentBranch.nHoles()++;
        }

        BranchStopperResult branchStopperResult =
            m_extensions.branchStopper(currentBranch, currentState);

        // Check the branch
        if (branchStopperResult == BranchStopperResult::Continue) {
          // Remembered the active branch and its state
        } else {
          // No branch on this surface
          if (branchStopperResult == BranchStopperResult::StopAndKeep) {
            storeLastActiveBranch(result);
          }
          // Remove the branch from list
          result.activeBranches.pop_back();

          // Branch on the surface has been stopped - reset
          ACTS_VERBOSE("Branch on surface " << surface->geometryId()
                                            << " has been stopped");

          reset(state, stepper, navigator, result);

          return Result<void>::success();
        }
      }

      auto currentState = currentBranch.outermostTrackState();

      if (currentState.typeFlags().test(TrackStateFlag::OutlierFlag)) {
        // We don't need to update the stepper given an outlier state
        ACTS_VERBOSE("Outlier state detected on surface "
                     << surface->geometryId());
      } else if (currentState.typeFlags().test(
                     TrackStateFlag::MeasurementFlag)) {
        // If there are measurement track states on this surface
        // Update stepping state using filtered parameters of last track
        // state on this surface
        stepper.update(state.stepping,
                       MultiTrajectoryHelpers::freeFiltered(
                           state.options.geoContext, currentState),
                       currentState.filtered(),
                       currentState.filteredCovariance(), *surface);
        ACTS_VERBOSE("Stepping state is updated with filtered parameter:");
        ACTS_VERBOSE("-> " << currentState.filtered().transpose()
                           << " of track state with tip = "
                           << currentState.index());
      }

      // Update state and stepper with post material effects
      materialInteractor(surface, state, stepper, navigator,
                         MaterialUpdateStage::PostUpdate);

      return Result<void>::success();
    }

    /// Process new, incompomplete track states and set the filtered state
    ///
    /// @note will process the given list of new states, run the updater
    ///     or share the predicted state for states flagged as outliers
    ///     and add them to the list of active branches
    ///
    /// @param gctx The geometry context for this track finding/fitting
    /// @param newTrackStateList index list of new track states
    /// @param result which contains among others the new states, and the list of active branches
    /// @return the number of newly added branches or an error
    Result<unsigned int> processNewTrackStates(
        const Acts::GeometryContext& gctx,
        const CkfTypes::BranchVector<TrackIndexType>& newTrackStateList,
        result_type& result) const {
      using PM = TrackStatePropMask;

      unsigned int nBranchesOnSurface = 0;

      auto rootBranch = result.activeBranches.back();

      // Build the new branches by forking the root branch. Reverse the order
      // to process the best candidate first
      CkfTypes::BranchVector<TrackProxy> newBranches;
      for (auto it = newTrackStateList.rbegin(); it != newTrackStateList.rend();
           ++it) {
        // Keep the root branch as the first branch, make a copy for the
        // others
        auto newBranch = (it == newTrackStateList.rbegin())
                             ? rootBranch
                             : rootBranch.shallowCopy();
        newBranch.tipIndex() = *it;
        newBranches.push_back(newBranch);
      }

      // Remove the root branch
      result.activeBranches.pop_back();

      // Update and select from the new branches
      for (TrackProxy newBranch : newBranches) {
        auto trackState = newBranch.outermostTrackState();
        TrackStateType typeFlags = trackState.typeFlags();

        if (typeFlags.test(TrackStateFlag::OutlierFlag)) {
          // No Kalman update for outlier
          // Set the filtered parameter index to be the same with predicted
          // parameter
          trackState.shareFrom(PM::Predicted, PM::Filtered);
          // Increment number of outliers
          newBranch.nOutliers()++;
        } else if (typeFlags.test(TrackStateFlag::MeasurementFlag)) {
          // Kalman update
          auto updateRes =
              m_extensions.updater(gctx, trackState, *updaterLogger);
          if (!updateRes.ok()) {
            ACTS_ERROR("Update step failed: " << updateRes.error());
            return updateRes.error();
          }
          ACTS_VERBOSE("Appended measurement track state with tip = "
                       << newBranch.tipIndex());
          // Set the measurement flag
          typeFlags.set(TrackStateFlag::MeasurementFlag);
          // Increment number of measurements
          newBranch.nMeasurements()++;
          newBranch.nDoF() += trackState.calibratedSize();
          newBranch.chi2() += trackState.chi2();
        } else {
          ACTS_WARNING("Cannot handle this track state flags");
          continue;
        }

        result.activeBranches.push_back(newBranch);

        BranchStopperResult branchStopperResult =
            m_extensions.branchStopper(newBranch, trackState);

        // Check if need to stop this branch
        if (branchStopperResult == BranchStopperResult::Continue) {
          // Record the number of branches on surface
          nBranchesOnSurface++;
        } else {
          // Record the number of stopped branches
          if (branchStopperResult == BranchStopperResult::StopAndKeep) {
            storeLastActiveBranch(result);
          }
          // Remove the branch from list
          result.activeBranches.pop_back();
        }
      }

      return nBranchesOnSurface;
    }

    /// @brief CombinatorialKalmanFilter actor operation: add a hole or material track state
    ///
    /// @param stateMask The bitmask that instructs which components to allocate
    /// @param boundState The bound state on current surface
    /// @param result is the mutable result state object and which to leave invalid
    /// @param isSensitive The surface is sensitive or passive
    /// @param prevTip The index of the previous state
    ///
    /// @return The tip of added state
    TrackIndexType addNonSourcelinkState(TrackStatePropMask stateMask,
                                         const BoundState& boundState,
                                         result_type& result, bool isSensitive,
                                         TrackIndexType prevTip) const {
      using PM = TrackStatePropMask;

      // Add a track state
      auto trackStateProxy =
          result.trackStates->makeTrackState(stateMask, prevTip);
      ACTS_VERBOSE("Create " << (isSensitive ? "Hole" : "Material")
                             << " output track state #"
                             << trackStateProxy.index()
                             << " with mask: " << stateMask);

      const auto& [boundParams, jacobian, pathLength] = boundState;
      // Fill the track state
      trackStateProxy.predicted() = boundParams.parameters();
      trackStateProxy.predictedCovariance() = boundParams.covariance().value();
      trackStateProxy.jacobian() = jacobian;
      trackStateProxy.pathLength() = pathLength;
      // Set the surface
      trackStateProxy.setReferenceSurface(
          boundParams.referenceSurface().getSharedPtr());

      // Set the track state flags
      auto typeFlags = trackStateProxy.typeFlags();
      if (trackStateProxy.referenceSurface().surfaceMaterial() != nullptr) {
        typeFlags.set(TrackStateFlag::MaterialFlag);
      }
      typeFlags.set(TrackStateFlag::ParameterFlag);
      if (isSensitive) {
        typeFlags.set(TrackStateFlag::HoleFlag);
      }

      // Set the filtered parameter index to be the same with predicted
      // parameter
      trackStateProxy.shareFrom(PM::Predicted, PM::Filtered);

      return trackStateProxy.index();
    }

    /// @brief CombinatorialKalmanFilter actor operation: material interaction
    ///
    /// @tparam propagator_state_t is the type of Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param surface The surface where the material interaction happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param navigator The navigator in use
    /// @param updateStage The material update stage
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void materialInteractor(const Surface* surface, propagator_state_t& state,
                            const stepper_t& stepper,
                            const navigator_t& navigator,
                            const MaterialUpdateStage& updateStage) const {
      if (surface == nullptr) {
        return;
      }

      // Indicator if having material
      bool hasMaterial = false;

      if (surface->surfaceMaterial() != nullptr) {
        // Prepare relevant input particle properties
        detail::PointwiseMaterialInteraction interaction(surface, state,
                                                         stepper);
        // Evaluate the material properties
        if (interaction.evaluateMaterialSlab(state, navigator, updateStage)) {
          // Surface has material at this stage
          hasMaterial = true;

          // Evaluate the material effects
          interaction.evaluatePointwiseMaterialInteraction(multipleScattering,
                                                           energyLoss);

          // Screen out material effects info
          ACTS_VERBOSE("Material effects on surface: "
                       << surface->geometryId()
                       << " at update stage: " << updateStage << " are :");
          ACTS_VERBOSE("eLoss = "
                       << interaction.Eloss * interaction.navDir << ", "
                       << "variancePhi = " << interaction.variancePhi << ", "
                       << "varianceTheta = " << interaction.varianceTheta
                       << ", "
                       << "varianceQoverP = " << interaction.varianceQoverP);

          // Update the state and stepper with material effects
          interaction.updateState(state, stepper, addNoise);
        }
      }

      if (!hasMaterial) {
        // Screen out message
        ACTS_VERBOSE("No material effects on surface: " << surface->geometryId()
                                                        << " at update stage: "
                                                        << updateStage);
      }
    }

    void storeLastActiveBranch(result_type& result) const {
      auto currentBranch = result.activeBranches.back();
      TrackIndexType currentTip = currentBranch.tipIndex();

      ACTS_VERBOSE("Storing track "
                   << currentBranch.index() << " with tip index " << currentTip
                   << ". nMeasurements = " << currentBranch.nMeasurements()
                   << ", nOutliers = " << currentBranch.nOutliers()
                   << ", nHoles = " << currentBranch.nHoles());

      result.collectedTracks.push_back(currentBranch);
    }

    CombinatorialKalmanFilterExtensions<track_container_t> m_extensions;

    /// The source link accessor
    source_link_accessor_t m_sourceLinkAccessor;

    using SourceLinkIterator =
        decltype(std::declval<decltype(m_sourceLinkAccessor(
                     *static_cast<const Surface*>(nullptr)))>()
                     .first);

    using TrackStateCandidateCreator =
        typename CombinatorialKalmanFilterOptions<
            SourceLinkIterator, track_container_t>::TrackStateCandidateCreator;

    /// the stateCandidator to be used
    /// @note will be set to a default trackStateCandidateCreator or the one
    //        provided via the extension
    TrackStateCandidateCreator trackStateCandidateCreator;

    /// End of world aborter
    EndOfWorldReached endOfWorldReached;

    /// Volume constraint aborter
    VolumeConstraintAborter volumeConstraintAborter;

    /// Actor logger instance
    const Logger* actorLogger{nullptr};
    /// Updater logger instance
    const Logger* updaterLogger{nullptr};

    const Logger& logger() const { return *actorLogger; }
  };

  /// Void path limit reached aborter to replace the default since the path
  /// limit is handled in the CKF actor internally.
  struct StubPathLimitReached {
    double internalLimit{};

    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    bool checkAbort(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/,
                    const Logger& /*logger*/) const {
      return false;
    }
  };

 public:
  /// Combinatorial Kalman Filter implementation, calls the Kalman filter
  ///
  /// @tparam source_link_iterator_t Type of the source link iterator
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param initialParameters The initial track parameters
  /// @param tfOptions CombinatorialKalmanFilterOptions steering the track
  ///                  finding
  /// @param trackContainer Track container in which to store the results
  /// @param rootBranch The track to be used as the root branch
  ///
  /// @note The input measurements are given in the form of @c SourceLinks.
  ///       It's @c calibrator_t's job to turn them into calibrated measurements
  ///       used in the track finding.
  ///
  /// @return a container of track finding result for all the initial track
  /// parameters
  template <typename source_link_iterator_t, typename start_parameters_t,
            typename parameters_t = BoundTrackParameters>
  auto findTracks(const start_parameters_t& initialParameters,
                  const CombinatorialKalmanFilterOptions<
                      source_link_iterator_t, track_container_t>& tfOptions,
                  track_container_t& trackContainer,
                  typename track_container_t::TrackProxy rootBranch) const
      -> Result<std::vector<
          typename std::decay_t<decltype(trackContainer)>::TrackProxy>> {
    using SourceLinkAccessor =
        SourceLinkAccessorDelegate<source_link_iterator_t>;

    // Create the ActorList
    using CombinatorialKalmanFilterActor =
        Actor<SourceLinkAccessor, parameters_t>;
    using Actors = ActorList<CombinatorialKalmanFilterActor>;

    // Create relevant options for the propagation options
    using PropagatorOptions = typename propagator_t::template Options<Actors>;
    PropagatorOptions propOptions(tfOptions.geoContext,
                                  tfOptions.magFieldContext);

    // Set the trivial propagator options
    propOptions.setPlainOptions(tfOptions.propagatorPlainOptions);

    // Catch the actor
    auto& combKalmanActor =
        propOptions.actorList.template get<CombinatorialKalmanFilterActor>();
    combKalmanActor.targetReached.surface = tfOptions.targetSurface;
    combKalmanActor.multipleScattering = tfOptions.multipleScattering;
    combKalmanActor.energyLoss = tfOptions.energyLoss;
    combKalmanActor.skipPrePropagationUpdate =
        tfOptions.skipPrePropagationUpdate;
    combKalmanActor.actorLogger = m_actorLogger.get();
    combKalmanActor.updaterLogger = m_updaterLogger.get();
    combKalmanActor.calibrationContextPtr = &tfOptions.calibrationContext.get();

    // copy source link accessor, calibrator and measurement selector
    combKalmanActor.m_sourceLinkAccessor = tfOptions.sourceLinkAccessor;
    combKalmanActor.m_extensions = tfOptions.extensions;
    combKalmanActor.trackStateCandidateCreator =
        tfOptions.trackStateCandidateCreator;
    DefaultTrackStateCreator defaultTrackStateCreator;
    // connect a default state candidate creator if no state candidate creator
    // was provided via the extension
    if (!combKalmanActor.trackStateCandidateCreator.connected()) {
      defaultTrackStateCreator.calibrator = tfOptions.extensions.calibrator;
      defaultTrackStateCreator.measurementSelector =
          tfOptions.extensions.measurementSelector;
      combKalmanActor.trackStateCandidateCreator.template connect<
          &DefaultTrackStateCreator::template createSourceLinkTrackStates<
              source_link_iterator_t>>(&defaultTrackStateCreator);
    }

    auto propState =
        m_propagator.template makeState<start_parameters_t, PropagatorOptions,
                                        StubPathLimitReached>(initialParameters,
                                                              propOptions);

    auto initResult =
        m_propagator
            .template initialize<decltype(propState), StubPathLimitReached>(
                propState);
    if (!initResult.ok()) {
      ACTS_ERROR("Propagation initialization failed: " << initResult.error());
      return initResult.error();
    }

    auto& r =
        propState
            .template get<CombinatorialKalmanFilterResult<track_container_t>>();
    r.tracks = &trackContainer;
    r.trackStates = &trackContainer.trackStateContainer();

    r.activeBranches.push_back(rootBranch);

    auto propagationResult = m_propagator.propagate(propState);

    auto result = m_propagator.makeResult(
        std::move(propState), propagationResult, propOptions, false);

    if (!result.ok()) {
      ACTS_ERROR("Propagation failed: " << result.error() << " "
                                        << result.error().message()
                                        << " with the initial parameters: \n"
                                        << initialParameters.parameters());
      return result.error();
    }

    auto& propRes = *result;

    // Get the result of the CombinatorialKalmanFilter
    auto combKalmanResult =
        std::move(propRes.template get<
                  CombinatorialKalmanFilterResult<track_container_t>>());

    Result<void> error = combKalmanResult.lastError;
    if (error.ok() && !combKalmanResult.finished) {
      error = Result<void>(
          CombinatorialKalmanFilterError::PropagationReachesMaxSteps);
    }
    if (!error.ok()) {
      ACTS_ERROR("CombinatorialKalmanFilter failed: "
                 << combKalmanResult.lastError.error() << " "
                 << combKalmanResult.lastError.error().message()
                 << " with the initial parameters: "
                 << initialParameters.parameters().transpose());
      return error.error();
    }

    return std::move(combKalmanResult.collectedTracks);
  }

  /// Combinatorial Kalman Filter implementation, calls the Kalman filter
  ///
  /// @tparam source_link_iterator_t Type of the source link iterator
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param initialParameters The initial track parameters
  /// @param tfOptions CombinatorialKalmanFilterOptions steering the track
  ///                  finding
  /// @param trackContainer Track container in which to store the results
  /// @note The input measurements are given in the form of @c SourceLinks.
  ///       It's @c calibrator_t's job to turn them into calibrated measurements
  ///       used in the track finding.
  ///
  /// @return a container of track finding result for all the initial track
  /// parameters
  template <typename source_link_iterator_t, typename start_parameters_t,
            typename parameters_t = BoundTrackParameters>
  auto findTracks(const start_parameters_t& initialParameters,
                  const CombinatorialKalmanFilterOptions<
                      source_link_iterator_t, track_container_t>& tfOptions,
                  track_container_t& trackContainer) const
      -> Result<std::vector<
          typename std::decay_t<decltype(trackContainer)>::TrackProxy>> {
    auto rootBranch = trackContainer.makeTrack();
    return findTracks(initialParameters, tfOptions, trackContainer, rootBranch);
  }
};

}  // namespace Acts
