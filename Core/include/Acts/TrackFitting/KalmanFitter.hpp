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

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackContainerFrontendConcept.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/TrackFitting/detail/KalmanUpdateHelpers.hpp"
#include "Acts/TrackFitting/detail/VoidFitterComponents.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"

#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <type_traits>

namespace Acts {

enum class KalmanFitterTargetSurfaceStrategy {
  /// Use the first trackstate to reach target surface
  first,
  /// Use the last trackstate to reach target surface
  last,
  /// Use the first or last trackstate to reach target surface depending on the
  /// distance
  firstOrLast,
};

/// Extension struct which holds delegates to customise the KF behavior
template <typename traj_t>
struct KalmanFitterExtensions {
  using TrackStateProxy = typename traj_t::TrackStateProxy;
  using ConstTrackStateProxy = typename traj_t::ConstTrackStateProxy;
  using Parameters = typename TrackStateProxy::Parameters;

  using Calibrator =
      Delegate<void(const GeometryContext&, const CalibrationContext&,
                    const SourceLink&, TrackStateProxy)>;

  using Smoother = Delegate<Result<void>(const GeometryContext&, traj_t&,
                                         std::size_t, const Logger&)>;

  using Updater = Delegate<Result<void>(const GeometryContext&, TrackStateProxy,
                                        const Logger&)>;

  using OutlierFinder = Delegate<bool(ConstTrackStateProxy)>;

  using ReverseFilteringLogic = Delegate<bool(ConstTrackStateProxy)>;

  /// The Calibrator is a dedicated calibration algorithm that allows
  /// to calibrate measurements using track information, this could be
  /// e.g. sagging for wires, module deformations, etc.
  Calibrator calibrator;

  /// The updater incorporates measurement information into the track parameters
  Updater updater;

  /// The smoother back-propagates measurement information along the track
  Smoother smoother;

  /// Determines whether a measurement is supposed to be considered as an
  /// outlier
  OutlierFinder outlierFinder;

  /// Decides whether the smoothing stage uses linearized transport or full
  /// reverse propagation
  ReverseFilteringLogic reverseFilteringLogic;

  /// Retrieves the associated surface from a source link
  SourceLinkSurfaceAccessor surfaceAccessor;

  /// Default constructor which connects the default void components
  KalmanFitterExtensions() {
    calibrator.template connect<&detail::voidFitterCalibrator<traj_t>>();
    updater.template connect<&detail::voidFitterUpdater<traj_t>>();
    smoother.template connect<&detail::voidFitterSmoother<traj_t>>();
    outlierFinder.template connect<&detail::voidOutlierFinder<traj_t>>();
    reverseFilteringLogic
        .template connect<&detail::voidReverseFilteringLogic<traj_t>>();
    surfaceAccessor.connect<&detail::voidSurfaceAccessor>();
  }
};

/// Combined options for the Kalman fitter.
///
/// @tparam traj_t The trajectory type
template <typename traj_t>
struct KalmanFitterOptions {
  /// PropagatorOptions with context.
  ///
  /// @param gctx The geometry context for this fit
  /// @param mctx The magnetic context for this fit
  /// @param cctx The calibration context for this fit
  /// @param extensions_ The KF extensions
  /// @param pOptions The plain propagator options
  /// @param rSurface The reference surface for the fit to be expressed at
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  /// @param rFiltering Whether to run filtering in reversed direction as smoothing
  /// @param rfScaling Scale factor for the covariance matrix before the backward filtering
  /// @param freeToBoundCorrection_ Correction for non-linearity effect during transform from free to bound
  KalmanFitterOptions(const GeometryContext& gctx,
                      const MagneticFieldContext& mctx,
                      std::reference_wrapper<const CalibrationContext> cctx,
                      KalmanFitterExtensions<traj_t> extensions_,
                      const PropagatorPlainOptions& pOptions,
                      const Surface* rSurface = nullptr,
                      bool mScattering = true, bool eLoss = true,
                      bool rFiltering = false, double rfScaling = 1.0,
                      const FreeToBoundCorrection& freeToBoundCorrection_ =
                          FreeToBoundCorrection(false))
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        extensions(extensions_),
        propagatorPlainOptions(pOptions),
        referenceSurface(rSurface),
        multipleScattering(mScattering),
        energyLoss(eLoss),
        reversedFiltering(rFiltering),
        reversedFilteringCovarianceScaling(rfScaling),
        freeToBoundCorrection(freeToBoundCorrection_) {}
  /// Contexts are required and the options must not be default-constructible.
  KalmanFitterOptions() = delete;

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  KalmanFitterExtensions<traj_t> extensions;

  /// The trivial propagator options
  PropagatorPlainOptions propagatorPlainOptions;

  /// The reference Surface
  const Surface* referenceSurface = nullptr;

  /// Strategy to propagate to reference surface
  KalmanFitterTargetSurfaceStrategy referenceSurfaceStrategy =
      KalmanFitterTargetSurfaceStrategy::firstOrLast;

  /// Whether to consider multiple scattering
  bool multipleScattering = true;

  /// Whether to consider energy loss
  bool energyLoss = true;

  /// Whether to run filtering in reversed direction overwrite the
  /// ReverseFilteringLogic
  bool reversedFiltering = false;

  /// Factor by which the covariance of the input of the reversed filtering is
  /// scaled. This is only used in the backwardfiltering (if reversedFiltering
  /// is true or if the ReverseFilteringLogic return true for the track of
  /// interest)
  double reversedFilteringCovarianceScaling = 1.0;

  /// Whether to include non-linear correction during global to local
  /// transformation
  FreeToBoundCorrection freeToBoundCorrection;
};

template <typename traj_t>
struct KalmanFitterResult {
  /// Fitted states that the actor has handled.
  traj_t* fittedStates{nullptr};

  /// This is the index of the 'tip' of the track stored in multitrajectory.
  /// This corresponds to the last measurement state in the multitrajectory.
  /// Since this KF only stores one trajectory, it is unambiguous.
  /// Acts::MultiTrajectoryTraits::kInvalid is the start of a trajectory.
  std::size_t lastMeasurementIndex = Acts::MultiTrajectoryTraits::kInvalid;

  /// This is the index of the 'tip' of the states stored in multitrajectory.
  /// This corresponds to the last state in the multitrajectory.
  /// Since this KF only stores one trajectory, it is unambiguous.
  /// Acts::MultiTrajectoryTraits::kInvalid is the start of a trajectory.
  std::size_t lastTrackIndex = Acts::MultiTrajectoryTraits::kInvalid;

  /// The optional Parameters at the provided surface
  std::optional<BoundTrackParameters> fittedParameters;

  /// Counter for states with non-outlier measurements
  std::size_t measurementStates = 0;

  /// Counter for measurements holes
  /// A hole correspond to a surface with an associated detector element with no
  /// associated measurement. Holes are only taken into account if they are
  /// between the first and last measurements.
  std::size_t measurementHoles = 0;

  /// Counter for handled states
  std::size_t processedStates = 0;

  /// Indicator if smoothing has been done.
  bool smoothed = false;

  /// Indicator if navigation direction has been reversed
  bool reversed = false;

  /// Indicator if track fitting has been done
  bool finished = false;

  /// Measurement surfaces without hits
  std::vector<const Surface*> missedActiveSurfaces;

  /// Measurement surfaces handled in both forward and backward filtering
  std::vector<const Surface*> passedAgainSurfaces;

  Result<void> result{Result<void>::success()};
};

/// Kalman fitter implementation.
///
/// @tparam propagator_t Type of the propagation class
///
/// The Kalman filter contains an Actor and a Sequencer sub-class.
/// The Sequencer has to be part of the Navigator of the Propagator
/// in order to initialize and provide the measurement surfaces.
///
/// The Actor is part of the Propagation call and does the Kalman update
/// and eventually the smoothing.  Updater, Smoother and Calibrator are
/// given to the Actor for further use:
/// - The Updater is the implemented kalman updater formalism, it
///   runs via a visitor pattern through the measurements.
/// - The Smoother is called at the end of the filtering by the Actor.
///
/// Measurements are not required to be ordered for the KalmanFilter,
/// measurement ordering needs to be figured out by the navigation of
/// the propagator.
///
/// The void components are provided mainly for unit testing.
template <typename propagator_t, typename traj_t>
class KalmanFitter {
  /// The navigator type
  using KalmanNavigator = typename propagator_t::Navigator;

  /// The navigator has DirectNavigator type or not
  static constexpr bool isDirectNavigator =
      std::is_same_v<KalmanNavigator, DirectNavigator>;

 public:
  KalmanFitter(propagator_t pPropagator,
               std::unique_ptr<const Logger> _logger =
                   getDefaultLogger("KalmanFitter", Logging::INFO))
      : m_propagator(std::move(pPropagator)),
        m_logger{std::move(_logger)},
        m_actorLogger{m_logger->cloneWithSuffix("Actor")} {}

 private:
  /// The propagator for the transport and material update
  propagator_t m_propagator;

  /// The logger instance
  std::unique_ptr<const Logger> m_logger;
  std::unique_ptr<const Logger> m_actorLogger;

  const Logger& logger() const { return *m_logger; }

  /// @brief Propagator Actor plugin for the KalmanFilter
  ///
  /// @tparam parameters_t The type of parameters used for "local" parameters.
  /// @tparam calibrator_t The type of calibrator
  /// @tparam outlier_finder_t Type of the outlier finder class
  ///
  /// The KalmanActor does not rely on the measurements to be
  /// sorted along the track.
  template <typename parameters_t>
  class Actor {
   public:
    /// Broadcast the result_type
    using result_type = KalmanFitterResult<traj_t>;

    /// The target surface aboter
    SurfaceReached targetReached{std::numeric_limits<double>::lowest()};

    /// Strategy to propagate to target surface
    KalmanFitterTargetSurfaceStrategy targetSurfaceStrategy =
        KalmanFitterTargetSurfaceStrategy::firstOrLast;

    /// Allows retrieving measurements for a surface
    const std::map<GeometryIdentifier, SourceLink>* inputMeasurements = nullptr;

    /// Whether to consider multiple scattering.
    bool multipleScattering = true;

    /// Whether to consider energy loss.
    bool energyLoss = true;

    /// Whether run reversed filtering
    bool reversedFiltering = false;

    /// Scale the covariance before the reversed filtering
    double reversedFilteringCovarianceScaling = 1.0;

    /// Whether to include non-linear correction during global to local
    /// transformation
    FreeToBoundCorrection freeToBoundCorrection;

    /// Input MultiTrajectory
    std::shared_ptr<traj_t> outputStates;

    /// The logger instance
    const Logger* actorLogger{nullptr};

    /// Logger helper
    const Logger& logger() const { return *actorLogger; }

    KalmanFitterExtensions<traj_t> extensions;

    /// Calibration context for the fit
    const CalibrationContext* calibrationContext{nullptr};

    /// @brief Kalman actor operation
    ///
    /// @tparam propagator_state_t is the type of Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param navigator The navigator in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void act(propagator_state_t& state, const stepper_t& stepper,
             const navigator_t& navigator, result_type& result,
             const Logger& /*logger*/) const {
      assert(result.fittedStates && "No MultiTrajectory set");

      if (result.finished) {
        return;
      }

      ACTS_VERBOSE("KalmanFitter step at pos: "
                   << stepper.position(state.stepping).transpose()
                   << " dir: " << stepper.direction(state.stepping).transpose()
                   << " momentum: "
                   << stepper.absoluteMomentum(state.stepping));

      // Update:
      // - Waiting for a current surface
      auto surface = navigator.currentSurface(state.navigation);
      std::string direction = state.options.direction.toString();
      if (surface != nullptr) {
        // Check if the surface is in the measurement map
        // -> Get the measurement / calibrate
        // -> Create the predicted state
        // -> Check outlier behavior, if non-outlier:
        // -> Perform the kalman update
        // -> Fill track state information & update stepper information

        if (!result.smoothed && !result.reversed) {
          ACTS_VERBOSE("Perform " << direction << " filter step");
          auto res = filter(surface, state, stepper, navigator, result);
          if (!res.ok()) {
            ACTS_ERROR("Error in " << direction << " filter: " << res.error());
            result.result = res.error();
          }
        }
        if (result.reversed) {
          ACTS_VERBOSE("Perform " << direction << " filter step");
          auto res = reversedFilter(surface, state, stepper, navigator, result);
          if (!res.ok()) {
            ACTS_ERROR("Error in " << direction << " filter: " << res.error());
            result.result = res.error();
          }
        }
      }

      // Finalization:
      // when all track states have been handled or the navigation is breaked,
      // reset navigation&stepping before run reversed filtering or
      // proceed to run smoothing
      if (!result.smoothed && !result.reversed) {
        if (result.measurementStates == inputMeasurements->size() ||
            (result.measurementStates > 0 &&
             navigator.navigationBreak(state.navigation))) {
          // Remove the missing surfaces that occur after the last measurement
          result.missedActiveSurfaces.resize(result.measurementHoles);
          // now get track state proxy for the smoothing logic
          auto trackStateProxy =
              result.fittedStates->getTrackState(result.lastMeasurementIndex);
          if (reversedFiltering ||
              extensions.reverseFilteringLogic(trackStateProxy)) {
            // Start to run reversed filtering:
            // Reverse navigation direction and reset navigation and stepping
            // state to last measurement
            ACTS_VERBOSE("Reverse navigation direction.");
            auto res = reverse(state, stepper, navigator, result);
            if (!res.ok()) {
              ACTS_ERROR("Error in reversing navigation: " << res.error());
              result.result = res.error();
            }
          } else {
            // --> Search the starting state to run the smoothing
            // --> Call the smoothing
            // --> Set a stop condition when all track states have been
            // handled
            ACTS_VERBOSE("Finalize/run smoothing");
            auto res = finalize(state, stepper, navigator, result);
            if (!res.ok()) {
              ACTS_ERROR("Error in finalize: " << res.error());
              result.result = res.error();
            }
          }
        }
      }

      // Post-finalization:
      // - Progress to target/reference surface and built the final track
      // parameters
      if (result.smoothed || result.reversed) {
        if (result.smoothed) {
          // Update state and stepper with material effects
          // Only for smoothed as reverse filtering will handle this separately
          materialInteractor(navigator.currentSurface(state.navigation), state,
                             stepper, navigator,
                             MaterialUpdateStage::FullUpdate);
        }

        if (targetReached.surface == nullptr) {
          // If no target surface provided:
          // -> Return an error when using reversed filtering mode
          // -> Fitting is finished here
          if (result.reversed) {
            ACTS_ERROR(
                "The target surface needed for aborting reversed propagation "
                "is not provided");
            result.result =
                Result<void>(KalmanFitterError::ReversePropagationFailed);
          } else {
            ACTS_VERBOSE(
                "No target surface set. Completing without fitted track "
                "parameter");
            // Remember the track fitting is done
            result.finished = true;
          }
        } else if (targetReached.checkAbort(state, stepper, navigator,
                                            logger())) {
          ACTS_VERBOSE("Completing with fitted track parameter");
          // Transport & bind the parameter to the final surface
          auto res = stepper.boundState(state.stepping, *targetReached.surface,
                                        true, freeToBoundCorrection);
          if (!res.ok()) {
            ACTS_ERROR("Error in " << direction << " filter: " << res.error());
            result.result = res.error();
            return;
          }
          auto& fittedState = *res;
          // Assign the fitted parameters
          result.fittedParameters = std::get<BoundTrackParameters>(fittedState);

          // Reset smoothed status of states missed in reversed filtering
          if (result.reversed) {
            result.fittedStates->applyBackwards(
                result.lastMeasurementIndex, [&](auto trackState) {
                  auto fSurface = &trackState.referenceSurface();
                  if (trackState.typeFlags().test(
                          TrackStateFlag::MeasurementFlag) &&
                      !rangeContainsValue(result.passedAgainSurfaces,
                                          fSurface)) {
                    // If reversed filtering missed this surface, then there is
                    // no smoothed parameter
                    trackState.unset(TrackStatePropMask::Smoothed);
                    trackState.typeFlags().set(TrackStateFlag::OutlierFlag);
                  }
                });
          }
          // Remember the track fitting is done
          result.finished = true;
        }
      }
    }

    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    bool checkAbort(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/, const result_type& result,
                    const Logger& /*logger*/) const {
      return (!result.result.ok() || result.finished);
    }

    /// @brief Kalman actor operation: reverse direction
    ///
    /// @tparam propagator_state_t is the type of Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param navigator The navigator in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    Result<void> reverse(propagator_state_t& state, const stepper_t& stepper,
                         const navigator_t& navigator,
                         result_type& result) const {
      // Check if there is a measurement on track
      if (result.lastMeasurementIndex ==
          Acts::MultiTrajectoryTraits::kInvalid) {
        ACTS_ERROR("No point to reverse for a track without measurements.");
        return KalmanFitterError::ReversePropagationFailed;
      }

      // Remember the navigation direction has been reversed
      result.reversed = true;

      // Reverse navigation direction
      state.options.direction = state.options.direction.invert();

      // Reset propagator options
      // TODO Not sure if reset of pathLimit during propagation makes any sense
      state.options.pathLimit =
          state.options.direction * std::abs(state.options.pathLimit);

      // Get the last measurement state and reset navigation&stepping state
      // based on information on this state
      auto st = result.fittedStates->getTrackState(result.lastMeasurementIndex);

      // Update the stepping state
      stepper.resetState(
          state.stepping, st.filtered(),
          reversedFilteringCovarianceScaling * st.filteredCovariance(),
          st.referenceSurface(), state.options.stepping.maxStepSize);

      // For the last measurement state, smoothed is filtered
      st.smoothed() = st.filtered();
      st.smoothedCovariance() = st.filteredCovariance();
      result.passedAgainSurfaces.push_back(&st.referenceSurface());

      // Reset navigation state
      // We do not need to specify a target here since this will be handled
      // separately in the KF actor
      auto navigationOptions = state.navigation.options;
      navigationOptions.startSurface = &st.referenceSurface();
      navigationOptions.targetSurface = nullptr;
      state.navigation = navigator.makeState(navigationOptions);
      navigator.initialize(state.navigation, stepper.position(state.stepping),
                           stepper.direction(state.stepping),
                           state.options.direction);

      // Update material effects for last measurement state in reversed
      // direction
      materialInteractor(navigator.currentSurface(state.navigation), state,
                         stepper, navigator, MaterialUpdateStage::FullUpdate);

      return Result<void>::success();
    }

    /// @brief Kalman actor operation: update
    ///
    /// @tparam propagator_state_t is the type of Propagator state
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
      const bool precedingMeasurementExists = result.measurementStates > 0;
      const bool surfaceIsSensitive =
          surface->associatedDetectorElement() != nullptr;
      const bool surfaceHasMaterial = surface->surfaceMaterial() != nullptr;

      // Try to find the surface in the measurement surfaces
      auto sourceLinkIt = inputMeasurements->find(surface->geometryId());
      if (sourceLinkIt != inputMeasurements->end()) {
        // Screen output message
        ACTS_VERBOSE("Measurement surface " << surface->geometryId()
                                            << " detected.");
        // Transport the covariance to the surface
        stepper.transportCovarianceToBound(state.stepping, *surface,
                                           freeToBoundCorrection);

        // Update state and stepper with pre material effects
        materialInteractor(surface, state, stepper, navigator,
                           MaterialUpdateStage::PreUpdate);

        // do the kalman update (no need to perform covTransport here, hence no
        // point in performing globalToLocal correction)
        auto trackStateProxyRes = detail::kalmanHandleMeasurement(
            *calibrationContext, state, stepper, extensions, *surface,
            sourceLinkIt->second, *result.fittedStates, result.lastTrackIndex,
            false, logger());

        if (!trackStateProxyRes.ok()) {
          return trackStateProxyRes.error();
        }

        const auto& trackStateProxy = *trackStateProxyRes;
        result.lastTrackIndex = trackStateProxy.index();

        // Update the stepper if it is not an outlier
        if (trackStateProxy.typeFlags().test(
                Acts::TrackStateFlag::MeasurementFlag)) {
          // Update the stepping state with filtered parameters
          ACTS_VERBOSE("Filtering step successful, updated parameters are:\n"
                       << trackStateProxy.filtered().transpose());
          // update stepping state using filtered parameters after kalman
          stepper.update(state.stepping,
                         MultiTrajectoryHelpers::freeFiltered(
                             state.options.geoContext, trackStateProxy),
                         trackStateProxy.filtered(),
                         trackStateProxy.filteredCovariance(), *surface);
          // We count the state with measurement
          ++result.measurementStates;
        }

        // Update state and stepper with post material effects
        materialInteractor(surface, state, stepper, navigator,
                           MaterialUpdateStage::PostUpdate);
        // We count the processed state
        ++result.processedStates;
        // Update the number of holes count only when encountering a
        // measurement
        result.measurementHoles = result.missedActiveSurfaces.size();
        // Since we encountered a measurement update the lastMeasurementIndex to
        // the lastTrackIndex.
        result.lastMeasurementIndex = result.lastTrackIndex;

      } else if ((precedingMeasurementExists && surfaceIsSensitive) ||
                 surfaceHasMaterial) {
        // We only create track states here if there is already measurement
        // detected or if the surface has material (no holes before the first
        // measurement)
        auto trackStateProxyRes = detail::kalmanHandleNoMeasurement(
            state, stepper, *surface, *result.fittedStates,
            result.lastTrackIndex, true, logger(), precedingMeasurementExists,
            freeToBoundCorrection);

        if (!trackStateProxyRes.ok()) {
          return trackStateProxyRes.error();
        }

        const auto& trackStateProxy = *trackStateProxyRes;
        result.lastTrackIndex = trackStateProxy.index();

        if (trackStateProxy.typeFlags().test(TrackStateFlag::HoleFlag)) {
          // Count the missed surface
          result.missedActiveSurfaces.push_back(surface);
        }

        ++result.processedStates;

        // Update state and stepper with (possible) material effects
        materialInteractor(surface, state, stepper, navigator,
                           MaterialUpdateStage::FullUpdate);
      }
      return Result<void>::success();
    }

    /// @brief Kalman actor operation: update in reversed direction
    ///
    /// @tparam propagator_state_t is the type of Propagator state
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
    Result<void> reversedFilter(const Surface* surface,
                                propagator_state_t& state,
                                const stepper_t& stepper,
                                const navigator_t& navigator,
                                result_type& result) const {
      // Try to find the surface in the measurement surfaces
      auto sourceLinkIt = inputMeasurements->find(surface->geometryId());
      if (sourceLinkIt != inputMeasurements->end()) {
        // Screen output message
        ACTS_VERBOSE("Measurement surface "
                     << surface->geometryId()
                     << " detected in reversed propagation.");

        // No reversed filtering for last measurement state, but still update
        // with material effects
        if (result.reversed &&
            surface == navigator.startSurface(state.navigation)) {
          materialInteractor(surface, state, stepper, navigator,
                             MaterialUpdateStage::FullUpdate);
          return Result<void>::success();
        }

        // Transport the covariance to the surface
        stepper.transportCovarianceToBound(state.stepping, *surface,
                                           freeToBoundCorrection);

        // Update state and stepper with pre material effects
        materialInteractor(surface, state, stepper, navigator,
                           MaterialUpdateStage::PreUpdate);

        auto fittedStates = *result.fittedStates;

        // Add a <mask> TrackState entry multi trajectory. This allocates
        // storage for all components, which we will set later.
        TrackStatePropMask mask =
            TrackStatePropMask::Predicted | TrackStatePropMask::Filtered |
            TrackStatePropMask::Smoothed | TrackStatePropMask::Jacobian |
            TrackStatePropMask::Calibrated;
        const std::size_t currentTrackIndex = fittedStates.addTrackState(
            mask, Acts::MultiTrajectoryTraits::kInvalid);

        // now get track state proxy back
        typename traj_t::TrackStateProxy trackStateProxy =
            fittedStates.getTrackState(currentTrackIndex);

        // Set the trackStateProxy components with the state from the ongoing
        // propagation
        {
          trackStateProxy.setReferenceSurface(surface->getSharedPtr());
          // Bind the transported state to the current surface
          auto res = stepper.boundState(state.stepping, *surface, false,
                                        freeToBoundCorrection);
          if (!res.ok()) {
            return res.error();
          }
          const auto& [boundParams, jacobian, pathLength] = *res;

          // Fill the track state
          trackStateProxy.predicted() = boundParams.parameters();
          trackStateProxy.predictedCovariance() = state.stepping.cov;

          trackStateProxy.jacobian() = jacobian;
          trackStateProxy.pathLength() = pathLength;
        }

        // We have predicted parameters, so calibrate the uncalibrated input
        // measurement
        extensions.calibrator(state.geoContext, *calibrationContext,
                              sourceLinkIt->second, trackStateProxy);

        // If the update is successful, set covariance and
        auto updateRes =
            extensions.updater(state.geoContext, trackStateProxy, logger());
        if (!updateRes.ok()) {
          ACTS_ERROR("Backward update step failed: " << updateRes.error());
          return updateRes.error();
        } else {
          // Update the stepping state with filtered parameters
          ACTS_VERBOSE(
              "Backward filtering step successful, updated parameters are:\n"
              << trackStateProxy.filtered().transpose());

          // Fill the smoothed parameter for the existing track state
          result.fittedStates->applyBackwards(
              result.lastMeasurementIndex, [&](auto trackState) {
                auto fSurface = &trackState.referenceSurface();
                if (fSurface == surface) {
                  result.passedAgainSurfaces.push_back(surface);
                  trackState.smoothed() = trackStateProxy.filtered();
                  trackState.smoothedCovariance() =
                      trackStateProxy.filteredCovariance();
                  return false;
                }
                return true;
              });

          // update stepping state using filtered parameters after kalman
          // update We need to (re-)construct a BoundTrackParameters instance
          // here, which is a bit awkward.
          stepper.update(state.stepping,
                         MultiTrajectoryHelpers::freeFiltered(
                             state.options.geoContext, trackStateProxy),
                         trackStateProxy.filtered(),
                         trackStateProxy.filteredCovariance(), *surface);

          // Update state and stepper with post material effects
          materialInteractor(surface, state, stepper, navigator,
                             MaterialUpdateStage::PostUpdate);
        }
      } else if (surface->associatedDetectorElement() != nullptr ||
                 surface->surfaceMaterial() != nullptr) {
        // Transport covariance
        if (surface->associatedDetectorElement() != nullptr) {
          ACTS_VERBOSE("Detected hole on " << surface->geometryId()
                                           << " in reversed filtering");
          if (state.stepping.covTransport) {
            stepper.transportCovarianceToBound(state.stepping, *surface);
          }
        } else if (surface->surfaceMaterial() != nullptr) {
          ACTS_VERBOSE("Detected in-sensitive surface "
                       << surface->geometryId() << " in reversed filtering");
          if (state.stepping.covTransport) {
            stepper.transportCovarianceToCurvilinear(state.stepping);
          }
        }
        // Not creating bound state here, so need manually reinitialize
        // jacobian
        stepper.setIdentityJacobian(state.stepping);
        if (surface->surfaceMaterial() != nullptr) {
          // Update state and stepper with material effects
          materialInteractor(surface, state, stepper, navigator,
                             MaterialUpdateStage::FullUpdate);
        }
      }

      return Result<void>::success();
    }

    /// @brief Kalman actor operation: material interaction
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
    ///
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void materialInteractor(const Surface* surface, propagator_state_t& state,
                            const stepper_t& stepper,
                            const navigator_t& navigator,
                            const MaterialUpdateStage& updateStage) const {
      // Protect against null surface
      if (!surface) {
        ACTS_VERBOSE(
            "Surface is nullptr. Cannot be used for material interaction");
        return;
      }

      // Indicator if having material
      bool hasMaterial = false;

      if (surface && surface->surfaceMaterial()) {
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
                       << interaction.Eloss << ", "
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

    /// @brief Kalman actor operation: finalize
    ///
    /// @tparam propagator_state_t is the type of Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param navigator The navigator in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    Result<void> finalize(propagator_state_t& state, const stepper_t& stepper,
                          const navigator_t& navigator,
                          result_type& result) const {
      // Remember you smoothed the track states
      result.smoothed = true;

      // Get the indices of the first states (can be either a measurement or
      // material);
      std::size_t firstStateIndex = result.lastMeasurementIndex;
      // Count track states to be smoothed
      std::size_t nStates = 0;
      result.fittedStates->applyBackwards(
          result.lastMeasurementIndex, [&](auto st) {
            bool isMeasurement =
                st.typeFlags().test(TrackStateFlag::MeasurementFlag);
            bool isOutlier = st.typeFlags().test(TrackStateFlag::OutlierFlag);
            // We are excluding non measurement states and outlier here. Those
            // can decrease resolution because only the smoothing corrected the
            // very first prediction as filtering is not possible.
            if (isMeasurement && !isOutlier) {
              firstStateIndex = st.index();
            }
            nStates++;
          });
      // Return error if the track has no measurement states (but this should
      // not happen)
      if (nStates == 0) {
        ACTS_ERROR("Smoothing for a track without measurements.");
        return KalmanFitterError::SmoothFailed;
      }
      // Screen output for debugging
      ACTS_VERBOSE("Apply smoothing on " << nStates
                                         << " filtered track states.");

      // Smooth the track states
      auto smoothRes =
          extensions.smoother(state.geoContext, *result.fittedStates,
                              result.lastMeasurementIndex, logger());
      if (!smoothRes.ok()) {
        ACTS_ERROR("Smoothing step failed: " << smoothRes.error());
        return smoothRes.error();
      }

      // Return in case no target surface
      if (targetReached.surface == nullptr) {
        return Result<void>::success();
      }

      // Obtain the smoothed parameters at first/last measurement state
      auto firstCreatedState =
          result.fittedStates->getTrackState(firstStateIndex);
      auto lastCreatedMeasurement =
          result.fittedStates->getTrackState(result.lastMeasurementIndex);

      // Lambda to get the intersection of the free params on the target surface
      auto target = [&](const FreeVector& freeVector) -> SurfaceIntersection {
        return targetReached.surface
            ->intersect(
                state.geoContext, freeVector.segment<3>(eFreePos0),
                state.options.direction * freeVector.segment<3>(eFreeDir0),
                BoundaryTolerance::None(), state.options.surfaceTolerance)
            .closest();
      };

      // The smoothed free params at the first/last measurement state.
      // (the first state can also be a material state)
      auto firstParams = MultiTrajectoryHelpers::freeSmoothed(
          state.options.geoContext, firstCreatedState);
      auto lastParams = MultiTrajectoryHelpers::freeSmoothed(
          state.options.geoContext, lastCreatedMeasurement);
      // Get the intersections of the smoothed free parameters with the target
      // surface
      const auto firstIntersection = target(firstParams);
      const auto lastIntersection = target(lastParams);

      // Update the stepping parameters - in order to progress to destination.
      // At the same time, reverse navigation direction for further stepping if
      // necessary.
      // @note The stepping parameters is updated to the smoothed parameters at
      // either the first measurement state or the last measurement state. It
      // assumes the target surface is not within the first and the last
      // smoothed measurement state. Also, whether the intersection is on
      // surface is not checked here.
      bool useFirstTrackState = true;
      switch (targetSurfaceStrategy) {
        case KalmanFitterTargetSurfaceStrategy::first:
          useFirstTrackState = true;
          break;
        case KalmanFitterTargetSurfaceStrategy::last:
          useFirstTrackState = false;
          break;
        case KalmanFitterTargetSurfaceStrategy::firstOrLast:
          useFirstTrackState = std::abs(firstIntersection.pathLength()) <=
                               std::abs(lastIntersection.pathLength());
          break;
        default:
          ACTS_ERROR("Unknown target surface strategy");
          return KalmanFitterError::SmoothFailed;
      }
      bool reverseDirection = false;
      if (useFirstTrackState) {
        stepper.resetState(state.stepping, firstCreatedState.smoothed(),
                           firstCreatedState.smoothedCovariance(),
                           firstCreatedState.referenceSurface(),
                           state.options.stepping.maxStepSize);
        reverseDirection = firstIntersection.pathLength() < 0;
      } else {
        stepper.resetState(state.stepping, lastCreatedMeasurement.smoothed(),
                           lastCreatedMeasurement.smoothedCovariance(),
                           lastCreatedMeasurement.referenceSurface(),
                           state.options.stepping.maxStepSize);
        reverseDirection = lastIntersection.pathLength() < 0;
      }
      // Reverse the navigation direction if necessary
      if (reverseDirection) {
        ACTS_VERBOSE(
            "Reverse navigation direction after smoothing for reaching the "
            "target surface");
        state.options.direction = state.options.direction.invert();
      }
      const auto& surface = useFirstTrackState
                                ? firstCreatedState.referenceSurface()
                                : lastCreatedMeasurement.referenceSurface();

      ACTS_VERBOSE(
          "Smoothing successful, updating stepping state to smoothed "
          "parameters at surface "
          << surface.geometryId() << ". Prepared to reach the target surface.");

      // Reset the navigation state to enable propagation towards the target
      // surface
      // Set targetSurface to nullptr as it is handled manually in the actor
      auto navigationOptions = state.navigation.options;
      navigationOptions.startSurface = &surface;
      navigationOptions.targetSurface = nullptr;
      state.navigation = navigator.makeState(navigationOptions);
      navigator.initialize(state.navigation, stepper.position(state.stepping),
                           stepper.direction(state.stepping),
                           state.options.direction);

      return Result<void>::success();
    }
  };

 public:
  /// Fit implementation of the forward filter, calls the
  /// the filter and smoother/reversed filter
  ///
  /// @tparam source_link_iterator_t Iterator type used to pass source links
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam parameters_t Type of parameters used for local parameters
  /// @tparam track_container_t Type of the track container
  ///
  /// @param it Begin iterator for the fittable uncalibrated measurements
  /// @param end End iterator for the fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param kfOptions KalmanOptions steering the fit
  /// @param trackContainer Input track container storage to append into
  /// @note The input measurements are given in the form of @c SourceLink s.
  /// It's the calibrators job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename source_link_iterator_t, typename start_parameters_t,
            typename parameters_t = BoundTrackParameters,
            TrackContainerFrontend track_container_t>
  Result<typename track_container_t::TrackProxy> fit(
      source_link_iterator_t it, source_link_iterator_t end,
      const start_parameters_t& sParameters,
      const KalmanFitterOptions<traj_t>& kfOptions,
      track_container_t& trackContainer) const
    requires(!isDirectNavigator)
  {
    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyway, so the map can own them.
    ACTS_VERBOSE("Preparing " << std::distance(it, end)
                              << " input measurements");
    std::map<GeometryIdentifier, SourceLink> inputMeasurements;
    // for (const auto& sl : sourceLinks) {
    for (; it != end; ++it) {
      SourceLink sl = *it;
      const Surface* surface = kfOptions.extensions.surfaceAccessor(sl);
      // @TODO: This can probably change over to surface pointers as keys
      auto geoId = surface->geometryId();
      inputMeasurements.emplace(geoId, std::move(sl));
    }

    // Create the ActorList
    using KalmanActor = Actor<parameters_t>;

    using KalmanResult = typename KalmanActor::result_type;
    using Actors = ActorList<KalmanActor>;
    using PropagatorOptions = typename propagator_t::template Options<Actors>;

    // Create relevant options for the propagation options
    PropagatorOptions propagatorOptions(kfOptions.geoContext,
                                        kfOptions.magFieldContext);

    // Set the trivial propagator options
    propagatorOptions.setPlainOptions(kfOptions.propagatorPlainOptions);

    // Add the measurement surface as external surface to navigator.
    // We will try to hit those surface by ignoring boundary checks.
    for (const auto& [surfaceId, _] : inputMeasurements) {
      propagatorOptions.navigation.insertExternalSurface(surfaceId);
    }

    // Catch the actor and set the measurements
    auto& kalmanActor = propagatorOptions.actorList.template get<KalmanActor>();
    kalmanActor.inputMeasurements = &inputMeasurements;
    kalmanActor.targetReached.surface = kfOptions.referenceSurface;
    kalmanActor.targetSurfaceStrategy = kfOptions.referenceSurfaceStrategy;
    kalmanActor.multipleScattering = kfOptions.multipleScattering;
    kalmanActor.energyLoss = kfOptions.energyLoss;
    kalmanActor.reversedFiltering = kfOptions.reversedFiltering;
    kalmanActor.reversedFilteringCovarianceScaling =
        kfOptions.reversedFilteringCovarianceScaling;
    kalmanActor.freeToBoundCorrection = kfOptions.freeToBoundCorrection;
    kalmanActor.calibrationContext = &kfOptions.calibrationContext.get();
    kalmanActor.extensions = kfOptions.extensions;
    kalmanActor.actorLogger = m_actorLogger.get();

    return fit_impl<start_parameters_t, PropagatorOptions, KalmanResult,
                    track_container_t>(sParameters, propagatorOptions,
                                       trackContainer);
  }

  /// Fit implementation of the forward filter, calls the
  /// the filter and smoother/reversed filter
  ///
  /// @tparam source_link_iterator_t Iterator type used to pass source links
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam parameters_t Type of parameters used for local parameters
  /// @tparam track_container_t Type of the track container
  ///
  /// @param it Begin iterator for the fittable uncalibrated measurements
  /// @param end End iterator for the fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param kfOptions KalmanOptions steering the fit
  /// @param sSequence surface sequence used to initialize a DirectNavigator
  /// @param trackContainer Input track container storage to append into
  /// @note The input measurements are given in the form of @c SourceLinks.
  /// It's
  /// @c calibrator_t's job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename source_link_iterator_t, typename start_parameters_t,
            typename parameters_t = BoundTrackParameters,
            TrackContainerFrontend track_container_t>
  Result<typename track_container_t::TrackProxy> fit(
      source_link_iterator_t it, source_link_iterator_t end,
      const start_parameters_t& sParameters,
      const KalmanFitterOptions<traj_t>& kfOptions,
      const std::vector<const Surface*>& sSequence,
      track_container_t& trackContainer) const
    requires(isDirectNavigator)
  {
    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyway, so the map can own them.
    ACTS_VERBOSE("Preparing " << std::distance(it, end)
                              << " input measurements");
    std::map<GeometryIdentifier, SourceLink> inputMeasurements;
    for (; it != end; ++it) {
      SourceLink sl = *it;
      const Surface* surface = kfOptions.extensions.surfaceAccessor(sl);
      // @TODO: This can probably change over to surface pointers as keys
      auto geoId = surface->geometryId();
      inputMeasurements.emplace(geoId, std::move(sl));
    }

    // Create the ActorList
    using KalmanActor = Actor<parameters_t>;

    using KalmanResult = typename KalmanActor::result_type;
    using Actors = ActorList<KalmanActor>;
    using PropagatorOptions = typename propagator_t::template Options<Actors>;

    // Create relevant options for the propagation options
    PropagatorOptions propagatorOptions(kfOptions.geoContext,
                                        kfOptions.magFieldContext);

    // Set the trivial propagator options
    propagatorOptions.setPlainOptions(kfOptions.propagatorPlainOptions);

    // Catch the actor and set the measurements
    auto& kalmanActor = propagatorOptions.actorList.template get<KalmanActor>();
    kalmanActor.inputMeasurements = &inputMeasurements;
    kalmanActor.targetReached.surface = kfOptions.referenceSurface;
    kalmanActor.targetSurfaceStrategy = kfOptions.referenceSurfaceStrategy;
    kalmanActor.multipleScattering = kfOptions.multipleScattering;
    kalmanActor.energyLoss = kfOptions.energyLoss;
    kalmanActor.reversedFiltering = kfOptions.reversedFiltering;
    kalmanActor.reversedFilteringCovarianceScaling =
        kfOptions.reversedFilteringCovarianceScaling;
    kalmanActor.extensions = kfOptions.extensions;
    kalmanActor.actorLogger = m_actorLogger.get();

    // Set the surface sequence
    propagatorOptions.navigation.surfaces = sSequence;

    return fit_impl<start_parameters_t, PropagatorOptions, KalmanResult,
                    track_container_t>(sParameters, propagatorOptions,
                                       trackContainer);
  }

 private:
  /// Common fit implementation
  ///
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam actor_list_t Type of the actor list
  /// @tparam aborter_list_t Type of the abort list
  /// @tparam kalman_result_t Type of the KF result
  /// @tparam track_container_t Type of the track container
  ///
  /// @param sParameters The initial track parameters
  /// @param propagatorOptions The Propagator Options
  /// @param trackContainer Input track container storage to append into
  ///
  /// @return the output as an output track
  template <typename start_parameters_t, typename propagator_options_t,
            typename kalman_result_t, TrackContainerFrontend track_container_t>
  auto fit_impl(const start_parameters_t& sParameters,
                const propagator_options_t& propagatorOptions,
                track_container_t& trackContainer) const
      -> Result<typename track_container_t::TrackProxy> {
    auto propagatorState =
        m_propagator.makeState(sParameters, propagatorOptions);

    auto& kalmanResult =
        propagatorState.template get<KalmanFitterResult<traj_t>>();
    kalmanResult.fittedStates = &trackContainer.trackStateContainer();

    // Run the fitter
    auto result = m_propagator.propagate(propagatorState);

    if (!result.ok()) {
      ACTS_ERROR("Propagation failed: " << result.error());
      return result.error();
    }

    /// It could happen that the fit ends in zero measurement states.
    /// The result gets meaningless so such case is regarded as fit failure.
    if (kalmanResult.result.ok() && !kalmanResult.measurementStates) {
      kalmanResult.result = Result<void>(KalmanFitterError::NoMeasurementFound);
    }

    if (!kalmanResult.result.ok()) {
      ACTS_ERROR("KalmanFilter failed: "
                 << kalmanResult.result.error() << ", "
                 << kalmanResult.result.error().message());
      return kalmanResult.result.error();
    }

    auto track = trackContainer.makeTrack();
    track.tipIndex() = kalmanResult.lastMeasurementIndex;
    if (kalmanResult.fittedParameters) {
      const auto& params = kalmanResult.fittedParameters.value();
      track.parameters() = params.parameters();
      track.covariance() = params.covariance().value();
      track.setReferenceSurface(params.referenceSurface().getSharedPtr());
    }

    calculateTrackQuantities(track);

    if (trackContainer.hasColumn(hashString("smoothed"))) {
      track.template component<bool, hashString("smoothed")>() =
          kalmanResult.smoothed;
    }

    if (trackContainer.hasColumn(hashString("reversed"))) {
      track.template component<bool, hashString("reversed")>() =
          kalmanResult.reversed;
    }

    // Return the converted Track
    return track;
  }
};

}  // namespace Acts
