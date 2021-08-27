// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/TrackFitting/detail/VoidKalmanComponents.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <functional>
#include <map>
#include <memory>

namespace Acts {

/// Combined options for the Kalman fitter.
///
/// @tparam calibrator_t Source link type, should be semiregular.
/// @tparam outlier_finder_t Outlier finder type, shoule be semiregular.
template <typename calibrator_t, typename outlier_finder_t>
struct KalmanFitterOptions {
  using Calibrator = calibrator_t;
  using OutlierFinder = outlier_finder_t;

  /// PropagatorOptions with context.
  ///
  /// @param gctx The goemetry context for this fit
  /// @param mctx The magnetic context for this fit
  /// @param cctx The calibration context for this fit
  /// @param calibrator_ The source link calibrator
  /// @param outlierFinder_ The outlier finder
  /// @param logger_ The logger wrapper
  /// @param pOptions The plain propagator options
  /// @param rSurface The reference surface for the fit to be expressed at
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  /// @param rFiltering Whether to run filtering in reversed direction as
  /// smoothing
  KalmanFitterOptions(const GeometryContext& gctx,
                      const MagneticFieldContext& mctx,
                      std::reference_wrapper<const CalibrationContext> cctx,
                      Calibrator calibrator_, OutlierFinder outlierFinder_,
                      LoggerWrapper logger_,
                      const PropagatorPlainOptions& pOptions,
                      const Surface* rSurface = nullptr,
                      bool mScattering = true, bool eLoss = true,
                      bool rFiltering = false)
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        calibrator(std::move(calibrator_)),
        outlierFinder(std::move(outlierFinder_)),
        propagatorPlainOptions(pOptions),
        referenceSurface(rSurface),
        multipleScattering(mScattering),
        energyLoss(eLoss),
        reversedFiltering(rFiltering),
        logger(logger_) {}
  /// Contexts are required and the options must not be default-constructible.
  KalmanFitterOptions() = delete;

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  /// The source link calibrator.
  Calibrator calibrator;

  /// The outlier finder.
  OutlierFinder outlierFinder;

  /// The trivial propagator options
  PropagatorPlainOptions propagatorPlainOptions;

  /// The reference Surface
  const Surface* referenceSurface = nullptr;

  /// Whether to consider multiple scattering
  bool multipleScattering = true;

  /// Whether to consider energy loss
  bool energyLoss = true;

  /// Whether to run filtering in reversed direction
  bool reversedFiltering = false;

  /// Logger
  LoggerWrapper logger;
};

template <typename source_link_t>
struct KalmanFitterResult {
  // Fitted states that the actor has handled.
  MultiTrajectory<source_link_t> fittedStates;

  // This is the index of the 'tip' of the track stored in multitrajectory.
  // This correspond to the last measurment state in the multitrajectory.
  // Since this KF only stores one trajectory, it is unambiguous.
  // SIZE_MAX is the start of a trajectory.
  size_t lastMeasurementIndex = SIZE_MAX;

  // This is the index of the 'tip' of the states stored in multitrajectory.
  // This correspond to the last state in the multitrajectory.
  // Since this KF only stores one trajectory, it is
  // unambiguous. SIZE_MAX is the start of a trajectory.
  size_t lastTrackIndex = SIZE_MAX;

  // The optional Parameters at the provided surface
  std::optional<BoundTrackParameters> fittedParameters;

  // Counter for states with non-outlier measurements
  size_t measurementStates = 0;

  // Counter for measurements holes
  // A hole correspond to a surface with an associated detector element with no
  // associated measurment. Holes are only taken into account if they are
  // between the first and last measurements.
  size_t measurementHoles = 0;

  // Counter for handled states
  size_t processedStates = 0;

  // Indicator if smoothing has been done.
  bool smoothed = false;

  // Indicator if navigation direction has been reversed
  bool reversed = false;

  // Indicator if track fitting has been done
  bool finished = false;

  // Measurement surfaces without hits
  std::vector<const Surface*> missedActiveSurfaces;

  // Measurement surfaces handled in both forward and backward filtering
  std::vector<const Surface*> passedAgainSurfaces;

  Result<void> result{Result<void>::success()};
};

/// Kalman fitter implementation.
///
/// @tparam propagator_t Type of the propagation class
/// @tparam updater_t Type of the kalman updater class
/// @tparam smoother_t Type of the kalman smoother class
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
template <typename propagator_t, typename updater_t = VoidKalmanUpdater,
          typename smoother_t = VoidKalmanSmoother>
class KalmanFitter {
  /// The navigator type
  using KalmanNavigator = typename propagator_t::Navigator;

  /// The navigator has DirectNavigator type or not
  static constexpr bool isDirectNavigator =
      std::is_same<KalmanNavigator, DirectNavigator>::value;

 public:
  KalmanFitter(propagator_t pPropagator)
      : m_propagator(std::move(pPropagator)) {}

 private:
  /// The propgator for the transport and material update
  propagator_t m_propagator;

  /// @brief Propagator Actor plugin for the KalmanFilter
  ///
  /// @tparam source_link_t is an type fulfilling the @c SourceLinkConcept
  /// @tparam parameters_t The type of parameters used for "local" paremeters.
  /// @tparam calibrator_t The type of calibrator
  /// @tparam outlier_finder_t Type of the outlier finder class
  ///
  /// The KalmanActor does not rely on the measurements to be
  /// sorted along the track.
  template <typename source_link_t, typename parameters_t,
            typename calibrator_t, typename outlier_finder_t>
  class Actor {
   public:
    /// Broadcast the result_type
    using result_type = KalmanFitterResult<source_link_t>;

    /// The target surface
    const Surface* targetSurface = nullptr;

    /// Allows retrieving measurements for a surface
    const std::map<GeometryIdentifier, source_link_t>* inputMeasurements =
        nullptr;

    /// Whether to consider multiple scattering.
    bool multipleScattering = true;

    /// Whether to consider energy loss.
    bool energyLoss = true;

    /// Whether run reversed filtering
    bool reversedFiltering = false;

    /// @brief Kalman actor operation
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    void operator()(propagator_state_t& state, const stepper_t& stepper,
                    result_type& result) const {
      const auto& logger = state.options.logger;

      if (result.finished) {
        return;
      }

      ACTS_VERBOSE("KalmanFitter step");

      // Add the measurement surface as external surface to navigator.
      // We will try to hit those surface by ignoring boundary checks.
      if constexpr (not isDirectNavigator) {
        if (result.processedStates == 0) {
          for (auto measurementIt = inputMeasurements->begin();
               measurementIt != inputMeasurements->end(); measurementIt++) {
            state.navigation.externalSurfaces.insert(
                std::pair<uint64_t, GeometryIdentifier>(
                    measurementIt->first.layer(), measurementIt->first));
          }
        }
      }

      // Update:
      // - Waiting for a current surface
      auto surface = state.navigation.currentSurface;
      std::string direction =
          (state.stepping.navDir == forward) ? "forward" : "backward";
      if (surface != nullptr) {
        // Check if the surface is in the measurement map
        // -> Get the measurement / calibrate
        // -> Create the predicted state
        // -> Check outlier behavior, if non-outlier:
        // -> Perform the kalman update
        // -> Fill strack state information & update stepper information

        if (not result.smoothed and not result.reversed) {
          ACTS_VERBOSE("Perform " << direction << " filter step");
          auto res = filter(surface, state, stepper, result);
          if (!res.ok()) {
            ACTS_ERROR("Error in " << direction << " filter: " << res.error());
            result.result = res.error();
          }
        }
        if (result.reversed) {
          ACTS_VERBOSE("Perform " << direction << " filter step");
          auto res = reversedFilter(surface, state, stepper, result);
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
      if (not result.smoothed and not result.reversed) {
        if (result.measurementStates == inputMeasurements->size() or
            (result.measurementStates > 0 and
             state.navigation.navigationBreak)) {
          // Remove the missing surfaces that occur after the last measurement
          result.missedActiveSurfaces.resize(result.measurementHoles);
          if (reversedFiltering) {
            // Start to run reversed filtering:
            // Reverse navigation direction and reset navigation and stepping
            // state to last measurement
            ACTS_VERBOSE("Reverse navigation direction.");
            auto res = reverse(state, stepper, result);
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
            auto res = finalize(state, stepper, result);
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
      if (result.smoothed or result.reversed) {
        if (targetSurface == nullptr) {
          // If no target surface provided:
          // -> Return an error when using reversed filtering mode
          // -> Fitting is finished here
          if (reversedFiltering) {
            ACTS_ERROR(
                "The target surface needed for aborting reversed propagation "
                "is not provided");
            result.result =
                Result<void>(KalmanFitterError::BackwardUpdateFailed);
          } else {
            ACTS_VERBOSE(
                "No target surface set. Completing without fitted track "
                "parameter");
            // Remember the track fitting is done
            result.finished = true;
          }
        } else if (targetReached(state, stepper, *targetSurface)) {
          ACTS_VERBOSE("Completing with fitted track parameter");
          // Transport & bind the parameter to the final surface
          auto res = stepper.boundState(state.stepping, *targetSurface);
          if (!res.ok()) {
            ACTS_ERROR("Error in " << direction << " filter: " << res.error());
            result.result = res.error();
            return;
          }
          auto& fittedState = *res;
          // Assign the fitted parameters
          result.fittedParameters = std::get<BoundTrackParameters>(fittedState);

          // Reset smoothed status of states missed in reversed filtering
          if (reversedFiltering) {
            result.fittedStates.applyBackwards(
                result.lastMeasurementIndex, [&](auto trackState) {
                  auto fSurface = &trackState.referenceSurface();
                  auto surface_it = std::find_if(
                      result.passedAgainSurfaces.begin(),
                      result.passedAgainSurfaces.end(),
                      [=](const Surface* s) { return s == fSurface; });
                  if (surface_it == result.passedAgainSurfaces.end()) {
                    // If reversed filtering missed this surface, then there is
                    // no smoothed parameter
                    trackState.data().ismoothed =
                        detail_lt::IndexData::kInvalid;
                  }
                });
          }
          // Remember the track fitting is done
          result.finished = true;
        }
      }
    }

    /// @brief Kalman actor operation : reverse direction
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state objecte
    template <typename propagator_state_t, typename stepper_t>
    Result<void> reverse(propagator_state_t& state, stepper_t& stepper,
                         result_type& result) const {
      const auto& logger = state.options.logger;

      // Check if there is a measurement on track
      if (result.lastMeasurementIndex == SIZE_MAX) {
        ACTS_ERROR("No point to reverse for a track without measurements.");
        return KalmanFitterError::ReverseNavigationFailed;
      }

      // Remember the navigation direciton has been reversed
      result.reversed = true;

      // Reverse navigation direction
      state.stepping.navDir =
          (state.stepping.navDir == forward) ? backward : forward;

      // Reset propagator options
      state.options.maxStepSize =
          state.stepping.navDir * std::abs(state.options.maxStepSize);
      // Not sure if reset of pathLimit during propagation makes any sense
      state.options.pathLimit =
          state.stepping.navDir * std::abs(state.options.pathLimit);

      // Get the last measurement state and reset navigation&stepping state
      // based on information on this state
      auto st = result.fittedStates.getTrackState(result.lastMeasurementIndex);

      // Update the stepping state
      stepper.resetState(state.stepping, st.filtered(), st.filteredCovariance(),
                         st.referenceSurface(), state.stepping.navDir,
                         state.options.maxStepSize);

      // For the last measurement state, smoothed is filtered
      st.smoothed() = st.filtered();
      st.smoothedCovariance() = st.filteredCovariance();
      result.passedAgainSurfaces.push_back(&st.referenceSurface());

      // Reset navigation state
      state.navigation.reset(state.geoContext, stepper.position(state.stepping),
                             stepper.direction(state.stepping),
                             state.stepping.navDir, &st.referenceSurface(),
                             targetSurface);

      // Update material effects for last measurement state in reversed
      // direction
      materialInteractor(state.navigation.currentSurface, state, stepper);

      return Result<void>::success();
    }

    /// @brief Kalman actor operation : update
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param surface The surface where the update happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result The mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    Result<void> filter(const Surface* surface, propagator_state_t& state,
                        const stepper_t& stepper, result_type& result) const {
      const auto& logger = state.options.logger;
      // Try to find the surface in the measurement surfaces
      auto sourcelink_it = inputMeasurements->find(surface->geometryId());
      if (sourcelink_it != inputMeasurements->end()) {
        // Screen output message
        ACTS_VERBOSE("Measurement surface " << surface->geometryId()
                                            << " detected.");
        // Transport the covariance to the surface
        stepper.transportCovarianceToBound(state.stepping, *surface);

        // Update state and stepper with pre material effects
        materialInteractor(surface, state, stepper, preUpdate);

        // Bind the transported state to the current surface
        auto res = stepper.boundState(state.stepping, *surface, false);
        if (!res.ok()) {
          return res.error();
        }
        auto& [boundParams, jacobian, pathLength] = *res;

        // add a full TrackState entry multi trajectory
        // (this allocates storage for all components, we will set them later)
        result.lastTrackIndex = result.fittedStates.addTrackState(
            TrackStatePropMask::All, result.lastTrackIndex);

        // now get track state proxy back
        auto trackStateProxy =
            result.fittedStates.getTrackState(result.lastTrackIndex);

        trackStateProxy.setReferenceSurface(surface->getSharedPtr());

        // assign the source link to the track state
        trackStateProxy.uncalibrated() = sourcelink_it->second;

        // Fill the track state
        trackStateProxy.predicted() = std::move(boundParams.parameters());
        if (boundParams.covariance().has_value()) {
          trackStateProxy.predictedCovariance() =
              std::move(*boundParams.covariance());
        }
        trackStateProxy.jacobian() = std::move(jacobian);
        trackStateProxy.pathLength() = std::move(pathLength);

        // We have predicted parameters, so calibrate the uncalibrated input
        // measuerement
        std::visit(
            [&](const auto& calibrated) {
              trackStateProxy.setCalibrated(calibrated);
            },
            m_calibrator(trackStateProxy.uncalibrated(),
                         trackStateProxy.predicted()));

        // Get and set the type flags
        auto& typeFlags = trackStateProxy.typeFlags();
        typeFlags.set(TrackStateFlag::ParameterFlag);
        if (surface->surfaceMaterial() != nullptr) {
          typeFlags.set(TrackStateFlag::MaterialFlag);
        }

        // Check if the state is an outlier.
        // If not, run Kalman update, tag it as a
        // measurement and update the stepping state. Otherwise, just tag it as
        // an outlier
        if (not m_outlierFinder(trackStateProxy)) {
          // Run Kalman update
          auto updateRes = m_updater(state.geoContext, trackStateProxy,
                                     state.stepping.navDir, logger);
          if (!updateRes.ok()) {
            ACTS_ERROR("Update step failed: " << updateRes.error());
            return updateRes.error();
          }
          // Set the measurement type flag
          typeFlags.set(TrackStateFlag::MeasurementFlag);
          // Update the stepping state with filtered parameters
          ACTS_VERBOSE("Filtering step successful, updated parameters are : \n"
                       << trackStateProxy.filtered().transpose());
          // update stepping state using filtered parameters after kalman
          stepper.update(state.stepping,
                         MultiTrajectoryHelpers::freeFiltered(
                             state.options.geoContext, trackStateProxy),
                         trackStateProxy.filtered(),
                         trackStateProxy.filteredCovariance(), *surface);
          // We count the state with measurement
          ++result.measurementStates;
        } else {
          ACTS_VERBOSE(
              "Filtering step successful. But measurement is deterimined "
              "to "
              "be an outlier. Stepping state is not updated.")
          // Set the outlier type flag
          typeFlags.set(TrackStateFlag::OutlierFlag);
          trackStateProxy.data().ifiltered = trackStateProxy.data().ipredicted;
        }

        // Update state and stepper with post material effects
        materialInteractor(surface, state, stepper, postUpdate);
        // We count the processed state
        ++result.processedStates;
        // Update the number of holes count only when encoutering a
        // measurement
        result.measurementHoles = result.missedActiveSurfaces.size();
        // Since we encountered a measurment update the lastMeasurementIndex to
        // the lastTrackIndex.
        result.lastMeasurementIndex = result.lastTrackIndex;

      } else if (surface->associatedDetectorElement() != nullptr ||
                 surface->surfaceMaterial() != nullptr) {
        // We only create track states here if there is already measurement
        // detected or if the surface has material (no holes before the first
        // measurement)
        if (result.measurementStates > 0 ||
            surface->surfaceMaterial() != nullptr) {
          // No source links on surface, add either hole or passive material
          // TrackState entry multi trajectory. No storage allocation for
          // uncalibrated/calibrated measurement and filtered parameter
          result.lastTrackIndex = result.fittedStates.addTrackState(
              ~(TrackStatePropMask::Uncalibrated |
                TrackStatePropMask::Calibrated | TrackStatePropMask::Filtered),
              result.lastTrackIndex);

          // now get track state proxy back
          auto trackStateProxy =
              result.fittedStates.getTrackState(result.lastTrackIndex);

          // Set the surface
          trackStateProxy.setReferenceSurface(surface->getSharedPtr());

          // Set the track state flags
          auto& typeFlags = trackStateProxy.typeFlags();
          typeFlags.set(TrackStateFlag::ParameterFlag);
          if (surface->surfaceMaterial() != nullptr) {
            typeFlags.set(TrackStateFlag::MaterialFlag);
          }
          if (surface->associatedDetectorElement() != nullptr) {
            ACTS_VERBOSE("Detected hole on " << surface->geometryId());
            // If the surface is sensitive, set the hole type flag
            typeFlags.set(TrackStateFlag::HoleFlag);

            // Count the missed surface
            result.missedActiveSurfaces.push_back(surface);
          } else if (surface->surfaceMaterial() != nullptr) {
            ACTS_VERBOSE("Detected in-sensitive surface "
                         << surface->geometryId());
          }

          // Transport & bind the state to the current surface
          auto res = stepper.boundState(state.stepping, *surface);
          if (!res.ok()) {
            ACTS_ERROR("Propagate to surface " << surface->geometryId()
                                               << " failed: " << res.error());
            return res.error();
          }
          auto& [boundParams, jacobian, pathLength] = *res;

          // Fill the track state
          trackStateProxy.predicted() = std::move(boundParams.parameters());
          if (boundParams.covariance().has_value()) {
            trackStateProxy.predictedCovariance() =
                std::move(*boundParams.covariance());
          }
          trackStateProxy.jacobian() = std::move(jacobian);
          trackStateProxy.pathLength() = std::move(pathLength);

          // Set the filtered parameter index to be the same with predicted
          // parameter
          trackStateProxy.data().ifiltered = trackStateProxy.data().ipredicted;
          // We count the processed state
          ++result.processedStates;
        }
        if (surface->surfaceMaterial() != nullptr) {
          // Update state and stepper with material effects
          materialInteractor(surface, state, stepper, fullUpdate);
        }
      }
      return Result<void>::success();
    }

    /// @brief Kalman actor operation : update in reversed direction
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param surface The surface where the update happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result The mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    Result<void> reversedFilter(const Surface* surface,
                                propagator_state_t& state,
                                const stepper_t& stepper,
                                result_type& result) const {
      const auto& logger = state.options.logger;
      // Try to find the surface in the measurement surfaces
      auto sourcelink_it = inputMeasurements->find(surface->geometryId());
      if (sourcelink_it != inputMeasurements->end()) {
        // Screen output message
        ACTS_VERBOSE("Measurement surface "
                     << surface->geometryId()
                     << " detected in reversed propagation.");

        // No reversed filtering for last measurement state, but still update
        // with material effects
        if (result.reversed and surface == state.navigation.startSurface) {
          materialInteractor(surface, state, stepper);
          return Result<void>::success();
        }

        // Transport the covariance to the surface
        stepper.transportCovarianceToBound(state.stepping, *surface);

        // Update state and stepper with pre material effects
        materialInteractor(surface, state, stepper, preUpdate);

        // Bind the transported state to the current surface
        auto res = stepper.boundState(state.stepping, *surface, false);
        if (!res.ok()) {
          return res.error();
        }

        auto& [boundParams, jacobian, pathLength] = *res;

        // Create a detached track state proxy
        auto tempTrackTip =
            result.fittedStates.addTrackState(TrackStatePropMask::All);

        // Get the detached track state proxy back
        auto trackStateProxy = result.fittedStates.getTrackState(tempTrackTip);

        trackStateProxy.setReferenceSurface(surface->getSharedPtr());

        // Assign the source link to the detached track state
        trackStateProxy.uncalibrated() = sourcelink_it->second;

        // Fill the track state
        trackStateProxy.predicted() = std::move(boundParams.parameters());
        if (boundParams.covariance().has_value()) {
          trackStateProxy.predictedCovariance() =
              std::move(*boundParams.covariance());
        }
        trackStateProxy.jacobian() = std::move(jacobian);
        trackStateProxy.pathLength() = std::move(pathLength);

        // We have predicted parameters, so calibrate the uncalibrated input
        // measuerement
        std::visit(
            [&](const auto& calibrated) {
              trackStateProxy.setCalibrated(calibrated);
            },
            m_calibrator(trackStateProxy.uncalibrated(),
                         trackStateProxy.predicted()));

        // If the update is successful, set covariance and
        auto updateRes = m_updater(state.geoContext, trackStateProxy,
                                   state.stepping.navDir, logger);
        if (!updateRes.ok()) {
          ACTS_ERROR("Backward update step failed: " << updateRes.error());
          return updateRes.error();
        } else {
          // Update the stepping state with filtered parameters
          ACTS_VERBOSE(
              "Backward Filtering step successful, updated parameters are : "
              "\n"
              << trackStateProxy.filtered().transpose());

          // Fill the smoothed parameter for the existing track state
          result.fittedStates.applyBackwards(
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
          materialInteractor(surface, state, stepper, postUpdate);
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
        state.stepping.jacobian = BoundMatrix::Identity();
        if (surface->surfaceMaterial() != nullptr) {
          // Update state and stepper with material effects
          materialInteractor(surface, state, stepper);
        }
      }

      return Result<void>::success();
    }

    /// @brief Kalman actor operation : material interaction
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param surface The surface where the material interaction happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param updateStage The materal update stage
    ///
    template <typename propagator_state_t, typename stepper_t>
    void materialInteractor(
        const Surface* surface, propagator_state_t& state, stepper_t& stepper,
        const MaterialUpdateStage& updateStage = fullUpdate) const {
      const auto& logger = state.options.logger;
      // Indicator if having material
      bool hasMaterial = false;

      if (surface and surface->surfaceMaterial()) {
        // Prepare relevant input particle properties
        detail::PointwiseMaterialInteraction interaction(surface, state,
                                                         stepper);
        // Evaluate the material properties
        if (interaction.evaluateMaterialSlab(state, updateStage)) {
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
          interaction.updateState(state, stepper);
        }
      }

      if (not hasMaterial) {
        // Screen out message
        ACTS_VERBOSE("No material effects on surface: " << surface->geometryId()
                                                        << " at update stage: "
                                                        << updateStage);
      }
    }

    /// @brief Kalman actor operation : finalize
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    Result<void> finalize(propagator_state_t& state, const stepper_t& stepper,
                          result_type& result) const {
      const auto& logger = state.options.logger;
      // Remember you smoothed the track states
      result.smoothed = true;

      // Get the indices of the first states (can be either a measurement or
      // material);
      size_t firstStateIndex = result.lastMeasurementIndex;
      // Count track states to be smoothed
      size_t nStates = 0;
      result.fittedStates.applyBackwards(
          result.lastMeasurementIndex, [&](auto st) {
            bool isMeasurement =
                st.typeFlags().test(TrackStateFlag::MeasurementFlag);
            bool isMaterial = st.typeFlags().test(TrackStateFlag::MaterialFlag);
            if (isMeasurement || isMaterial) {
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
      if (logger().doPrint(Logging::VERBOSE)) {
        ACTS_VERBOSE("Apply smoothing on " << nStates
                                           << " filtered track states.");
      }

      // Smooth the track states
      auto smoothRes = m_smoother(state.geoContext, result.fittedStates,
                                  result.lastMeasurementIndex, logger);
      if (!smoothRes.ok()) {
        ACTS_ERROR("Smoothing step failed: " << smoothRes.error());
        return smoothRes.error();
      }

      // Return in case no target surface
      if (targetSurface == nullptr) {
        return Result<void>::success();
      }

      // Obtain the smoothed parameters at first/last measurement state
      auto firstCreatedState =
          result.fittedStates.getTrackState(firstStateIndex);
      auto lastCreatedMeasurement =
          result.fittedStates.getTrackState(result.lastMeasurementIndex);

      // Lambda to get the intersection of the free params on the target surface
      auto target = [&](const FreeVector& freeVector) -> SurfaceIntersection {
        return targetSurface->intersect(
            state.geoContext, freeVector.segment<3>(eFreePos0),
            state.stepping.navDir * freeVector.segment<3>(eFreeDir0), true);
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
      // At the same time, reverse navigation direction for further
      // stepping if necessary.
      // @note The stepping parameters is updated to the smoothed parameters at
      // either the first measurement state or the last measurement state. It
      // assumes the target surface is not within the first and the last
      // smoothed measurement state. Also, whether the intersection is on
      // surface is not checked here.
      bool reverseDirection = false;
      bool closerTofirstCreatedState =
          (std::abs(firstIntersection.intersection.pathLength) <=
           std::abs(lastIntersection.intersection.pathLength));
      if (closerTofirstCreatedState) {
        stepper.update(state.stepping, firstParams,
                       firstCreatedState.smoothed(),
                       firstCreatedState.smoothedCovariance(),
                       firstCreatedState.referenceSurface());
        reverseDirection = (firstIntersection.intersection.pathLength < 0);
      } else {
        stepper.update(state.stepping, lastParams,
                       lastCreatedMeasurement.smoothed(),
                       lastCreatedMeasurement.smoothedCovariance(),
                       lastCreatedMeasurement.referenceSurface());
        reverseDirection = (lastIntersection.intersection.pathLength < 0);
      }
      const auto& surface = closerTofirstCreatedState
                                ? firstCreatedState.referenceSurface()
                                : lastCreatedMeasurement.referenceSurface();
      ACTS_VERBOSE(
          "Smoothing successful, updating stepping state to smoothed "
          "parameters at surface "
          << surface.geometryId() << ". Prepared to reach the target surface.");

      // Reverse the navigation direction if necessary
      if (reverseDirection) {
        ACTS_VERBOSE(
            "Reverse navigation direction after smoothing for reaching the "
            "target surface");
        state.stepping.navDir =
            (state.stepping.navDir == forward) ? backward : forward;
      }
      // Reset the step size
      state.stepping.stepSize = ConstrainedStep(
          state.stepping.navDir * std::abs(state.options.maxStepSize));
      // Set accumulatd path to zero before targeting surface
      state.stepping.pathAccumulated = 0.;

      return Result<void>::success();
    }

    /// The Kalman updater
    updater_t m_updater;

    /// The Kalman smoother
    smoother_t m_smoother;

    /// The measurement calibrator
    calibrator_t m_calibrator;

    /// The outlier finder
    outlier_finder_t m_outlierFinder;

    /// The Surface beeing
    SurfaceReached targetReached;
  };

  template <typename source_link_t, typename parameters_t,
            typename calibrator_t, typename outlier_finder_t>
  class Aborter {
   public:
    /// Broadcast the result_type
    using action_type =
        Actor<source_link_t, parameters_t, calibrator_t, outlier_finder_t>;

    template <typename propagator_state_t, typename stepper_t,
              typename result_t>
    bool operator()(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const result_t& result) const {
      if (!result.result.ok() or result.finished) {
        return true;
      }
      return false;
    }
  };

 public:
  /// Fit implementation of the foward filter, calls the
  /// the filter and smoother/reversed filter
  ///
  /// @tparam source_link_t Type of the source link
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam calibrator_t Type of the source link calibrator
  /// @tparam outlier_finder_t Type of the outlier finder
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param kfOptions KalmanOptions steering the fit
  /// @note The input measurements are given in the form of @c SourceLinks.
  /// It's
  /// @c calibrator_t's job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename source_link_t, typename start_parameters_t,
            typename calibrator_t, typename outlier_finder_t,
            typename parameters_t = BoundTrackParameters>
  auto fit(const std::vector<source_link_t>& sourcelinks,
           const start_parameters_t& sParameters,
           const KalmanFitterOptions<calibrator_t, outlier_finder_t>& kfOptions)
      const -> std::enable_if_t<!isDirectNavigator,
                                Result<KalmanFitterResult<source_link_t>>> {
    const auto& logger = kfOptions.logger;

    static_assert(SourceLinkConcept<source_link_t>,
                  "Source link does not fulfill SourceLinkConcept");

    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyways, so the map can own them.
    ACTS_VERBOSE("Preparing " << sourcelinks.size() << " input measurements");
    std::map<GeometryIdentifier, source_link_t> inputMeasurements;
    for (const auto& sl : sourcelinks) {
      inputMeasurements.emplace(sl.geometryId(), sl);
    }

    // Create the ActionList and AbortList
    using KalmanAborter =
        Aborter<source_link_t, parameters_t, calibrator_t, outlier_finder_t>;
    using KalmanActor =
        Actor<source_link_t, parameters_t, calibrator_t, outlier_finder_t>;
    using KalmanResult = typename KalmanActor::result_type;
    using Actors = ActionList<KalmanActor>;
    using Aborters = AbortList<KalmanAborter>;

    // Create relevant options for the propagation options
    PropagatorOptions<Actors, Aborters> kalmanOptions(
        kfOptions.geoContext, kfOptions.magFieldContext, logger);

    // Set the trivial propagator options
    kalmanOptions.setPlainOptions(kfOptions.propagatorPlainOptions);

    // Catch the actor and set the measurements
    auto& kalmanActor = kalmanOptions.actionList.template get<KalmanActor>();
    kalmanActor.inputMeasurements = &inputMeasurements;
    kalmanActor.targetSurface = kfOptions.referenceSurface;
    kalmanActor.multipleScattering = kfOptions.multipleScattering;
    kalmanActor.energyLoss = kfOptions.energyLoss;
    kalmanActor.reversedFiltering = kfOptions.reversedFiltering;
    kalmanActor.m_calibrator = kfOptions.calibrator;
    kalmanActor.m_outlierFinder = kfOptions.outlierFinder;

    // Run the fitter
    auto result = m_propagator.template propagate(sParameters, kalmanOptions);

    if (!result.ok()) {
      ACTS_ERROR("Propapation failed: " << result.error());
      return result.error();
    }

    const auto& propRes = *result;

    /// Get the result of the fit
    auto kalmanResult = propRes.template get<KalmanResult>();

    /// It could happen that the fit ends in zero processed states.
    /// The result gets meaningless so such case is regarded as fit failure.
    if (kalmanResult.result.ok() and not kalmanResult.processedStates) {
      kalmanResult.result = Result<void>(KalmanFitterError::NoMeasurementFound);
    }

    if (!kalmanResult.result.ok()) {
      ACTS_ERROR("KalmanFilter failed: "
                 << kalmanResult.result.error() << ", "
                 << kalmanResult.result.error().message());
      return kalmanResult.result.error();
    }

    // Return the converted Track
    return kalmanResult;
  }

  /// Fit implementation of the foward filter, calls the
  /// the filter and smoother/reversed filter
  ///
  /// @tparam source_link_t Type of the source link
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam calibrator_t Type of the source link calibrator
  /// @tparam outlier_finder_t Type of the outlier finder
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param kfOptions KalmanOptions steering the fit
  /// @param sSequence surface sequence used to initialize a DirectNavigator
  /// @note The input measurements are given in the form of @c SourceLinks.
  /// It's
  /// @c calibrator_t's job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename source_link_t, typename start_parameters_t,
            typename calibrator_t, typename outlier_finder_t,
            typename parameters_t = BoundTrackParameters>
  auto fit(const std::vector<source_link_t>& sourcelinks,
           const start_parameters_t& sParameters,
           const KalmanFitterOptions<calibrator_t, outlier_finder_t>& kfOptions,
           const std::vector<const Surface*>& sSequence) const
      -> std::enable_if_t<isDirectNavigator,
                          Result<KalmanFitterResult<source_link_t>>> {
    const auto& logger = kfOptions.logger;
    static_assert(SourceLinkConcept<source_link_t>,
                  "Source link does not fulfill SourceLinkConcept");

    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyways, so the map can own them.
    ACTS_VERBOSE("Preparing " << sourcelinks.size() << " input measurements");
    std::map<GeometryIdentifier, source_link_t> inputMeasurements;
    for (const auto& sl : sourcelinks) {
      inputMeasurements.emplace(sl.geometryId(), sl);
    }

    // Create the ActionList and AbortList
    using KalmanAborter =
        Aborter<source_link_t, parameters_t, calibrator_t, outlier_finder_t>;
    using KalmanActor =
        Actor<source_link_t, parameters_t, calibrator_t, outlier_finder_t>;
    using KalmanResult = typename KalmanActor::result_type;
    using Actors = ActionList<DirectNavigatorInitializer, KalmanActor>;
    using Aborters = AbortList<KalmanAborter>;

    // Create relevant options for the propagation options
    PropagatorOptions<Actors, Aborters> kalmanOptions(
        kfOptions.geoContext, kfOptions.magFieldContext, logger);

    // Set the trivial propagator options
    kalmanOptions.setPlainOptions(kfOptions.propagatorPlainOptions);

    // Catch the actor and set the measurements
    auto& kalmanActor = kalmanOptions.actionList.template get<KalmanActor>();
    kalmanActor.inputMeasurements = &inputMeasurements;
    kalmanActor.targetSurface = kfOptions.referenceSurface;
    kalmanActor.multipleScattering = kfOptions.multipleScattering;
    kalmanActor.energyLoss = kfOptions.energyLoss;
    kalmanActor.reversedFiltering = kfOptions.reversedFiltering;
    kalmanActor.m_calibrator = kfOptions.calibrator;
    // Set config for outlier finder
    kalmanActor.m_outlierFinder = kfOptions.outlierFinder;

    // Set the surface sequence
    auto& dInitializer =
        kalmanOptions.actionList.template get<DirectNavigatorInitializer>();
    dInitializer.navSurfaces = sSequence;

    // Run the fitter
    auto result = m_propagator.template propagate(sParameters, kalmanOptions);

    if (!result.ok()) {
      ACTS_ERROR("Propapation failed: " << result.error());
      return result.error();
    }

    const auto& propRes = *result;

    /// Get the result of the fit
    auto kalmanResult = propRes.template get<KalmanResult>();

    /// It could happen that the fit ends in zero processed states.
    /// The result gets meaningless so such case is regarded as fit failure.
    if (kalmanResult.result.ok() and not kalmanResult.processedStates) {
      kalmanResult.result = Result<void>(KalmanFitterError::NoMeasurementFound);
    }

    if (!kalmanResult.result.ok()) {
      ACTS_ERROR("KalmanFilter failed: " << kalmanResult.result.error());
      return kalmanResult.result.error();
    }

    // Return the converted Track
    return kalmanResult;
  }
};

}  // namespace Acts
