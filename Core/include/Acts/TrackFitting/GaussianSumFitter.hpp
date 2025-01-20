// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/MultiStepperAborters.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/TrackFitting/GsfError.hpp"
#include "Acts/TrackFitting/GsfOptions.hpp"
#include "Acts/TrackFitting/detail/GsfActor.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"

namespace Acts {

namespace detail {

/// Type trait to identify if a type is a MultiComponentBoundTrackParameters and
/// to inspect its charge representation if not TODO this probably gives an ugly
/// error message if detectCharge does not compile
template <typename T>
struct IsMultiComponentBoundParameters : public std::false_type {};

template <>
struct IsMultiComponentBoundParameters<MultiComponentBoundTrackParameters>
    : public std::true_type {};

}  // namespace detail

/// Gaussian Sum Fitter implementation.
/// @tparam propagator_t The propagator type on which the algorithm is built on
/// @tparam bethe_heitler_approx_t The type of the Bethe-Heitler-Approximation
/// @tparam traj_t The MultiTrajectory type (backend)
///
/// @note This GSF implementation tries to be as compatible to the KalmanFitter
/// as possible. However, strict compatibility is not garantueed.
/// @note Currently there is no possibility to export the states of the
/// individual components from the GSF, the only information returned in the
/// MultiTrajectory are the means of the states. Therefore, also NO dedicated
/// component smoothing is performed as described e.g. by R. Fruewirth.
template <typename propagator_t, typename bethe_heitler_approx_t,
          typename traj_t>
struct GaussianSumFitter {
  GaussianSumFitter(propagator_t&& propagator, bethe_heitler_approx_t&& bha,
                    std::unique_ptr<const Logger> _logger =
                        getDefaultLogger("GSF", Logging::INFO))
      : m_propagator(std::move(propagator)),
        m_betheHeitlerApproximation(std::move(bha)),
        m_logger{std::move(_logger)},
        m_actorLogger(m_logger->cloneWithSuffix("Actor")) {}

  /// The propagator instance used by the fit function
  propagator_t m_propagator;

  /// The fitter holds the instance of the bethe heitler approx
  bethe_heitler_approx_t m_betheHeitlerApproximation;

  /// The logger
  std::unique_ptr<const Logger> m_logger;
  std::unique_ptr<const Logger> m_actorLogger;

  const Logger& logger() const { return *m_logger; }

  /// The navigator type
  using GsfNavigator = typename propagator_t::Navigator;

  /// The actor type
  using GsfActor = detail::GsfActor<bethe_heitler_approx_t, traj_t>;

  /// @brief The fit function for the Direct navigator
  template <typename source_link_it_t, typename start_parameters_t,
            TrackContainerFrontend track_container_t>
  auto fit(source_link_it_t begin, source_link_it_t end,
           const start_parameters_t& sParameters,
           const GsfOptions<traj_t>& options,
           const std::vector<const Surface*>& sSequence,
           track_container_t& trackContainer) const {
    // Check if we have the correct navigator
    static_assert(
        std::is_same_v<DirectNavigator, typename propagator_t::Navigator>);

    // Initialize the forward propagation with the DirectNavigator
    auto fwdPropInitializer = [&sSequence, this](const auto& opts) {
      using Actors = ActorList<GsfActor>;
      using PropagatorOptions = typename propagator_t::template Options<Actors>;

      PropagatorOptions propOptions(opts.geoContext, opts.magFieldContext);

      propOptions.setPlainOptions(opts.propagatorPlainOptions);

      propOptions.navigation.surfaces = sSequence;
      propOptions.actorList.template get<GsfActor>()
          .m_cfg.bethe_heitler_approx = &m_betheHeitlerApproximation;

      return propOptions;
    };

    // Initialize the backward propagation with the DirectNavigator
    auto bwdPropInitializer = [&sSequence, this](const auto& opts) {
      using Actors = ActorList<GsfActor>;
      using PropagatorOptions = typename propagator_t::template Options<Actors>;

      PropagatorOptions propOptions(opts.geoContext, opts.magFieldContext);

      propOptions.setPlainOptions(opts.propagatorPlainOptions);

      propOptions.navigation.surfaces = sSequence;
      propOptions.actorList.template get<GsfActor>()
          .m_cfg.bethe_heitler_approx = &m_betheHeitlerApproximation;

      return propOptions;
    };

    return fit_impl(begin, end, sParameters, options, fwdPropInitializer,
                    bwdPropInitializer, trackContainer);
  }

  /// @brief The fit function for the standard navigator
  template <typename source_link_it_t, typename start_parameters_t,
            TrackContainerFrontend track_container_t>
  auto fit(source_link_it_t begin, source_link_it_t end,
           const start_parameters_t& sParameters,
           const GsfOptions<traj_t>& options,
           track_container_t& trackContainer) const {
    // Check if we have the correct navigator
    static_assert(std::is_same_v<Navigator, typename propagator_t::Navigator>);

    // Initialize the forward propagation with the DirectNavigator
    auto fwdPropInitializer = [this](const auto& opts) {
      using Actors = ActorList<GsfActor, EndOfWorldReached>;
      using PropagatorOptions = typename propagator_t::template Options<Actors>;

      PropagatorOptions propOptions(opts.geoContext, opts.magFieldContext);

      propOptions.setPlainOptions(opts.propagatorPlainOptions);

      propOptions.actorList.template get<GsfActor>()
          .m_cfg.bethe_heitler_approx = &m_betheHeitlerApproximation;

      return propOptions;
    };

    // Initialize the backward propagation with the DirectNavigator
    auto bwdPropInitializer = [this](const auto& opts) {
      using Actors = ActorList<GsfActor, EndOfWorldReached>;
      using PropagatorOptions = typename propagator_t::template Options<Actors>;

      PropagatorOptions propOptions(opts.geoContext, opts.magFieldContext);

      propOptions.setPlainOptions(opts.propagatorPlainOptions);

      propOptions.actorList.template get<GsfActor>()
          .m_cfg.bethe_heitler_approx = &m_betheHeitlerApproximation;

      return propOptions;
    };

    return fit_impl(begin, end, sParameters, options, fwdPropInitializer,
                    bwdPropInitializer, trackContainer);
  }

  /// The generic implementation of the fit function.
  /// TODO check what this function does with the referenceSurface is e.g. the
  /// first measurementSurface
  template <typename source_link_it_t, typename start_parameters_t,
            typename fwd_prop_initializer_t, typename bwd_prop_initializer_t,
            TrackContainerFrontend track_container_t>
  Acts::Result<typename track_container_t::TrackProxy> fit_impl(
      source_link_it_t begin, source_link_it_t end,
      const start_parameters_t& sParameters, const GsfOptions<traj_t>& options,
      const fwd_prop_initializer_t& fwdPropInitializer,
      const bwd_prop_initializer_t& bwdPropInitializer,
      track_container_t& trackContainer) const {
    // return or abort utility
    auto return_error_or_abort = [&](auto error) {
      if (options.abortOnError) {
        std::abort();
      }
      return error;
    };

    // Define directions based on input propagation direction. This way we can
    // refer to 'forward' and 'backward' regardless of the actual direction.
    const auto gsfForward = options.propagatorPlainOptions.direction;
    const auto gsfBackward = gsfForward.invert();

    // Check if the start parameters are on the start surface
    auto intersectionStatusStartSurface =
        sParameters.referenceSurface()
            .intersect(GeometryContext{},
                       sParameters.position(GeometryContext{}),
                       sParameters.direction(), BoundaryTolerance::None())
            .closest()
            .status();

    if (intersectionStatusStartSurface != IntersectionStatus::onSurface) {
      ACTS_DEBUG(
          "Surface intersection of start parameters WITH bound-check failed");
    }

    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyway, so the map can own them.
    ACTS_VERBOSE("Preparing " << std::distance(begin, end)
                              << " input measurements");
    std::map<GeometryIdentifier, SourceLink> inputMeasurements;
    for (auto it = begin; it != end; ++it) {
      SourceLink sl = *it;
      inputMeasurements.emplace(
          options.extensions.surfaceAccessor(sl)->geometryId(), std::move(sl));
    }

    ACTS_VERBOSE(
        "Gsf: Final measurement map size: " << inputMeasurements.size());

    if (sParameters.covariance() == std::nullopt) {
      return GsfError::StartParametersHaveNoCovariance;
    }

    /////////////////
    // Forward pass
    /////////////////
    ACTS_VERBOSE("+-----------------------------+");
    ACTS_VERBOSE("| Gsf: Do forward propagation |");
    ACTS_VERBOSE("+-----------------------------+");

    auto fwdResult = [&]() {
      auto fwdPropOptions = fwdPropInitializer(options);

      // Catch the actor and set the measurements
      auto& actor = fwdPropOptions.actorList.template get<GsfActor>();
      actor.setOptions(options);
      actor.m_cfg.inputMeasurements = &inputMeasurements;
      actor.m_cfg.numberMeasurements = inputMeasurements.size();
      actor.m_cfg.inReversePass = false;
      actor.m_cfg.logger = m_actorLogger.get();

      fwdPropOptions.direction = gsfForward;

      // If necessary convert to MultiComponentBoundTrackParameters
      using IsMultiParameters =
          detail::IsMultiComponentBoundParameters<start_parameters_t>;

      // dirty optional because parameters are not default constructible
      std::optional<MultiComponentBoundTrackParameters> params;

      // This allows the initialization with single- and multicomponent start
      // parameters
      if constexpr (!IsMultiParameters::value) {
        params = MultiComponentBoundTrackParameters(
            sParameters.referenceSurface().getSharedPtr(),
            sParameters.parameters(), *sParameters.covariance(),
            sParameters.particleHypothesis());
      } else {
        params = sParameters;
      }

      auto state = m_propagator.makeState(*params, fwdPropOptions);

      auto& r = state.template get<typename GsfActor::result_type>();
      r.fittedStates = &trackContainer.trackStateContainer();

      auto propagationResult = m_propagator.propagate(state);

      return m_propagator.makeResult(std::move(state), propagationResult,
                                     fwdPropOptions, false);
    }();

    if (!fwdResult.ok()) {
      return return_error_or_abort(fwdResult.error());
    }

    const auto& fwdGsfResult =
        fwdResult->template get<typename GsfActor::result_type>();

    if (!fwdGsfResult.result.ok()) {
      return return_error_or_abort(fwdGsfResult.result.error());
    }

    if (fwdGsfResult.measurementStates == 0) {
      return return_error_or_abort(GsfError::NoMeasurementStatesCreatedForward);
    }

    ACTS_VERBOSE("Finished forward propagation");
    ACTS_VERBOSE("- visited surfaces: " << fwdGsfResult.visitedSurfaces.size());
    ACTS_VERBOSE("- processed states: " << fwdGsfResult.processedStates);
    ACTS_VERBOSE("- measurement states: " << fwdGsfResult.measurementStates);

    std::size_t nInvalidBetheHeitler = fwdGsfResult.nInvalidBetheHeitler.val();
    double maxPathXOverX0 = fwdGsfResult.maxPathXOverX0.val();

    //////////////////
    // Backward pass
    //////////////////
    ACTS_VERBOSE("+------------------------------+");
    ACTS_VERBOSE("| Gsf: Do backward propagation |");
    ACTS_VERBOSE("+------------------------------+");

    auto bwdResult = [&]() {
      auto bwdPropOptions = bwdPropInitializer(options);

      auto& actor = bwdPropOptions.actorList.template get<GsfActor>();
      actor.setOptions(options);
      actor.m_cfg.inputMeasurements = &inputMeasurements;
      actor.m_cfg.inReversePass = true;
      actor.m_cfg.logger = m_actorLogger.get();

      bwdPropOptions.direction = gsfBackward;

      const Surface& target = options.referenceSurface
                                  ? *options.referenceSurface
                                  : sParameters.referenceSurface();

      const auto& params = *fwdGsfResult.lastMeasurementState;
      std::vector<
          std::tuple<double, BoundVector, std::optional<BoundSquareMatrix>>>
          inflatedParamVector;
      for (auto& [w, p, cov] : params.components()) {
        inflatedParamVector.emplace_back(
            w, p, cov.value() * options.reverseFilteringCovarianceScaling);
      }

      MultiComponentBoundTrackParameters inflatedParams(
          params.referenceSurface().getSharedPtr(),
          std::move(inflatedParamVector), params.particleHypothesis());

      auto state =
          m_propagator.template makeState<MultiComponentBoundTrackParameters,
                                          decltype(bwdPropOptions),
                                          MultiStepperSurfaceReached>(
              inflatedParams, target, bwdPropOptions);

      assert(
          (fwdGsfResult.lastMeasurementTip != MultiTrajectoryTraits::kInvalid &&
           "tip is invalid"));

      auto proxy = trackContainer.trackStateContainer().getTrackState(
          fwdGsfResult.lastMeasurementTip);
      proxy.shareFrom(TrackStatePropMask::Filtered,
                      TrackStatePropMask::Smoothed);

      auto& r = state.template get<typename GsfActor::result_type>();
      r.fittedStates = &trackContainer.trackStateContainer();
      r.currentTip = fwdGsfResult.lastMeasurementTip;
      r.visitedSurfaces.push_back(&proxy.referenceSurface());
      r.surfacesVisitedBwdAgain.push_back(&proxy.referenceSurface());
      r.measurementStates++;
      r.processedStates++;

      auto propagationResult = m_propagator.propagate(state);

      return m_propagator.makeResult(std::move(state), propagationResult,
                                     target, bwdPropOptions);
    }();

    if (!bwdResult.ok()) {
      return return_error_or_abort(bwdResult.error());
    }

    auto& bwdGsfResult =
        bwdResult->template get<typename GsfActor::result_type>();

    if (!bwdGsfResult.result.ok()) {
      return return_error_or_abort(bwdGsfResult.result.error());
    }

    if (bwdGsfResult.measurementStates == 0) {
      return return_error_or_abort(
          GsfError::NoMeasurementStatesCreatedBackward);
    }

    // For the backward pass we want the counters at in end (= at the
    // interaction point) and not at the last measurement surface
    bwdGsfResult.nInvalidBetheHeitler.update();
    bwdGsfResult.maxPathXOverX0.update();
    bwdGsfResult.sumPathXOverX0.update();
    nInvalidBetheHeitler += bwdGsfResult.nInvalidBetheHeitler.val();
    maxPathXOverX0 =
        std::max(maxPathXOverX0, bwdGsfResult.maxPathXOverX0.val());

    if (nInvalidBetheHeitler > 0) {
      ACTS_WARNING("Encountered " << nInvalidBetheHeitler
                                  << " cases where x/X0 exceeds the range "
                                     "of the Bethe-Heitler-Approximation. The "
                                     "maximum x/X0 encountered was "
                                  << maxPathXOverX0
                                  << ". Enable DEBUG output "
                                     "for more information.");
    }

    ////////////////////////////////////
    // Create Kalman Result
    ////////////////////////////////////
    ACTS_VERBOSE("Gsf - States summary:");
    ACTS_VERBOSE("- Fwd measurement states: " << fwdGsfResult.measurementStates
                                              << ", holes: "
                                              << fwdGsfResult.measurementHoles);
    ACTS_VERBOSE("- Bwd measurement states: " << bwdGsfResult.measurementStates
                                              << ", holes: "
                                              << bwdGsfResult.measurementHoles);

    // TODO should this be warning level? it happens quite often... Investigate!
    if (bwdGsfResult.measurementStates != fwdGsfResult.measurementStates) {
      ACTS_DEBUG("Fwd and bwd measurement states do not match");
    }

    // Go through the states and assign outliers / unset smoothed if surface not
    // passed in backward pass
    const auto& foundBwd = bwdGsfResult.surfacesVisitedBwdAgain;
    std::size_t measurementStatesFinal = 0;

    for (auto state : fwdGsfResult.fittedStates->reverseTrackStateRange(
             fwdGsfResult.currentTip)) {
      const bool found =
          rangeContainsValue(foundBwd, &state.referenceSurface());
      if (!found && state.typeFlags().test(MeasurementFlag)) {
        state.typeFlags().set(OutlierFlag);
        state.typeFlags().reset(MeasurementFlag);
        state.unset(TrackStatePropMask::Smoothed);
      }

      measurementStatesFinal +=
          static_cast<std::size_t>(state.typeFlags().test(MeasurementFlag));
    }

    if (measurementStatesFinal == 0) {
      return return_error_or_abort(GsfError::NoMeasurementStatesCreatedFinal);
    }

    auto track = trackContainer.makeTrack();
    track.tipIndex() = fwdGsfResult.lastMeasurementTip;

    if (options.referenceSurface) {
      const auto& params = *bwdResult->endParameters;

      const auto [finalPars, finalCov] = detail::mergeGaussianMixture(
          params.components(), params.referenceSurface(),
          options.componentMergeMethod, [](auto& t) {
            return std::tie(std::get<0>(t), std::get<1>(t), *std::get<2>(t));
          });

      track.parameters() = finalPars;
      track.covariance() = finalCov;

      track.setReferenceSurface(params.referenceSurface().getSharedPtr());

      if (trackContainer.hasColumn(
              hashString(GsfConstants::kFinalMultiComponentStateColumn))) {
        ACTS_DEBUG("Add final multi-component state to track");
        track.template component<GsfConstants::FinalMultiComponentState>(
            GsfConstants::kFinalMultiComponentStateColumn) = std::move(params);
      }
    }

    if (trackContainer.hasColumn(
            hashString(GsfConstants::kFwdMaxMaterialXOverX0))) {
      track.template component<double>(GsfConstants::kFwdMaxMaterialXOverX0) =
          fwdGsfResult.maxPathXOverX0.val();
    }
    if (trackContainer.hasColumn(
            hashString(GsfConstants::kFwdSumMaterialXOverX0))) {
      track.template component<double>(GsfConstants::kFwdSumMaterialXOverX0) =
          fwdGsfResult.sumPathXOverX0.val();
    }

    calculateTrackQuantities(track);

    return track;
  }
};

}  // namespace Acts
