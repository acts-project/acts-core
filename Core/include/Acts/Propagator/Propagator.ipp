// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackParametersConcept.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/PropagatorError.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <concepts>

template <typename S, typename N>
template <typename propagator_state_t>
auto Acts::Propagator<S, N>::propagate(propagator_state_t& state) const
    -> Result<void> {
  ACTS_VERBOSE("Entering propagation.");

  state.stage = PropagatorStage::prePropagation;

  // Pre-Propagation: call to the actor list, abort condition check
  state.options.actorList.act(state, m_stepper, m_navigator, logger());

  if (state.options.actorList.checkAbort(state, m_stepper, m_navigator,
                                         logger())) {
    ACTS_VERBOSE("Propagation terminated without going into stepping loop.");

    state.stage = PropagatorStage::postPropagation;

    state.options.actorList.act(state, m_stepper, m_navigator, logger());

    return Result<void>::success();
  }

  auto getNextTarget = [&]() -> Result<NavigationTarget> {
    for (unsigned int i = 0; i < state.options.maxTargetSkipping; ++i) {
      NavigationTarget nextTarget = m_navigator.nextTarget(
          state.navigation, state.position, state.direction);
      if (nextTarget.isNone()) {
        return NavigationTarget::None();
      }
      IntersectionStatus preStepSurfaceStatus = m_stepper.updateSurfaceStatus(
          state.stepping, *nextTarget.surface,
          nextTarget.surfaceIntersectionIndex, state.options.direction,
          nextTarget.boundaryTolerance, state.options.surfaceTolerance,
          ConstrainedStep::Type::Navigator, logger());
      if (preStepSurfaceStatus == IntersectionStatus::reachable ||
          preStepSurfaceStatus == IntersectionStatus::onSurface) {
        return nextTarget;
      }
    }

    ACTS_ERROR("getNextTarget failed to find a valid target surface after "
               << state.options.maxTargetSkipping << " attempts.");
    return Result<NavigationTarget>::failure(
        PropagatorError::NextTargetLimitReached);
  };

  // priming error condition
  bool terminatedNormally = false;

  // Pre-Stepping: target setting
  state.stage = PropagatorStage::preStep;

  Result<NavigationTarget> nextTargetResult = getNextTarget();
  if (!nextTargetResult.ok()) {
    return nextTargetResult.error();
  }
  NavigationTarget nextTarget = *nextTargetResult;

  ACTS_VERBOSE("Starting stepping loop.");

  // Stepping loop
  for (; state.steps < state.options.maxSteps; ++state.steps) {
    // Perform a step
    Result<double> res =
        m_stepper.step(state.stepping, state.options.direction,
                       m_navigator.currentVolumeMaterial(state.navigation));
    if (!res.ok()) {
      ACTS_ERROR("Step failed with " << res.error() << ": "
                                     << res.error().message());
      // pass error to caller
      return res.error();
    }
    // Accumulate the path length
    state.pathLength += *res;
    // Update the position and direction
    state.position = m_stepper.position(state.stepping);
    state.direction =
        state.options.direction * m_stepper.direction(state.stepping);

    ACTS_VERBOSE("Step with size " << *res << " performed. We are now at "
                                   << state.position.transpose()
                                   << " with direction "
                                   << state.direction.transpose());

    // release actor and aborter constrains after step was performed
    m_stepper.releaseStepSize(state.stepping, ConstrainedStep::Type::Navigator);
    m_stepper.releaseStepSize(state.stepping, ConstrainedStep::Type::Actor);

    // Post-stepping: check target status, call actors, check abort conditions
    state.stage = PropagatorStage::postStep;

    if (!nextTarget.isNone()) {
      IntersectionStatus postStepSurfaceStatus = m_stepper.updateSurfaceStatus(
          state.stepping, *nextTarget.surface,
          nextTarget.surfaceIntersectionIndex, state.options.direction,
          nextTarget.boundaryTolerance, state.options.surfaceTolerance,
          ConstrainedStep::Type::Navigator, logger());
      if (postStepSurfaceStatus == IntersectionStatus::onSurface) {
        m_navigator.handleSurfaceReached(state.navigation, state.position,
                                         state.direction, *nextTarget.surface);
      }
      if (postStepSurfaceStatus != IntersectionStatus::reachable) {
        nextTarget = NavigationTarget::None();
      }
    }

    state.options.actorList.act(state, m_stepper, m_navigator, logger());

    if (state.options.actorList.checkAbort(state, m_stepper, m_navigator,
                                           logger())) {
      terminatedNormally = true;
      break;
    }

    // Update the position and direction because actors might have changed it
    state.position = m_stepper.position(state.stepping);
    state.direction =
        state.options.direction * m_stepper.direction(state.stepping);

    // Pre-Stepping: target setting
    state.stage = PropagatorStage::preStep;

    if (!nextTarget.isNone() &&
        !m_navigator.checkTargetValid(state.navigation, state.position,
                                      state.direction)) {
      ACTS_VERBOSE("Target is not valid anymore.");
      nextTarget = NavigationTarget::None();
    }

    if (nextTarget.isNone()) {
      // navigator step constraint is not valid anymore
      m_stepper.releaseStepSize(state.stepping,
                                ConstrainedStep::Type::Navigator);

      nextTargetResult = getNextTarget();
      if (!nextTargetResult.ok()) {
        return nextTargetResult.error();
      }
      nextTarget = *nextTargetResult;
    }
  }  // end of stepping loop

  // check if we didn't terminate normally via aborters
  if (!terminatedNormally) {
    ACTS_ERROR("Propagation reached the step count limit of "
               << state.options.maxSteps << " (did " << state.steps
               << " steps)");
    return PropagatorError::StepCountLimitReached;
  }

  ACTS_VERBOSE("Stepping loop done.");

  state.stage = PropagatorStage::postPropagation;

  // Post-stepping call to the actor list
  state.options.actorList.act(state, m_stepper, m_navigator, logger());

  return Result<void>::success();
}

template <typename S, typename N>
template <typename parameters_t, typename propagator_options_t,
          typename path_aborter_t>
auto Acts::Propagator<S, N>::propagate(const parameters_t& start,
                                       const propagator_options_t& options,
                                       bool makeCurvilinear) const
    -> Result<
        actor_list_t_result_t<StepperCurvilinearTrackParameters,
                              typename propagator_options_t::actor_list_type>> {
  static_assert(std::copy_constructible<StepperCurvilinearTrackParameters>,
                "return track parameter type must be copy-constructible");

  auto state = makeState(start, options);

  // Perform the actual propagation
  auto propagationResult = propagate(state);

  return makeResult(std::move(state), propagationResult, options,
                    makeCurvilinear);
}

template <typename S, typename N>
template <typename parameters_t, typename propagator_options_t,
          typename target_aborter_t, typename path_aborter_t>
auto Acts::Propagator<S, N>::propagate(
    const parameters_t& start, const Surface& target,
    const propagator_options_t& options) const
    -> Result<
        actor_list_t_result_t<StepperBoundTrackParameters,
                              typename propagator_options_t::actor_list_type>> {
  static_assert(BoundTrackParametersConcept<parameters_t>,
                "Parameters do not fulfill bound parameters concept.");

  auto state = makeState<parameters_t, propagator_options_t, target_aborter_t,
                         path_aborter_t>(start, target, options);

  // Perform the actual propagation
  auto propagationResult = propagate(state);

  return makeResult(std::move(state), propagationResult, target, options);
}

template <typename S, typename N>
template <typename parameters_t, typename propagator_options_t,
          typename path_aborter_t>
auto Acts::Propagator<S, N>::makeState(
    const parameters_t& start, const propagator_options_t& options) const {
  static_assert(BoundTrackParametersConcept<parameters_t>,
                "Parameters do not fulfill bound parameters concept.");

  // Type of track parameters produced by the propagation
  using ReturnParameterType = StepperCurvilinearTrackParameters;

  static_assert(std::copy_constructible<ReturnParameterType>,
                "return track parameter type must be copy-constructible");

  // Expand the abort list with a path aborter
  path_aborter_t pathAborter;
  pathAborter.internalLimit = options.pathLimit;

  auto actorList = options.actorList.append(pathAborter);

  // The expanded options (including path limit)
  auto eOptions = options.extend(actorList);
  eOptions.navigation.startSurface = &start.referenceSurface();
  eOptions.navigation.targetSurface = nullptr;
  using OptionsType = decltype(eOptions);
  using StateType =
      actor_list_t_state_t<OptionsType,
                           typename propagator_options_t::actor_list_type>;
  // Initialize the internal propagator state
  StateType state{eOptions, m_stepper.makeState(eOptions.stepping, start),
                  m_navigator.makeState(eOptions.navigation)};

  initialize<StateType, path_aborter_t>(state);

  return state;
}

template <typename S, typename N>
template <typename parameters_t, typename propagator_options_t,
          typename target_aborter_t, typename path_aborter_t>
auto Acts::Propagator<S, N>::makeState(
    const parameters_t& start, const Surface& target,
    const propagator_options_t& options) const {
  static_assert(BoundTrackParametersConcept<parameters_t>,
                "Parameters do not fulfill bound parameters concept.");

  // Type of provided options
  target_aborter_t targetAborter;
  targetAborter.surface = &target;
  path_aborter_t pathAborter;
  pathAborter.internalLimit = options.pathLimit;
  auto actorList = options.actorList.append(targetAborter, pathAborter);

  // Create the extended options and declare their type
  auto eOptions = options.extend(actorList);
  eOptions.navigation.startSurface = &start.referenceSurface();
  eOptions.navigation.targetSurface = &target;
  using OptionsType = decltype(eOptions);

  // Initialize the internal propagator state
  using StateType =
      actor_list_t_state_t<OptionsType,
                           typename propagator_options_t::actor_list_type>;
  StateType state{eOptions, m_stepper.makeState(eOptions.stepping, start),
                  m_navigator.makeState(eOptions.navigation)};

  initialize<StateType, path_aborter_t>(state);

  return state;
}

template <typename S, typename N>
template <typename propagator_state_t, typename propagator_options_t>
auto Acts::Propagator<S, N>::makeResult(propagator_state_t state,
                                        Result<void> propagationResult,
                                        const propagator_options_t& /*options*/,
                                        bool makeCurvilinear) const
    -> Result<
        actor_list_t_result_t<StepperCurvilinearTrackParameters,
                              typename propagator_options_t::actor_list_type>> {
  // Type of track parameters produced by the propagation
  using ReturnParameterType = StepperCurvilinearTrackParameters;

  static_assert(std::copy_constructible<ReturnParameterType>,
                "return track parameter type must be copy-constructible");

  // Type of the full propagation result, including output from actors
  using ResultType =
      actor_list_t_result_t<ReturnParameterType,
                            typename propagator_options_t::actor_list_type>;

  if (!propagationResult.ok()) {
    return propagationResult.error();
  }

  ResultType result{};
  moveStateToResult(state, result);

  if (makeCurvilinear) {
    if (!m_stepper.prepareCurvilinearState(state.stepping)) {
      // information to compute curvilinearState is incomplete.
      return propagationResult.error();
    }
    /// Convert into return type and fill the result object
    auto curvState = m_stepper.curvilinearState(state.stepping);
    // Fill the end parameters
    result.endParameters =
        std::get<StepperCurvilinearTrackParameters>(curvState);
    // Only fill the transport jacobian when covariance transport was done
    if (state.stepping.covTransport) {
      result.transportJacobian = std::get<Jacobian>(curvState);
    }
  }

  return Result<ResultType>::success(std::move(result));
}

template <typename S, typename N>
template <typename propagator_state_t, typename propagator_options_t>
auto Acts::Propagator<S, N>::makeResult(
    propagator_state_t state, Result<void> propagationResult,
    const Surface& target, const propagator_options_t& /*options*/) const
    -> Result<
        actor_list_t_result_t<StepperBoundTrackParameters,
                              typename propagator_options_t::actor_list_type>> {
  // Type of track parameters produced at the end of the propagation
  using ReturnParameterType = StepperBoundTrackParameters;

  static_assert(std::copy_constructible<ReturnParameterType>,
                "return track parameter type must be copy-constructible");

  // Type of the full propagation result, including output from actors
  using ResultType =
      actor_list_t_result_t<ReturnParameterType,
                            typename propagator_options_t::actor_list_type>;

  if (!propagationResult.ok()) {
    return propagationResult.error();
  }

  ResultType result{};
  moveStateToResult(state, result);

  // Compute the final results and mark the propagation as successful
  auto bsRes = m_stepper.boundState(state.stepping, target);
  if (!bsRes.ok()) {
    return bsRes.error();
  }
  const auto& bs = *bsRes;

  // Fill the end parameters
  result.endParameters = std::get<StepperBoundTrackParameters>(bs);
  // Only fill the transport jacobian when covariance transport was done
  if (state.stepping.covTransport) {
    result.transportJacobian = std::get<Jacobian>(bs);
  }
  return Result<ResultType>::success(std::move(result));
}

template <typename S, typename N>
template <typename propagator_state_t, typename path_aborter_t>
void Acts::Propagator<S, N>::initialize(propagator_state_t& state) const {
  state.position = m_stepper.position(state.stepping);
  state.direction =
      state.options.direction * m_stepper.direction(state.stepping);

  // Navigator initialize state call
  m_navigator.initialize(state.navigation, state.position, state.direction,
                         state.options.direction);

  // Apply the loop protection - it resets the internal path limit
  detail::setupLoopProtection(
      state, m_stepper, state.options.actorList.template get<path_aborter_t>(),
      false, logger());
}

template <typename S, typename N>
template <typename propagator_state_t, typename result_t>
void Acts::Propagator<S, N>::moveStateToResult(propagator_state_t& state,
                                               result_t& result) const {
  result.tuple() = std::move(state.tuple());

  result.steps = state.steps;
  result.pathLength = state.pathLength;

  result.statistics.stepping = state.stepping.statistics;
  result.statistics.navigation = state.navigation.statistics;
}

template <typename derived_t>
Acts::Result<Acts::BoundTrackParameters>
Acts::detail::BasePropagatorHelper<derived_t>::propagateToSurface(
    const BoundTrackParameters& start, const Surface& target,
    const Options& options) const {
  using ResultType = Result<typename derived_t::template actor_list_t_result_t<
      BoundTrackParameters, ActorList<>>>;
  using DerivedOptions = typename derived_t::template Options<>;

  DerivedOptions derivedOptions(options);

  // dummy initialization
  ResultType res = ResultType::failure(PropagatorError::Failure);

  // Due to the geometry of the perigee surface the overstepping tolerance
  // is sometimes not met.
  if (target.type() == Surface::SurfaceType::Perigee) {
    res = static_cast<const derived_t*>(this)
              ->template propagate<BoundTrackParameters, DerivedOptions,
                                   ForcedSurfaceReached, PathLimitReached>(
                  start, target, derivedOptions);
  } else {
    res = static_cast<const derived_t*>(this)
              ->template propagate<BoundTrackParameters, DerivedOptions,
                                   SurfaceReached, PathLimitReached>(
                  start, target, derivedOptions);
  }

  if (res.ok()) {
    // Without errors we can expect a valid endParameters when propagating to a
    // target surface
    assert((*res).endParameters);
    return std::move((*res).endParameters.value());
  } else {
    return res.error();
  }
}
