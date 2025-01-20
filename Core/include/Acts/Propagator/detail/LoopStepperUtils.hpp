// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp>
#include <Acts/Propagator/detail/CovarianceEngine.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Result.hpp>

namespace Acts::detail {

/// A helper type for providinig a propagation state which can be used with
/// functions expecting single-component steppers and states
template <typename stepping_t, typename navigation_t, typename options_t,
          typename geoctx_t>
struct SinglePropState {
  stepping_t& stepping;
  navigation_t& navigation;
  options_t& options;
  geoctx_t& geoContext;

  SinglePropState(stepping_t& s, navigation_t& n, options_t& o, geoctx_t& g)
      : stepping(s), navigation(n), options(o), geoContext(g) {}
};

/// A template class which contains all const member functions, that should be
/// available both in the mutable ComponentProxy and the ConstComponentProxy.
/// @tparam component_t Must be a const or mutable State::Component.
template <typename component_t, typename loop_stepper_t>
struct LoopComponentProxyBase {
  using SingleStepper = typename loop_stepper_t::SingleStepper;
  using SingleState = typename loop_stepper_t::SingleState;

  static_assert(std::is_same_v<std::remove_const_t<component_t>,
                               typename loop_stepper_t::State::Component>);

  component_t& cmp;

  LoopComponentProxyBase(component_t& c) : cmp(c) {}

  // These are the const accessors, which are shared between the mutable
  // ComponentProxy and the ConstComponentProxy
  const auto& state() const { return cmp.state; }
  auto status() const { return cmp.status; }
  auto weight() const { return cmp.weight; }
  auto pathAccumulated() const { return cmp.state.pathAccumulated; }
  const auto& pars() const { return cmp.state.pars; }
  const auto& derivative() const { return cmp.state.derivative; }
  const auto& jacTransport() const { return cmp.state.jacTransport; }
  const auto& cov() const { return cmp.state.cov; }
  const auto& jacobian() const { return cmp.state.jacobian; }
  const auto& jacToGlobal() const { return cmp.state.jacToGlobal; }

  template <typename propagator_state_t>
  auto singleState(const propagator_state_t& state) const {
    using DeducedStepping = decltype(state.stepping.components.front().state);
    static_assert(std::is_same_v<SingleState, DeducedStepping>);

    return SinglePropState<const SingleState, const decltype(state.navigation),
                           const decltype(state.options),
                           const decltype(state.geoContext)>(
        cmp.state, state.navigation, state.options, state.geoContext);
  }

  const auto& singleStepper(const loop_stepper_t& stepper) const {
    return static_cast<const SingleStepper&>(stepper);
  }
};

/// A proxy struct which allows access to a single component of the
/// multi-component state. It has the semantics of a mutable reference, i.e.
/// it requires a mutable reference of the single-component state it
/// represents
template <typename component_t, typename loop_stepper_t>
struct LoopComponentProxy
    : LoopComponentProxyBase<component_t, loop_stepper_t> {
  using State = typename loop_stepper_t::State;
  using Base = LoopComponentProxyBase<component_t, loop_stepper_t>;

  using SingleState = typename loop_stepper_t::SingleState;
  using SingleStepper = typename loop_stepper_t::SingleStepper;
  using Covariance = typename loop_stepper_t::Covariance;

  // Import the const accessors from ComponentProxyBase
  using Base::cmp;
  using Base::cov;
  using Base::derivative;
  using Base::jacobian;
  using Base::jacToGlobal;
  using Base::jacTransport;
  using Base::pars;
  using Base::pathAccumulated;
  using Base::singleState;
  using Base::singleStepper;
  using Base::state;
  using Base::status;
  using Base::weight;

  // The multi-component state of the stepper
  const State& all_state;

  LoopComponentProxy(typename State::Component& c, const State& s)
      : Base(c), all_state(s) {}

  // These are the mutable accessors, the const ones are inherited from the
  // ComponentProxyBase
  auto& state() { return cmp.state; }
  auto& status() { return cmp.status; }
  auto& weight() { return cmp.weight; }
  auto& pathAccumulated() { return cmp.state.pathAccumulated; }
  auto& pars() { return cmp.state.pars; }
  auto& derivative() { return cmp.state.derivative; }
  auto& jacTransport() { return cmp.state.jacTransport; }
  auto& cov() { return cmp.state.cov; }
  auto& jacobian() { return cmp.state.jacobian; }
  auto& jacToGlobal() { return cmp.state.jacToGlobal; }

  template <typename propagator_state_t>
  auto singleState(propagator_state_t& state) {
    using DeducedStepping = decltype(state.stepping.components.front().state);
    static_assert(std::is_same_v<SingleState, DeducedStepping>);

    return SinglePropState<SingleState, decltype(state.navigation),
                           decltype(state.options), decltype(state.geoContext)>(
        cmp.state, state.navigation, state.options, state.geoContext);
  }

  Result<typename SingleStepper::BoundState> boundState(
      const Surface& surface, bool transportCov,
      const FreeToBoundCorrection& freeToBoundCorrection) {
    return detail::boundState(
        all_state.options.geoContext, surface, cov(), jacobian(),
        jacTransport(), derivative(), jacToGlobal(), pars(),
        all_state.particleHypothesis, all_state.covTransport && transportCov,
        cmp.state.pathAccumulated, freeToBoundCorrection);
  }

  void update(const FreeVector& freeParams, const BoundVector& boundParams,
              const Covariance& covariance, const Surface& surface) {
    cmp.state.pars = freeParams;
    cmp.state.cov = covariance;
    cmp.state.jacToGlobal =
        surface.boundToFreeJacobian(all_state.geoContext, boundParams);
  }
};

}  // namespace Acts::detail
