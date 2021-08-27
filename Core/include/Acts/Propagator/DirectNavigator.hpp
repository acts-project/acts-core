// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/GuidedNavigator.hpp"

namespace Acts {
namespace detail {

struct DirectGuider {
  template <typename propagator_state_t, typename stepper_t>
  void updateOnSurface(propagator_state_t& state, const stepper_t&) const {
    ++state.navigation.navSurfaceIter;
  }
};

}  // namespace detail

using DirectNavigator = GuidedNavigator<detail::DirectGuider>;

struct DirectNavigatorInitializer {
  /// The Surface sequence
  DirectNavigator::SurfaceSequence navSurfaces = {};

  /// Actor result / state
  struct this_result {
    bool initialized = false;
  };
  using result_type = this_result;

  /// Defaulting the constructor
  DirectNavigatorInitializer() = default;

  /// Actor operator call
  /// @tparam statet Type of the full propagator state
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param state the entire propagator state
  /// @param r the result of this Actor
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& /*unused*/,
                  result_type& r) const {
    // Only act once
    if (not r.initialized) {
      // Initialize the surface sequence
      state.navigation.navSurfaces = navSurfaces;
      state.navigation.navSurfaceIter = state.navigation.navSurfaces.begin();
      r.initialized = true;
    }
  }

  /// Actor operator call - resultless, unused
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& /*unused*/,
                  const stepper_t& /*unused*/) const {}
};

}  // namespace Acts
