// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <vector>

namespace Acts {

class Surface;

namespace detail {

/// @brief The step information for recording
///
/// The surface object could be a temporarily created object
/// and as the Step vector is collected and written out at a
/// later stage, the surface is referenced counted here.
struct Step {
  ConstrainedStep stepSize;
  Direction navDir = Direction::Forward();
  Vector3 position = Vector3(0., 0., 0.);
  Vector3 momentum = Vector3(0., 0., 0.);
  std::shared_ptr<const Surface> surface = nullptr;
  GeometryIdentifier geoID;
  /// Note that this is the total number of trials including the previous steps
  std::size_t nTotalTrials = 0;
};

/// @brief a step length logger for debugging the stepping
///
/// It simply logs the constrained step length per step
struct SteppingLogger {
  /// Simple result struct to be returned
  struct this_result {
    std::vector<Step> steps;
  };

  using result_type = this_result;

  /// Set the Logger to sterile
  bool sterile = false;

  /// SteppingLogger action for the ActionList of the Propagator
  ///
  /// @tparam propagator_state_t is the type of Propagator state
  /// @tparam stepper_t is the type of the Stepper
  /// @tparam navigator_t is the type of the Navigator
  ///
  /// @param [in,out] state is the mutable stepper state object
  /// @param [in] stepper the stepper in use
  /// @param [in] navigator the navigator in use
  /// @param [in,out] result is the mutable result object
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  void act(propagator_state_t& state, const stepper_t& stepper,
           const navigator_t& navigator, result_type& result,
           const Logger& /*logger*/) const {
    // Don't log if you have reached the target or are sterile
    if (sterile || state.stage == PropagatorStage::postPropagation) {
      return;
    }

    // Record the propagation state
    Step step;
    step.stepSize = state.stepping.stepSize;
    step.navDir = state.options.direction;
    step.position = stepper.position(state.stepping);
    step.momentum = stepper.momentum(state.stepping);
    step.nTotalTrials = state.stepping.nStepTrials;

    // Record the information about the surface
    if (navigator.currentSurface(state.navigation) != nullptr) {
      step.surface = navigator.currentSurface(state.navigation)->getSharedPtr();
      step.geoID = step.surface->geometryId();
    } else if (navigator.currentVolume(state.navigation) != nullptr) {
      // If there's no surface but a volume, this sets the geoID
      step.geoID = navigator.currentVolume(state.navigation)->geometryId();
    }
    result.steps.push_back(std::move(step));
  }
};

}  // namespace detail
}  // namespace Acts
