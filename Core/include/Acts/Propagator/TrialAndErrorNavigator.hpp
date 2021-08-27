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

template <typename propagator_state_t, typename stepper_t>
void trialAndErrorSurfaceUpdate(
    propagator_state_t &state, const stepper_t &stepper,
    const std::vector<const Surface *> &possibleSurfaces) {
  // Prepare the logger
  const auto &logger = state.options.logger;

  // Important parameters
  const FreeVector &params = state.stepping.pars;
  const GeometryContext &gctx = state.geoContext;
  const double oLimit = stepper.overstepLimit(state.stepping);
  const double pLimit = state.options.pathLimit;

  ACTS_VERBOSE("updateOnSurface at "
               << params.segment<3>(eFreePos0).transpose() << " with dir "
               << params.segment<3>(eFreeDir0).transpose());

  // Make intersection objects
  std::vector<SurfaceIntersection> sfis(possibleSurfaces.size());
  std::transform(possibleSurfaces.begin(), possibleSurfaces.end(), sfis.begin(),
                 [&](auto s) {
                   return s->intersect(
                       gctx, params.segment<3>(eFreePos0),
                       state.stepping.navDir * params.segment<3>(eFreeDir0),
                       true);
                 });

  ACTS_VERBOSE("have " << sfis.size() << " candidate intersections");

  // Filter out intersections which are not valid
  auto isIntersectionValid = [&](const SurfaceIntersection &sfi) {
    const double sfiPath = sfi.intersection.pathLength;
    return sfi.intersection.status == Intersection3D::Status::reachable &&
           sfiPath > oLimit && sfiPath * sfiPath <= pLimit * pLimit;
  };

  sfis.erase(std::remove_if(
                 sfis.begin(), sfis.end(),
                 [&](const auto &sfi) { return not isIntersectionValid(sfi); }),
             sfis.end());

  ACTS_VERBOSE("after remove " << sfis.size() << " candidates remain");

  std::sort(sfis.begin(), sfis.end());
  for (auto &sfi : sfis) {
    sfi.intersection.pathLength *= std::copysign(1., state.stepping.navDir);
  }

  // Transform back to vector of surface pointers
  std::vector<const Surface *> finalSurfaces;
  finalSurfaces.reserve(sfis.size());
  for (const auto &si : sfis) {
    finalSurfaces.push_back(si.representation);
  }

  state.navigation.navSurfaces = std::move(finalSurfaces);
  state.navigation.navSurfaceIter = state.navigation.navSurfaces.begin();
}

struct TrialAndErrorGuider {
  std::vector<const Surface *> possibleSurfaces;

  template <typename propagator_state_t, typename stepper_t>
  void updateOnSurface(propagator_state_t &state,
                       const stepper_t &stepper) const {
    trialAndErrorSurfaceUpdate(state, stepper, possibleSurfaces);
  }
};

}  // namespace detail

using TrialAndErrorNavigator = GuidedNavigator<detail::TrialAndErrorGuider>;

struct TrialAndErrorNavigatorInitializer {
  /// Pointer to the surface list
  std::vector<const Surface *> *surfaces = nullptr;

  /// Actor result / state
  struct this_result {
    bool initialized = false;
  };

  using result_type = this_result;

  /// Defaulting the constructor
  TrialAndErrorNavigatorInitializer() = default;

  /// Actor operator call
  /// @tparam statet Type of the full propagator state
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param state the entire propagator state
  /// @param r the result of this Actor
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t &state, const stepper_t &stepper,
                  result_type &r) const {
    // Only act once
    if (not surfaces) {
      const auto &logger = state.options.logger;
      ACTS_ERROR(
          "No surfaces are provided to the TrialAndErrorNavigatorInitializer");
    }

    if (not r.initialized) {
      trialAndErrorSurfaceUpdate(state, stepper, *surfaces);
      r.initialized = true;
    }
  }

  /// Actor operator call - resultless, unused
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t & /*unused*/,
                  const stepper_t & /*unused*/) const {}
};

}  // namespace Acts
