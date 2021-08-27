// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <iomanip>
#include <iterator>
#include <sstream>
#include <string>

#include <boost/algorithm/string.hpp>

namespace Acts {

/// DirectNavigator class
///
/// This is a fully guided navigator that progresses through
/// a pre-given sequence of surfaces.
///
/// This can either be used as a validation tool, for truth
/// tracking, or track refitting
class DirectNavigator {
 public:
  /// The sequentially crossed surfaces
  using SurfaceSequence = std::vector<const Surface*>;
  using SurfaceIter = std::vector<const Surface*>::iterator;

  /// Defaulted Constructed
  DirectNavigator() = default;

  /// The tolerance used to define "surface reached"
  double tolerance = s_onSurfaceTolerance;

  /// Nested Actor struct, called Initializer
  ///
  /// This is needed for the initialization of the
  /// surface sequence

  /// Nested State struct
  ///
  /// It acts as an internal state which is
  /// created for every propagation/extrapolation step
  /// and keep thread-local navigation information
  struct State {
    /// Externally provided surfaces - expected to be ordered
    /// along the path
    SurfaceSequence navSurfaces = {};

    /// Iterator the the next surface
    SurfaceIter navSurfaceIter = navSurfaces.begin();

    /// Navigation state - external interface: the start surface
    const Surface* startSurface = nullptr;
    /// Navigation state - external interface: the current surface
    const Surface* currentSurface = nullptr;
    /// Navigation state - external interface: the target surface
    const Surface* targetSurface = nullptr;
    /// Navigation state - starting layer
    const Layer* startLayer = nullptr;
    /// Navigation state - target layer
    const Layer* targetLayer = nullptr;
    /// Navigation state: the start volume
    const TrackingVolume* startVolume = nullptr;
    /// Navigation state: the current volume
    const TrackingVolume* currentVolume = nullptr;
    /// Navigation state: the target volume
    const TrackingVolume* targetVolume = nullptr;

    /// Navigation state - external interface: target is reached
    bool targetReached = false;
    /// Navigation state - external interface: a break has been detected
    bool navigationBreak = false;

    /// Reset state
    ///
    /// @param ssurface is the new starting surface
    /// @param tsurface is the target surface
    void reset(const GeometryContext& /*geoContext*/, const Vector3& /*pos*/,
               const Vector3& /*dir*/, NavigationDirection /*navDir*/,
               const Surface* ssurface, const Surface* tsurface) {
      // Reset everything except the navSurfaces
      State newState = State();
      newState.navSurfaces = this->navSurfaces;
      *this = newState;

      // Reset others
      navSurfaceIter =
          std::find(navSurfaces.begin(), navSurfaces.end(), ssurface);
      startSurface = ssurface;
      currentSurface = ssurface;
      targetSurface = tsurface;
    }
  };

  /// @brief Navigator status call
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void status(propagator_state_t& state, const stepper_t& stepper) const {
    const auto& logger = state.options.logger;
    // Screen output
    ACTS_VERBOSE("Entering navigator::status.");

    // Navigator status always resets the current surface
    state.navigation.currentSurface = nullptr;
    // Output the position in the sequence
    ACTS_VERBOSE(std::distance(state.navigation.navSurfaceIter,
                               state.navigation.navSurfaces.end())
                 << " out of " << state.navigation.navSurfaces.size()
                 << " surfaces remain to try.");

    // Check if we are on surface
    if (state.navigation.navSurfaceIter != state.navigation.navSurfaces.end()) {
      // Establish the surface status
      auto surfaceStatus = stepper.updateSurfaceStatus(
          state.stepping, **state.navigation.navSurfaceIter, false);
      if (surfaceStatus == Intersection3D::Status::onSurface) {
        // Set the current surface
        state.navigation.currentSurface = *state.navigation.navSurfaceIter;
        ACTS_VERBOSE("Current surface set to  "
                     << state.navigation.currentSurface->geometryId())
        // Move the sequence to the next surface
        ++state.navigation.navSurfaceIter;
        if (state.navigation.navSurfaceIter !=
            state.navigation.navSurfaces.end()) {
          ACTS_VERBOSE("Next surface candidate is  "
                       << (*state.navigation.navSurfaceIter)->geometryId());
          stepper.releaseStepSize(state.stepping);
        }
      } else if (surfaceStatus == Intersection3D::Status::reachable) {
        ACTS_VERBOSE("Next surface reachable at distance  "
                     << stepper.outputStepSize(state.stepping));
      }
    }
  }

  /// @brief Navigator target call
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void target(propagator_state_t& state, const stepper_t& stepper) const {
    const auto& logger = state.options.logger;
    // Screen output
    ACTS_VERBOSE("Entering navigator::target.");

    // Navigator target always resets the current surface
    state.navigation.currentSurface = nullptr;
    // Output the position in the sequence
    ACTS_VERBOSE(std::distance(state.navigation.navSurfaceIter,
                               state.navigation.navSurfaces.end())
                 << " out of " << state.navigation.navSurfaces.size()
                 << " surfaces remain to try.");

    if (state.navigation.navSurfaceIter != state.navigation.navSurfaces.end()) {
      // Establish & update the surface status
      auto surfaceStatus = stepper.updateSurfaceStatus(
          state.stepping, **state.navigation.navSurfaceIter, false);
      if (surfaceStatus == Intersection3D::Status::unreachable) {
        ACTS_VERBOSE(
            "Surface not reachable anymore, switching to next one in "
            "sequence");
        // Move the sequence to the next surface
        ++state.navigation.navSurfaceIter;
      } else {
        ACTS_VERBOSE("Navigation stepSize set to "
                     << stepper.outputStepSize(state.stepping));
      }
    } else {
      // Set the navigation break
      state.navigation.navigationBreak = true;
      // If no externally provided target is given, the target is reached
      if (state.navigation.targetSurface == nullptr) {
        state.navigation.targetReached = true;
        // Announce it then
        ACTS_VERBOSE("No target Surface, job done.");
      }
    }
  }
};

}  // namespace Acts
