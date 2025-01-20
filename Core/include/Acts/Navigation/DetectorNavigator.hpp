// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Propagator/NavigatorStatistics.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <sstream>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/container/small_vector.hpp>

namespace Acts::Experimental {

class DetectorNavigator {
 public:
  struct Config {
    /// Detector for this Navigation
    const Detector* detector = nullptr;

    /// Configuration for this Navigator
    /// stop at every sensitive surface (whether it has material or not)
    bool resolveSensitive = true;
    /// stop at every material surface (whether it is passive or not)
    bool resolveMaterial = true;
    /// stop at every surface regardless what it is
    bool resolvePassive = false;
  };

  struct Options : public NavigatorPlainOptions {
    explicit Options(const GeometryContext& gctx)
        : NavigatorPlainOptions(gctx) {}

    void setPlainOptions(const NavigatorPlainOptions& options) {
      static_cast<NavigatorPlainOptions&>(*this) = options;
    }
  };

  /// Nested State struct
  ///
  /// It acts as an internal state which is
  /// created for every propagation/extrapolation step
  /// and keep thread-local navigation information
  struct State : public NavigationState {
    explicit State(const Options& options_) : options(options_) {}

    Options options;

    /// Navigation state - external state: the current surface
    const Surface* currentSurface = nullptr;
    /// Navigation state : a break has been detected
    bool navigationBreak = false;

    /// Navigation statistics
    NavigatorStatistics statistics;
  };

  /// Constructor with configuration object
  ///
  /// @param cfg The navigator configuration
  /// @param _logger a logger instance
  explicit DetectorNavigator(Config cfg,
                             std::shared_ptr<const Logger> _logger =
                                 getDefaultLogger("DetectorNavigator",
                                                  Logging::Level::INFO))
      : m_cfg{cfg}, m_logger{std::move(_logger)} {}

  State makeState(const Options& options) const {
    State state(options);
    return state;
  }

  const Surface* currentSurface(const State& state) const {
    return state.currentSurface;
  }

  const DetectorVolume* currentVolume(const State& state) const {
    return state.currentVolume;
  }

  const IVolumeMaterial* currentVolumeMaterial(const State& state) const {
    return state.currentVolume->volumeMaterial();
  }

  const Surface* startSurface(const State& state) const {
    return state.options.startSurface;
  }

  const Surface* targetSurface(const State& state) const {
    return state.options.targetSurface;
  }

  bool endOfWorldReached(State& state) const {
    return state.currentVolume == nullptr;
  }

  bool navigationBreak(const State& state) const {
    return state.navigationBreak;
  }

  Result<void> initialize(State& state, const Vector3& position,
                          const Vector3& direction,
                          Direction propagationDirection) const {
    (void)propagationDirection;

    ACTS_VERBOSE(volInfo(state) << posInfo(state, position) << "initialize");

    if (state.currentDetector == nullptr) {
      ACTS_VERBOSE("Assigning detector from the config.");
      state.currentDetector = m_cfg.detector;
    }
    if (state.currentDetector == nullptr) {
      throw std::invalid_argument("DetectorNavigator: no detector assigned");
    }

    fillNavigationState(position, direction, state);
    if (state.currentVolume == nullptr) {
      state.currentVolume = state.currentDetector->findDetectorVolume(
          state.options.geoContext, state.position);
    }
    if (state.currentVolume == nullptr) {
      throw std::invalid_argument("DetectorNavigator: no current volume found");
    }
    updateCandidateSurfaces(state, position);

    return Result<void>::success();
  }

  NavigationTarget nextTarget(State& state, const Vector3& position,
                              const Vector3& direction) const {
    ACTS_VERBOSE(volInfo(state)
                 << posInfo(state, position) << "Entering navigator::preStep.");

    if (inactive()) {
      ACTS_VERBOSE(volInfo(state)
                   << posInfo(state, position) << "navigator inactive");
      return NavigationTarget::None();
    }

    fillNavigationState(position, direction, state);

    if (state.currentSurface != nullptr) {
      ACTS_VERBOSE(volInfo(state)
                   << posInfo(state, position) << "stepping through surface");
    }

    if (state.surfaceCandidateIndex == state.surfaceCandidates.size()) {
      ACTS_VERBOSE(volInfo(state)
                   << posInfo(state, position) << "no surface candidates");
      return NavigationTarget::None();
    }

    // Screen output how much is left to try
    ACTS_VERBOSE(volInfo(state) << posInfo(state, position)
                                << (state.surfaceCandidates.size() -
                                    state.surfaceCandidateIndex)
                                << " out of " << state.surfaceCandidates.size()
                                << " surfaces remain to try.");
    // Take the surface
    const auto& candidate = state.surfaceCandidate();
    const auto& surface = (candidate.surface != nullptr)
                              ? (*candidate.surface)
                              : (candidate.portal->surface());
    // Screen output which surface you are on
    ACTS_VERBOSE(volInfo(state)
                 << posInfo(state, position)
                 << "next surface candidate will be " << surface.geometryId()
                 << " (" << surface.center(state.options.geoContext).transpose()
                 << ")");

    state.currentSurface = nullptr;
    state.currentPortal = nullptr;

    return NavigationTarget(surface, candidate.objectIntersection.index(),
                            candidate.boundaryTolerance);
  }

  bool checkTargetValid(const State& state, const Vector3& position,
                        const Vector3& direction) const {
    (void)state;
    (void)position;
    (void)direction;

    return true;
  }

  void handleSurfaceReached(State& state, const Vector3& position,
                            const Vector3& direction,
                            const Surface& surface) const {
    (void)surface;

    ACTS_VERBOSE(volInfo(state) << posInfo(state, position)
                                << "Entering navigator::handleSurfaceReached.");

    fillNavigationState(position, direction, state);

    if (inactive()) {
      ACTS_VERBOSE(volInfo(state)
                   << posInfo(state, position) << "navigator inactive");
      return;
    }

    if (state.surfaceCandidateIndex == state.surfaceCandidates.size()) {
      ACTS_VERBOSE(volInfo(state)
                   << posInfo(state, position)
                   << "no surface candidates - waiting for target call");
      return;
    }

    const Portal* nextPortal = nullptr;
    const Surface* nextSurface = nullptr;
    bool isPortal = false;

    if (state.surfaceCandidate().surface != nullptr) {
      nextSurface = state.surfaceCandidate().surface;
    } else if (state.surfaceCandidate().portal != nullptr) {
      nextPortal = state.surfaceCandidate().portal;
      nextSurface = &nextPortal->surface();
      isPortal = true;
    } else {
      std::string msg = "DetectorNavigator: " + volInfo(state) +
                        posInfo(state, position) +
                        "panic: not a surface not a portal - what is it?";
      throw std::runtime_error(msg);
    }

    ACTS_VERBOSE(volInfo(state)
                 << posInfo(state, position) << "landed on surface");

    if (isPortal) {
      ACTS_VERBOSE(volInfo(state)
                   << posInfo(state, position)
                   << "this is a portal, updating to new volume.");
      state.currentPortal = nextPortal;
      state.currentSurface = &nextPortal->surface();
      state.surfaceCandidates.clear();
      state.surfaceCandidateIndex = 0;

      state.currentPortal->updateDetectorVolume(state.options.geoContext,
                                                state);

      // If no Volume is found, we are at the end of the world
      if (state.currentVolume == nullptr) {
        ACTS_VERBOSE(volInfo(state)
                     << posInfo(state, position)
                     << "no volume after Portal update, end of world.");
        state.navigationBreak = true;
        return;
      }

      // Switched to a new volume
      // Update candidate surfaces
      updateCandidateSurfaces(state, position);

      ACTS_VERBOSE(volInfo(state)
                   << posInfo(state, position) << "current portal set to "
                   << state.currentPortal->surface().geometryId());
    } else {
      ACTS_VERBOSE(volInfo(state) << posInfo(state, position)
                                  << "this is a surface, storing it.");

      // If we are on the surface pointed at by the iterator, we can make
      // it the current one to pass it to the other actors
      state.currentSurface = nextSurface;
      ACTS_VERBOSE(volInfo(state)
                   << posInfo(state, position) << "current surface set to "
                   << state.currentSurface->geometryId());
      ++state.surfaceCandidateIndex;
    }
  }

 private:
  Config m_cfg;

  std::shared_ptr<const Logger> m_logger;

  std::string volInfo(const State& state) const {
    return (state.currentVolume != nullptr ? state.currentVolume->name()
                                           : "No Volume") +
           " | ";
  }

  std::string posInfo(const State& /*state*/, const Vector3& position) const {
    std::stringstream ss;
    ss << position.transpose();
    ss << " | ";
    return ss.str();
  }

  const Logger& logger() const { return *m_logger; }

  /// This checks if a navigation break had been triggered or navigator
  /// is misconfigured
  ///
  /// @return true if the navigator is inactive
  bool inactive() const {
    if (m_cfg.detector == nullptr) {
      return true;
    }

    if (!m_cfg.resolveSensitive && !m_cfg.resolveMaterial &&
        !m_cfg.resolvePassive) {
      return true;
    }

    return false;
  }

  /// @brief Navigation (re-)initialisation for the target
  ///
  /// @note This is only called a few times every propagation/extrapolation
  ///
  /// As a straight line estimate can lead you to the wrong destination
  /// Volume, this will be called at:
  /// - initialization
  /// - attempted volume switch
  /// Target finding by association will not be done again
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] position is the current position
  void updateCandidateSurfaces(State& state, const Vector3& position) const {
    ACTS_VERBOSE(volInfo(state)
                 << posInfo(state, position) << "initialize target");

    // Here we get the candidate surfaces
    state.currentVolume->updateNavigationState(state.options.geoContext, state);

    // Sort properly the surface candidates
    auto& nCandidates = state.surfaceCandidates;
    std::ranges::sort(nCandidates, {}, [](const auto& c) {
      return c.objectIntersection.pathLength();
    });
    // Set the surface candidate
    state.surfaceCandidateIndex = 0;
  }

  void fillNavigationState(const Vector3& position, const Vector3& direction,
                           State& state) const {
    state.position = position;
    state.direction = direction;
  }
};

}  // namespace Acts::Experimental
