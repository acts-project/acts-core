// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Propagator/NavigatorStatistics.hpp"
#include "Acts/Propagator/detail/NavigationHelpers.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

namespace Acts {

/// @brief Captures the common functionality of the try-all navigators
///
/// This class is not meant to be used directly, but to be inherited by the
/// actual navigator implementations.
///
class TryAllNavigatorBase {
 public:
  /// @brief Configuration for this Navigator
  struct Config {
    /// Tracking Geometry for this Navigator
    std::shared_ptr<const TrackingGeometry> trackingGeometry;

    /// stop at every sensitive surface (whether it has material or not)
    bool resolveSensitive = true;
    /// stop at every material surface (whether it is passive or not)
    bool resolveMaterial = true;
    /// stop at every surface regardless what it is
    bool resolvePassive = false;

    /// Which boundary checks to perform for surface approach
    BoundaryTolerance boundaryToleranceSurfaceApproach =
        BoundaryTolerance::None();
  };

  /// @brief Options for this Navigator
  struct Options : public NavigatorPlainOptions {
    explicit Options(const GeometryContext& gctx)
        : NavigatorPlainOptions(gctx) {}

    /// The surface tolerance
    double surfaceTolerance = s_onSurfaceTolerance;

    /// The near limit to resolve surfaces
    double nearLimit = s_onSurfaceTolerance;

    /// The far limit to resolve surfaces
    double farLimit = std::numeric_limits<double>::max();

    void setPlainOptions(const NavigatorPlainOptions& options) {
      static_cast<NavigatorPlainOptions&>(*this) = options;
    }
  };

  /// @brief Nested State struct
  ///
  /// It acts as an internal state which is created for every propagation and
  /// meant to keep thread-local navigation information.
  struct State {
    explicit State(const Options& options_) : options(options_) {}

    Options options;

    // Starting geometry information of the navigation which should only be set
    // while initialization. NOTE: This information is mostly used by actors to
    // check if we are on the starting surface (e.g. MaterialInteraction).
    const Surface* startSurface = nullptr;

    // Target geometry information of the navigation which should only be set
    // while initialization. NOTE: This information is mostly used by actors to
    // check if we are on the target surface (e.g. MaterialInteraction).
    const Surface* targetSurface = nullptr;

    // Current geometry information of the navigation which is set during
    // initialization and potentially updated after each step.
    const Surface* currentSurface = nullptr;
    const TrackingVolume* currentVolume = nullptr;

    /// The vector of navigation candidates to work through
    std::vector<detail::NavigationObjectCandidate> navigationCandidates;

    /// If a break has been detected
    bool navigationBreak = false;

    /// Navigation statistics
    NavigatorStatistics statistics;
  };

  /// Constructor with configuration object
  ///
  /// @param cfg The navigator configuration
  /// @param logger a logger instance
  TryAllNavigatorBase(Config cfg, std::unique_ptr<const Logger> logger)
      : m_cfg(std::move(cfg)), m_logger{std::move(logger)} {}

  const Surface* currentSurface(const State& state) const {
    return state.currentSurface;
  }

  const TrackingVolume* currentVolume(const State& state) const {
    return state.currentVolume;
  }

  const IVolumeMaterial* currentVolumeMaterial(const State& state) const {
    if (state.currentVolume == nullptr) {
      return nullptr;
    }
    return state.currentVolume->volumeMaterial();
  }

  const Surface* startSurface(const State& state) const {
    return state.startSurface;
  }

  const Surface* targetSurface(const State& state) const {
    return state.targetSurface;
  }

  bool endOfWorldReached(State& state) const {
    return state.currentVolume == nullptr;
  }

  bool navigationBreak(const State& state) const {
    return state.navigationBreak;
  }

  /// @brief Initialize the navigator
  ///
  /// This method initializes the navigator for a new propagation. It sets the
  /// current volume and surface to the start volume and surface, respectively.
  ///
  /// @param state The navigation state
  /// @param position The starting position
  /// @param direction The starting direction
  /// @param propagationDirection The propagation direction
  Result<void> initialize(State& state, const Vector3& position,
                          const Vector3& direction,
                          Direction propagationDirection) const {
    (void)propagationDirection;

    ACTS_VERBOSE("initialize");

    const TrackingVolume* startVolume = nullptr;

    if (state.startSurface != nullptr &&
        state.startSurface->associatedLayer() != nullptr) {
      ACTS_VERBOSE(
          "Fast start initialization through association from Surface.");
      const auto* startLayer = state.startSurface->associatedLayer();
      startVolume = startLayer->trackingVolume();
    } else {
      ACTS_VERBOSE("Slow start initialization through search.");
      ACTS_VERBOSE("Starting from position " << toString(position)
                                             << " and direction "
                                             << toString(direction));
      startVolume = m_cfg.trackingGeometry->lowestTrackingVolume(
          state.options.geoContext, position);
    }

    // Initialize current volume, layer and surface
    {
      state.currentVolume = startVolume;
      if (state.currentVolume != nullptr) {
        ACTS_VERBOSE(volInfo(state) << "Start volume resolved.");
      } else {
        ACTS_ERROR("Start volume not resolved.");
      }

      state.currentSurface = state.startSurface;
      if (state.currentSurface != nullptr) {
        ACTS_VERBOSE(volInfo(state) << "Current surface set to start surface "
                                    << state.currentSurface->geometryId());
      } else {
        ACTS_VERBOSE(volInfo(state) << "No start surface set.");
      }
    }

    return Result<void>::success();
  }

 protected:
  /// Helper method to initialize navigation candidates for the current volume.
  void initializeVolumeCandidates(State& state) const {
    const TrackingVolume* volume = state.currentVolume;
    ACTS_VERBOSE(volInfo(state) << "Initialize volume");

    if (volume == nullptr) {
      state.navigationBreak = true;
      ACTS_VERBOSE(volInfo(state) << "No volume set. Good luck.");
      return;
    }

    emplaceAllVolumeCandidates(
        state.navigationCandidates, *volume, m_cfg.resolveSensitive,
        m_cfg.resolveMaterial, m_cfg.resolvePassive,
        m_cfg.boundaryToleranceSurfaceApproach, logger());
  }

  std::string volInfo(const State& state) const {
    return (state.currentVolume != nullptr ? state.currentVolume->volumeName()
                                           : "No Volume") +
           " | ";
  }

  const Logger& logger() const { return *m_logger; }

  Config m_cfg;

  std::unique_ptr<const Logger> m_logger;
};

/// @brief Alternative @c Navigator which tries all possible intersections
///
/// See @c Navigator for more general information about the Navigator concept.
///
/// This Navigator tries all possible intersections with all surfaces in the
/// current volume. It does not use any information about the geometry to
/// optimize the search. It is therefore very slow, but can be used as a
/// reference implementation.
///
class TryAllNavigator : public TryAllNavigatorBase {
 public:
  using Config = TryAllNavigatorBase::Config;
  using Options = TryAllNavigatorBase::Options;

  /// @brief Nested State struct
  struct State : public TryAllNavigatorBase::State {
    explicit State(const Options& options_)
        : TryAllNavigatorBase::State(options_) {}

    std::vector<detail::IntersectedNavigationObject> currentCandidates;
  };

  /// Constructor with configuration object
  ///
  /// @param cfg The navigator configuration
  /// @param logger a logger instance
  TryAllNavigator(Config cfg,
                  std::unique_ptr<const Logger> logger =
                      getDefaultLogger("TryAllNavigator", Logging::INFO))
      : TryAllNavigatorBase(std::move(cfg), std::move(logger)) {}

  State makeState(const Options& options) const {
    State state(options);
    state.startSurface = options.startSurface;
    state.targetSurface = options.targetSurface;

    return state;
  }

  using TryAllNavigatorBase::currentSurface;
  using TryAllNavigatorBase::currentVolume;
  using TryAllNavigatorBase::currentVolumeMaterial;
  using TryAllNavigatorBase::endOfWorldReached;
  using TryAllNavigatorBase::navigationBreak;
  using TryAllNavigatorBase::startSurface;
  using TryAllNavigatorBase::targetSurface;

  /// @brief Initialize the navigator
  ///
  /// This method initializes the navigator for a new propagation. It sets the
  /// current volume and surface to the start volume and surface, respectively.
  ///
  /// @param state The navigation state
  /// @param position The starting position
  /// @param direction The starting direction
  /// @param propagationDirection The propagation direction
  Result<void> initialize(State& state, const Vector3& position,
                          const Vector3& direction,
                          Direction propagationDirection) const {
    auto baseRes = TryAllNavigatorBase::initialize(state, position, direction,
                                                   propagationDirection);
    if (!baseRes.ok()) {
      return baseRes.error();
    }

    // Initialize navigation candidates for the start volume
    reinitializeCandidates(state);

    return Result<void>::success();
  }

  /// @brief Get the next target surface
  ///
  /// This method gets the next target surface based on the current
  /// position and direction. It returns a none target if no target can be
  /// found.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  ///
  /// @return The next target surface
  NavigationTarget nextTarget(State& state, const Vector3& position,
                              const Vector3& direction) const {
    // Check if the navigator is inactive
    if (state.navigationBreak) {
      return NavigationTarget::None();
    }

    ACTS_VERBOSE(volInfo(state) << "nextTarget");

    // Navigator preStep always resets the current surface
    state.currentSurface = nullptr;

    double nearLimit = state.options.nearLimit;
    double farLimit = state.options.farLimit;

    // handle overstepping
    if (!state.currentCandidates.empty()) {
      const detail::IntersectedNavigationObject& previousCandidate =
          state.currentCandidates.front();

      const Surface& surface = *previousCandidate.intersection.object();
      std::uint8_t index = previousCandidate.intersection.index();
      BoundaryTolerance boundaryTolerance = previousCandidate.boundaryTolerance;

      auto intersection = surface.intersect(
          state.options.geoContext, position, direction, boundaryTolerance,
          state.options.surfaceTolerance)[index];

      if (intersection.pathLength() < 0) {
        nearLimit = std::min(nearLimit, intersection.pathLength() -
                                            state.options.surfaceTolerance);
        farLimit = -state.options.surfaceTolerance;

        ACTS_VERBOSE(volInfo(state)
                     << "handle overstepping with nearLimit " << nearLimit
                     << " and farLimit " << farLimit);
      }
    }

    std::vector<detail::IntersectedNavigationObject> intersectionCandidates;

    // Find intersections with all candidates
    for (const auto& candidate : state.navigationCandidates) {
      auto intersections =
          candidate.intersect(state.options.geoContext, position, direction,
                              state.options.surfaceTolerance);
      for (const auto& intersection : intersections.first.split()) {
        // exclude invalid intersections
        if (!intersection.isValid() ||
            !detail::checkPathLength(intersection.pathLength(), nearLimit,
                                     farLimit)) {
          continue;
        }
        // store candidate
        intersectionCandidates.emplace_back(intersection, intersections.second,
                                            candidate.boundaryTolerance);
      }
    }

    std::ranges::sort(intersectionCandidates,
                      detail::IntersectedNavigationObject::forwardOrder);

    ACTS_VERBOSE(volInfo(state) << "found " << intersectionCandidates.size()
                                << " intersections");

    NavigationTarget nextTarget = NavigationTarget::None();
    state.currentCandidates.clear();

    for (const auto& candidate : intersectionCandidates) {
      const auto& intersection = candidate.intersection;
      const Surface& surface = *intersection.object();
      BoundaryTolerance boundaryTolerance = candidate.boundaryTolerance;

      if (intersection.status() == IntersectionStatus::onSurface) {
        ACTS_ERROR(volInfo(state)
                   << "We are on surface " << surface.geometryId()
                   << " before trying to reach it. This should not happen. "
                      "Good luck.");
        continue;
      }

      if (intersection.status() == IntersectionStatus::reachable) {
        nextTarget =
            NavigationTarget(surface, intersection.index(), boundaryTolerance);
        break;
      }
    }

    state.currentCandidates = std::move(intersectionCandidates);

    if (nextTarget.isNone()) {
      ACTS_VERBOSE(volInfo(state) << "no target found");
    } else {
      ACTS_VERBOSE(volInfo(state)
                   << "next target is " << nextTarget.surface->geometryId());
    }

    return nextTarget;
  }

  /// @brief Check if the target is still valid
  ///
  /// This method checks if the target is valid based on the current position
  /// and direction. It returns true if the target is still valid.
  ///
  /// For the TryAllNavigator, the target is always invalid since we do not want
  /// to assume any specific surface sequence over multiple steps.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  ///
  /// @return True if the target is still valid
  bool checkTargetValid(const State& state, const Vector3& position,
                        const Vector3& direction) const {
    (void)state;
    (void)position;
    (void)direction;

    return false;
  }

  /// @brief Handle the surface reached
  ///
  /// This method is called when a surface is reached. It sets the current
  /// surface in the navigation state and updates the navigation candidates.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  void handleSurfaceReached(State& state, const Vector3& position,
                            const Vector3& direction,
                            const Surface& /*surface*/) const {
    // Check if the navigator is inactive
    if (state.navigationBreak) {
      return;
    }

    ACTS_VERBOSE(volInfo(state) << "handleSurfaceReached");

    if (state.currentCandidates.empty()) {
      ACTS_VERBOSE(volInfo(state) << "No current candidate set.");
      return;
    }

    assert(state.currentSurface == nullptr && "Current surface must be reset.");

    // handle multiple surface intersections due to increased bounds

    std::vector<detail::IntersectedNavigationObject> hitCandidates;

    for (const auto& candidate : state.currentCandidates) {
      const Surface& surface = *candidate.intersection.object();
      std::uint8_t index = candidate.intersection.index();
      BoundaryTolerance boundaryTolerance = BoundaryTolerance::None();

      auto intersection = surface.intersect(
          state.options.geoContext, position, direction, boundaryTolerance,
          state.options.surfaceTolerance)[index];

      if (intersection.status() == IntersectionStatus::onSurface) {
        hitCandidates.emplace_back(candidate);
      }
    }

    state.currentCandidates.clear();

    ACTS_VERBOSE(volInfo(state)
                 << "Found " << hitCandidates.size()
                 << " intersections on surface with bounds check.");

    if (hitCandidates.empty()) {
      ACTS_VERBOSE(volInfo(state) << "No hit candidates found.");
      return;
    }

    if (hitCandidates.size() > 1) {
      ACTS_VERBOSE(volInfo(state)
                   << "Only using first intersection within bounds.");
    }

    // we can only handle a single surface hit so we pick the first one
    const auto candidate = hitCandidates.front();
    const auto& intersection = candidate.intersection;
    const Surface& surface = *intersection.object();

    ACTS_VERBOSE(volInfo(state) << "Surface " << surface.geometryId()
                                << " successfully hit, storing it.");
    state.currentSurface = &surface;

    if (candidate.template checkType<Surface>()) {
      ACTS_VERBOSE(volInfo(state) << "This is a surface");
    } else if (candidate.template checkType<Layer>()) {
      ACTS_VERBOSE(volInfo(state) << "This is a layer");
    } else if (candidate.template checkType<BoundarySurface>()) {
      ACTS_VERBOSE(volInfo(state)
                   << "This is a boundary. Reinitialize navigation");

      const auto& boundary = *candidate.template object<BoundarySurface>();

      state.currentVolume = boundary.attachedVolume(state.options.geoContext,
                                                    position, direction);

      ACTS_VERBOSE(volInfo(state) << "Switched volume");

      reinitializeCandidates(state);
    } else {
      ACTS_ERROR(volInfo(state) << "Unknown intersection type");
    }
  }

 private:
  /// Helper method to reset and reinitialize the navigation candidates.
  void reinitializeCandidates(State& state) const {
    state.navigationCandidates.clear();
    state.currentCandidates.clear();

    initializeVolumeCandidates(state);
  }
};

/// @brief Alternative @c Navigator which tries all possible intersections
///
/// See @c Navigator for more general information about the Navigator concept.
///
/// This Navigator tries all possible intersections with all surfaces in the
/// current volume. It does not use any information about the geometry to
/// optimize the search. It is therefore very slow, but can be used as a
/// reference implementation.
///
/// Different to @c TryAllNavigator, this Navigator discovers intersections by
/// stepping forward blindly and then checking for intersections with all
/// surfaces in the current volume. This is slower, but more robust against
/// bent tracks.
///
class TryAllOverstepNavigator : public TryAllNavigatorBase {
 public:
  using Config = TryAllNavigatorBase::Config;

  using Options = TryAllNavigatorBase::Options;

  /// @brief Nested State struct
  ///
  /// It acts as an internal state which is created for every propagation and
  /// meant to keep thread-local navigation information.
  struct State : public TryAllNavigatorBase::State {
    explicit State(const Options& options_)
        : TryAllNavigatorBase::State(options_) {}

    /// The vector of active intersection candidates to work through
    std::vector<detail::IntersectedNavigationObject> activeCandidates;
    /// The current active candidate index of the navigation state
    int activeCandidateIndex = -1;

    /// The position before the last step
    std::optional<Vector3> lastPosition;

    /// Provides easy access to the active intersection candidate
    const detail::IntersectedNavigationObject& activeCandidate() const {
      return activeCandidates.at(activeCandidateIndex);
    }

    bool endOfCandidates() const {
      return activeCandidateIndex >= static_cast<int>(activeCandidates.size());
    }
  };

  /// Constructor with configuration object
  ///
  /// @param cfg The navigator configuration
  /// @param logger a logger instance
  TryAllOverstepNavigator(Config cfg,
                          std::unique_ptr<const Logger> logger =
                              getDefaultLogger("TryAllOverstepNavigator",
                                               Logging::INFO))
      : TryAllNavigatorBase(std::move(cfg), std::move(logger)) {}

  State makeState(const Options& options) const {
    State state(options);
    state.startSurface = options.startSurface;
    state.targetSurface = options.targetSurface;

    return state;
  }

  using TryAllNavigatorBase::currentSurface;
  using TryAllNavigatorBase::currentVolume;
  using TryAllNavigatorBase::currentVolumeMaterial;
  using TryAllNavigatorBase::endOfWorldReached;
  using TryAllNavigatorBase::navigationBreak;
  using TryAllNavigatorBase::startSurface;
  using TryAllNavigatorBase::targetSurface;

  /// @brief Initialize the navigator
  ///
  /// This method initializes the navigator for a new propagation. It sets the
  /// current volume and surface to the start volume and surface, respectively.
  ///
  /// @param state The navigation state
  /// @param position The starting position
  /// @param direction The starting direction
  /// @param propagationDirection The propagation direction
  Result<void> initialize(State& state, const Vector3& position,
                          const Vector3& direction,
                          Direction propagationDirection) const {
    auto baseRes = TryAllNavigatorBase::initialize(state, position, direction,
                                                   propagationDirection);
    if (!baseRes.ok()) {
      return baseRes.error();
    }

    // Initialize navigation candidates for the start volume
    reinitializeCandidates(state);

    state.lastPosition.reset();

    return Result<void>::success();
  }

  /// @brief Get the next target surface
  ///
  /// This method gets the next target surface based on the current
  /// position and direction. It returns an invalid target if no target can be
  /// found.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  ///
  /// @return The next target surface
  NavigationTarget nextTarget(State& state, const Vector3& position,
                              const Vector3& direction) const {
    (void)direction;

    if (state.navigationBreak) {
      return NavigationTarget::None();
    }

    ACTS_VERBOSE(volInfo(state) << "nextTarget");

    state.currentSurface = nullptr;

    // We cannot do anything without a last position
    if (!state.lastPosition.has_value() && state.endOfCandidates()) {
      ACTS_VERBOSE(
          volInfo(state)
          << "Initial position, nothing to do, blindly stepping forward.");
      state.lastPosition = position;
      return NavigationTarget::None();
    }

    if (state.endOfCandidates()) {
      ACTS_VERBOSE(volInfo(state) << "evaluate blind step");

      Vector3 stepStart = state.lastPosition.value();
      Vector3 stepEnd = position;
      Vector3 step = stepEnd - stepStart;
      double stepDistance = step.norm();
      if (stepDistance < std::numeric_limits<double>::epsilon()) {
        ACTS_ERROR(volInfo(state) << "Step distance is zero. " << stepDistance);
      }
      Vector3 stepDirection = step.normalized();

      double nearLimit = -stepDistance + state.options.surfaceTolerance;
      double farLimit = 0;

      state.lastPosition.reset();
      state.activeCandidates.clear();
      state.activeCandidateIndex = -1;

      // Find intersections with all candidates
      for (const auto& candidate : state.navigationCandidates) {
        auto intersections =
            candidate.intersect(state.options.geoContext, stepEnd,
                                stepDirection, state.options.surfaceTolerance);
        for (const auto& intersection : intersections.first.split()) {
          // exclude invalid intersections
          if (!intersection.isValid() ||
              !detail::checkPathLength(intersection.pathLength(), nearLimit,
                                       farLimit)) {
            continue;
          }
          // store candidate
          state.activeCandidates.emplace_back(
              intersection, intersections.second, candidate.boundaryTolerance);
        }
      }

      std::ranges::sort(state.activeCandidates,
                        detail::IntersectedNavigationObject::forwardOrder);

      ACTS_VERBOSE(volInfo(state) << "Found " << state.activeCandidates.size()
                                  << " intersections");

      for (const auto& candidate : state.activeCandidates) {
        ACTS_VERBOSE("found candidate "
                     << candidate.intersection.object()->geometryId());
      }
    }

    ++state.activeCandidateIndex;

    if (state.endOfCandidates()) {
      ACTS_VERBOSE(volInfo(state)
                   << "No target found, blindly stepping forward.");
      state.lastPosition = position;
      return NavigationTarget::None();
    }

    ACTS_VERBOSE(volInfo(state) << "handle active candidates");

    ACTS_VERBOSE(volInfo(state)
                 << (state.activeCandidates.size() - state.activeCandidateIndex)
                 << " out of " << state.activeCandidates.size()
                 << " surfaces remain to try.");

    const auto& candidate = state.activeCandidate();
    const auto& intersection = candidate.intersection;
    const Surface& surface = *intersection.object();
    BoundaryTolerance boundaryTolerance = candidate.boundaryTolerance;

    ACTS_VERBOSE(volInfo(state)
                 << "Next surface candidate will be " << surface.geometryId());

    return NavigationTarget(surface, intersection.index(), boundaryTolerance);
  }

  /// @brief Check if the target is still valid
  ///
  /// This method checks if the target is valid based on the current position
  /// and direction. It returns true if the target is still valid.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  ///
  /// @return True if the target is still valid
  bool checkTargetValid(const State& state, const Vector3& position,
                        const Vector3& direction) const {
    (void)state;
    (void)position;
    (void)direction;

    return true;
  }

  /// @brief Handle the surface reached
  ///
  /// This method is called when a surface is reached. It sets the current
  /// surface in the navigation state and updates the navigation candidates.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  void handleSurfaceReached(State& state, const Vector3& position,
                            const Vector3& direction,
                            const Surface& /*surface*/) const {
    if (state.navigationBreak) {
      return;
    }

    ACTS_VERBOSE(volInfo(state) << "handleSurfaceReached");

    assert(state.currentSurface == nullptr && "Current surface must be reset.");

    if (state.endOfCandidates()) {
      ACTS_VERBOSE(volInfo(state) << "No active candidate set.");
      return;
    }

    std::vector<detail::IntersectedNavigationObject> hitCandidates;

    while (!state.endOfCandidates()) {
      const auto& candidate = state.activeCandidate();
      const auto& intersection = candidate.intersection;
      const Surface& surface = *intersection.object();
      BoundaryTolerance boundaryTolerance = candidate.boundaryTolerance;

      // first with boundary tolerance
      IntersectionStatus surfaceStatus =
          surface
              .intersect(state.options.geoContext, position, direction,
                         boundaryTolerance,
                         state.options.surfaceTolerance)[intersection.index()]
              .status();

      if (surfaceStatus != IntersectionStatus::onSurface) {
        break;
      }

      // now without boundary tolerance
      boundaryTolerance = BoundaryTolerance::None();
      surfaceStatus =
          surface
              .intersect(state.options.geoContext, position, direction,
                         boundaryTolerance,
                         state.options.surfaceTolerance)[intersection.index()]
              .status();

      if (surfaceStatus == IntersectionStatus::onSurface) {
        hitCandidates.emplace_back(candidate);
      }

      ++state.activeCandidateIndex;
      ACTS_VERBOSE("skip candidate " << surface.geometryId());
    }

    // we increased the candidate index one too many times
    --state.activeCandidateIndex;

    ACTS_VERBOSE(volInfo(state)
                 << "Found " << hitCandidates.size()
                 << " intersections on surface with bounds check.");

    if (hitCandidates.empty()) {
      ACTS_VERBOSE(volInfo(state)
                   << "Surface successfully hit, but outside bounds.");
      return;
    }

    if (hitCandidates.size() > 1) {
      ACTS_VERBOSE(volInfo(state)
                   << "Only using first intersection within bounds.");
    }

    // we can only handle a single surface hit so we pick the first one
    const auto& candidate = hitCandidates.front();
    const auto& intersection = candidate.intersection;
    const Surface& surface = *intersection.object();

    ACTS_VERBOSE(volInfo(state) << "Surface successfully hit, storing it.");
    state.currentSurface = &surface;

    if (state.currentSurface != nullptr) {
      ACTS_VERBOSE(volInfo(state) << "Current surface set to surface "
                                  << surface.geometryId());
    }

    if (candidate.template checkType<Surface>()) {
      ACTS_VERBOSE(volInfo(state) << "This is a surface");
    } else if (candidate.template checkType<Layer>()) {
      ACTS_VERBOSE(volInfo(state) << "This is a layer");
    } else if (candidate.template checkType<BoundarySurface>()) {
      ACTS_VERBOSE(volInfo(state)
                   << "This is a boundary. Reinitialize navigation");

      const auto& boundary = *candidate.template object<BoundarySurface>();

      state.currentVolume = boundary.attachedVolume(state.options.geoContext,
                                                    position, direction);

      ACTS_VERBOSE(volInfo(state) << "Switched volume");

      reinitializeCandidates(state);
    } else {
      ACTS_ERROR(volInfo(state) << "Unknown intersection type");
    }
  }

 private:
  /// Helper method to reset and reinitialize the navigation candidates.
  void reinitializeCandidates(State& state) const {
    state.navigationCandidates.clear();
    state.activeCandidates.clear();
    state.activeCandidateIndex = -1;

    initializeVolumeCandidates(state);
  }
};

}  // namespace Acts
