// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Propagator/NavigatorStatistics.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <algorithm>
#include <optional>
#include <sstream>
#include <string>

#include <boost/container/small_vector.hpp>

namespace Acts {

/// @brief The navigation options for the tracking geometry
///
/// @tparam object_t Type of the object for navigation to check against
template <typename object_t>
struct NavigationOptions {
  /// The boundary check directive
  BoundaryTolerance boundaryTolerance = BoundaryTolerance::None();

  // How to resolve the geometry
  /// Always look for sensitive
  bool resolveSensitive = true;
  /// Always look for material
  bool resolveMaterial = true;
  /// always look for passive
  bool resolvePassive = false;

  /// Hint for start object
  const object_t* startObject = nullptr;
  /// Hint for end object
  const object_t* endObject = nullptr;

  /// External surface identifier for which the boundary check is ignored
  std::vector<GeometryIdentifier> externalSurfaces = {};

  /// The minimum distance for a surface to be considered
  double nearLimit = 0;
  /// The maximum distance for a surface to be considered
  double farLimit = std::numeric_limits<double>::max();
};

/// @brief Steers the propagation through the geometry by providing the next
/// surface to be targeted.
///
/// The Navigator is part of the propagation and responsible for steering
/// the surface sequence to encounter all the relevant surfaces which are
/// intersected by the trajectory.
///
/// The current navigation stage is cached in the state struct and updated
/// when necessary. If any surface in the extrapolation flow is hit, it is
/// set to the navigation state, such that other actors can deal with it.
///
/// The current target surface is referenced by an index which points into
/// the navigation candidates. The navigation candidates are ordered by the
/// path length to the surface. If a surface is hit, the
/// `state.currentSurface` pointer is set. This actors to observe
/// that we are on a surface.
///
class Navigator {
 public:
  using NavigationSurfaces =
      boost::container::small_vector<SurfaceIntersection, 10>;

  using NavigationLayers =
      boost::container::small_vector<LayerIntersection, 10>;

  using NavigationBoundaries =
      boost::container::small_vector<BoundaryIntersection, 4>;

  using ExternalSurfaces = std::multimap<std::uint64_t, GeometryIdentifier>;

  /// The navigation stage
  enum struct Stage : int {
    initial = 0,
    surfaceTarget = 1,
    layerTarget = 2,
    boundaryTarget = 3,
  };

  /// The navigator configuration
  struct Config {
    /// Tracking Geometry for this Navigator
    std::shared_ptr<const TrackingGeometry> trackingGeometry{nullptr};

    /// stop at every sensitive surface (whether it has material or not)
    bool resolveSensitive = true;
    /// stop at every material surface (whether it is passive or not)
    bool resolveMaterial = true;
    /// stop at every surface regardless what it is
    bool resolvePassive = false;
  };

  /// The navigator options
  struct Options : public NavigatorPlainOptions {
    explicit Options(const GeometryContext& gctx)
        : NavigatorPlainOptions(gctx) {}

    /// The surface tolerance
    double surfaceTolerance = s_onSurfaceTolerance;

    /// The near limit to resolve surfaces
    double nearLimit = s_onSurfaceTolerance;

    /// The far limit to resolve surfaces
    double farLimit = std::numeric_limits<double>::max();

    /// Externally provided surfaces - these are tried to be hit
    ExternalSurfaces externalSurfaces = {};

    void insertExternalSurface(GeometryIdentifier geoid) {
      externalSurfaces.insert(
          std::pair<std::uint64_t, GeometryIdentifier>(geoid.layer(), geoid));
    }

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

    // Navigation on surface level
    /// the vector of navigation surfaces to work through
    NavigationSurfaces navSurfaces = {};
    /// the current surface index of the navigation state
    std::optional<std::size_t> navSurfaceIndex;

    // Navigation on layer level
    /// the vector of navigation layers to work through
    NavigationLayers navLayers = {};
    /// the current layer index of the navigation state
    std::optional<std::size_t> navLayerIndex;

    // Navigation on volume level
    /// the vector of boundary surfaces to work through
    NavigationBoundaries navBoundaries = {};
    /// the current boundary index of the navigation state
    std::optional<std::size_t> navBoundaryIndex;

    SurfaceIntersection& navSurface() {
      return navSurfaces.at(navSurfaceIndex.value());
    }
    LayerIntersection& navLayer() {
      return navLayers.at(navLayerIndex.value());
    }
    BoundaryIntersection& navBoundary() {
      return navBoundaries.at(navBoundaryIndex.value());
    }

    const TrackingVolume* startVolume = nullptr;
    const Layer* startLayer = nullptr;
    const Surface* startSurface = nullptr;
    const TrackingVolume* currentVolume = nullptr;
    const Layer* currentLayer = nullptr;
    const Surface* currentSurface = nullptr;
    const Surface* targetSurface = nullptr;

    bool navigationBreak = false;
    Stage navigationStage = Stage::initial;

    NavigatorStatistics statistics;

    void resetAfterLayerSwitch() {
      navSurfaces.clear();
      navSurfaceIndex.reset();
    }

    void resetAfterVolumeSwitch() {
      resetAfterLayerSwitch();

      navLayers.clear();
      navLayerIndex.reset();
      navBoundaries.clear();
      navBoundaryIndex.reset();

      currentLayer = nullptr;
    }

    void reset() {
      resetAfterVolumeSwitch();

      currentVolume = nullptr;
      currentSurface = nullptr;

      navigationBreak = false;
      navigationStage = Stage::initial;
    }
  };

  /// Constructor with configuration object
  ///
  /// @param cfg The navigator configuration
  /// @param _logger a logger instance
  explicit Navigator(Config cfg,
                     std::shared_ptr<const Logger> _logger =
                         getDefaultLogger("Navigator", Logging::Level::INFO))
      : m_cfg{std::move(cfg)}, m_logger{std::move(_logger)} {}

  State makeState(const Options& options) const {
    State state(options);
    state.startSurface = options.startSurface;
    state.targetSurface = options.targetSurface;
    return state;
  }

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

  bool endOfWorldReached(const State& state) const {
    return state.currentVolume == nullptr;
  }

  bool navigationBreak(const State& state) const {
    return state.navigationBreak;
  }

  /// @brief Initialize the navigator
  ///
  /// This function initializes the navigator for a new propagation.
  ///
  /// @param state The navigation state
  /// @param position The start position
  /// @param direction The start direction
  /// @param propagationDirection The propagation direction
  void initialize(State& state, const Vector3& position,
                  const Vector3& direction,
                  Direction propagationDirection) const {
    (void)propagationDirection;

    ACTS_VERBOSE(volInfo(state) << "Initialization.");

    state.reset();

    // Fast Navigation initialization for start condition:
    // - short-cut through object association, saves navigation in the
    // - geometry and volume tree search for the lowest volume
    if (state.startSurface != nullptr &&
        state.startSurface->associatedLayer() != nullptr) {
      ACTS_VERBOSE(
          volInfo(state)
          << "Fast start initialization through association from Surface.");

      state.startLayer = state.startSurface->associatedLayer();
      state.startVolume = state.startLayer->trackingVolume();
    } else if (state.startVolume != nullptr) {
      ACTS_VERBOSE(
          volInfo(state)
          << "Fast start initialization through association from Volume.");

      state.startLayer = state.startVolume->associatedLayer(
          state.options.geoContext, position);
    } else {
      ACTS_VERBOSE(volInfo(state)
                   << "Slow start initialization through search.");
      ACTS_VERBOSE(volInfo(state)
                   << "Starting from position " << toString(position)
                   << " and direction " << toString(direction));

      // current volume and layer search through global search
      state.startVolume = m_cfg.trackingGeometry->lowestTrackingVolume(
          state.options.geoContext, position);

      if (state.startVolume != nullptr) {
        state.startLayer = state.startVolume->associatedLayer(
            state.options.geoContext, position);
      } else {
        ACTS_ERROR(volInfo(state)
                   << "No start volume resolved. Nothing left to do.");
        state.navigationBreak = true;
      }
    }

    state.currentVolume = state.startVolume;
    state.currentLayer = state.startLayer;
    state.currentSurface = state.startSurface;

    if (state.currentVolume != nullptr) {
      ACTS_VERBOSE(volInfo(state) << "Start volume resolved "
                                  << state.currentVolume->geometryId());
      assert(state.currentVolume->inside(position,
                                         state.options.surfaceTolerance) &&
             "We did not end up inside the volume.");
    }
    if (state.currentLayer != nullptr) {
      ACTS_VERBOSE(volInfo(state) << "Start layer resolved "
                                  << state.currentLayer->geometryId());
    }
    if (state.currentSurface != nullptr) {
      ACTS_VERBOSE(volInfo(state) << "Start surface resolved "
                                  << state.currentSurface->geometryId());
      assert(state.currentSurface->isOnSurface(
                 state.options.geoContext, position, direction,
                 BoundaryTolerance::Infinite(),
                 state.options.surfaceTolerance) &&
             "Stepper not on surface");
    }
  }

  /// @brief Get the next target surface
  ///
  /// This function gets the next target surface for the propagation.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  ///
  /// @return The next target surface
  NavigationTarget nextTarget(State& state, const Vector3& position,
                              const Vector3& direction) const {
    if (inactive(state)) {
      return NavigationTarget::None();
    }

    ACTS_VERBOSE(volInfo(state) << "Entering Navigator::nextTarget.");

    // Reset the current surface
    state.currentSurface = nullptr;

    auto tryGetNextTarget = [&]() -> NavigationTarget {
      // Try targeting the surfaces - then layers - then boundaries

      if (state.navigationStage == Stage::initial) {
        ACTS_VERBOSE(volInfo(state) << "Target surfaces.");
        state.navigationStage = Stage::surfaceTarget;
      }

      if (state.navigationStage == Stage::surfaceTarget) {
        if (!state.navSurfaceIndex.has_value()) {
          // First time, resolve the surfaces
          resolveSurfaces(state, position, direction);
          state.navSurfaceIndex = 0;
        } else {
          ++state.navSurfaceIndex.value();
        }
        if (state.navSurfaceIndex.value() < state.navSurfaces.size()) {
          ACTS_VERBOSE(volInfo(state) << "Target set to next surface.");
          return NavigationTarget(*state.navSurface().object(),
                                  state.navSurface().index(),
                                  BoundaryTolerance::None());
        } else {
          // This was the last surface, switch to layers
          ACTS_VERBOSE(volInfo(state) << "Target layers.");
          state.navigationStage = Stage::layerTarget;
        }
      }

      if (state.navigationStage == Stage::layerTarget) {
        if (!state.navLayerIndex.has_value()) {
          // First time, resolve the layers
          resolveLayers(state, position, direction);
          state.navLayerIndex = 0;
        } else {
          ++state.navLayerIndex.value();
        }
        if (state.navLayerIndex.value() < state.navLayers.size()) {
          ACTS_VERBOSE(volInfo(state) << "Target set to next layer.");
          return NavigationTarget(*state.navLayer().first.object(),
                                  state.navLayer().first.index(),
                                  BoundaryTolerance::None());
        } else {
          // This was the last layer, switch to boundaries
          ACTS_VERBOSE(volInfo(state) << "Target boundaries.");
          state.navigationStage = Stage::boundaryTarget;
        }
      }

      if (state.navigationStage == Stage::boundaryTarget) {
        if (!state.navBoundaryIndex.has_value()) {
          // First time, resolve the boundaries
          resolveBoundaries(state, position, direction);
          state.navBoundaryIndex = 0;
        } else {
          ++state.navBoundaryIndex.value();
        }
        if (state.navBoundaryIndex.value() < state.navBoundaries.size()) {
          ACTS_VERBOSE(volInfo(state) << "Target set to next boundary.");
          return NavigationTarget(*state.navBoundary().first.object(),
                                  state.navBoundary().first.index(),
                                  BoundaryTolerance::None());
        } else {
          // This was the last boundary, we have to leave the volume somehow,
          // renavigate
          ACTS_VERBOSE(volInfo(state)
                       << "Boundary targets exhausted. Renavigate.");
        }
      }

      ACTS_VERBOSE(volInfo(state)
                   << "Unknown state. No target found. Renavigate.");
      return NavigationTarget::None();
    };

    NavigationTarget nextTarget = tryGetNextTarget();
    if (!nextTarget.isNone()) {
      return nextTarget;
    }

    state.reset();
    ++state.statistics.nRenavigations;

    // We might have punched through a boundary and entered another volume
    // so we have to reinitialize
    state.currentVolume = m_cfg.trackingGeometry->lowestTrackingVolume(
        state.options.geoContext, position);

    if (state.currentVolume == nullptr) {
      ACTS_VERBOSE(volInfo(state) << "No volume found, stop navigation.");
      state.navigationBreak = true;
      return NavigationTarget::None();
    }

    state.currentLayer = state.currentVolume->associatedLayer(
        state.options.geoContext, position);

    ACTS_VERBOSE(volInfo(state) << "Resolved volume and layer.");

    // Rerun the targeting
    nextTarget = tryGetNextTarget();
    if (!nextTarget.isNone()) {
      return nextTarget;
    }

    ACTS_VERBOSE(
        volInfo(state)
        << "No targets found again, we got really lost! Stop navigation.");
    state.navigationBreak = true;
    return NavigationTarget::None();
  }

  /// @brief Check if the current target is still valid
  ///
  /// This function checks if the target is valid.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  ///
  /// @return True if the target is valid
  bool checkTargetValid(const State& state, const Vector3& position,
                        const Vector3& direction) const {
    (void)position;
    (void)direction;

    return state.navigationStage != Stage::initial;
  }

  /// @brief Handle the surface reached
  ///
  /// This function handles the surface reached.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  /// @param surface The surface reached
  void handleSurfaceReached(State& state, const Vector3& position,
                            const Vector3& direction,
                            const Surface& surface) const {
    if (inactive(state)) {
      return;
    }

    ACTS_VERBOSE(volInfo(state) << "Entering Navigator::handleSurfaceReached.");

    state.currentSurface = &surface;

    if (state.navigationStage == Stage::surfaceTarget &&
        state.navSurface().object() == &surface) {
      ACTS_VERBOSE(volInfo(state) << "Handling surface status.");

      return;
    }

    if (state.navigationStage == Stage::layerTarget &&
        state.navLayer().first.object() == &surface) {
      ACTS_VERBOSE(volInfo(state) << "Handling layer status.");

      // Switch to the next layer
      state.currentLayer = state.navLayer().second;
      state.navigationStage = Stage::surfaceTarget;

      // partial reset
      state.resetAfterLayerSwitch();

      return;
    }

    if (state.navigationStage == Stage::boundaryTarget &&
        state.navBoundary().first.object() == &surface) {
      ACTS_VERBOSE(volInfo(state) << "Handling boundary status.");

      // Switch to the next volume using the boundary
      const BoundarySurface* boundary = state.navBoundary().second;
      state.currentVolume = boundary->attachedVolume(state.options.geoContext,
                                                     position, direction);

      // partial reset
      state.resetAfterVolumeSwitch();

      if (state.currentVolume != nullptr) {
        ACTS_VERBOSE(volInfo(state) << "Volume updated.");
        state.navigationStage = Stage::layerTarget;
      } else {
        ACTS_VERBOSE(volInfo(state)
                     << "No more volume to progress to, stopping navigation.");
        state.navigationBreak = true;
      }

      return;
    }

    ACTS_ERROR(volInfo(state) << "Surface reached but unknown state.");
  }

 private:
  /// @brief Resolve compatible surfaces
  ///
  /// This function resolves the compatible surfaces for the navigation.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  void resolveSurfaces(State& state, const Vector3& position,
                       const Vector3& direction) const {
    ACTS_VERBOSE(volInfo(state) << "Searching for compatible surfaces.");

    const Layer* currentLayer = state.currentLayer;

    if (currentLayer == nullptr) {
      ACTS_VERBOSE(volInfo(state) << "No layer to resolve surfaces.");
      return;
    }

    const Surface* layerSurface = &currentLayer->surfaceRepresentation();

    NavigationOptions<Surface> navOpts;
    navOpts.resolveSensitive = m_cfg.resolveSensitive;
    navOpts.resolveMaterial = m_cfg.resolveMaterial;
    navOpts.resolvePassive = m_cfg.resolvePassive;
    navOpts.startObject = state.currentSurface;
    navOpts.endObject = state.targetSurface;
    navOpts.nearLimit = state.options.nearLimit;
    navOpts.farLimit = state.options.farLimit;

    if (!state.options.externalSurfaces.empty()) {
      auto layerId = layerSurface->geometryId().layer();
      auto externalSurfaceRange =
          state.options.externalSurfaces.equal_range(layerId);
      navOpts.externalSurfaces.reserve(
          state.options.externalSurfaces.count(layerId));
      for (auto itSurface = externalSurfaceRange.first;
           itSurface != externalSurfaceRange.second; itSurface++) {
        navOpts.externalSurfaces.push_back(itSurface->second);
      }
    }

    // Request the compatible surfaces
    state.navSurfaces = currentLayer->compatibleSurfaces(
        state.options.geoContext, position, direction, navOpts);
    // Sort the surfaces by path length.
    // Special care is taken for the external surfaces which should always come
    // first, so they are preferred to be targeted and hit first.
    std::ranges::sort(
        state.navSurfaces, [&state, &navOpts](const SurfaceIntersection& a,
                                              const SurfaceIntersection& b) {
          // Prefer to sort by path length. We assume surfaces are at the same
          // distance if the difference is smaller than the tolerance.
          if (std::abs(a.pathLength() - b.pathLength()) >
              state.options.surfaceTolerance) {
            return SurfaceIntersection::pathLengthOrder(a, b);
          }
          // If the path length is practically the same, sort by geometry. First
          // we check if one of the surfaces is external.
          bool aIsExternal = rangeContainsValue(navOpts.externalSurfaces,
                                                a.object()->geometryId());
          bool bIsExternal = rangeContainsValue(navOpts.externalSurfaces,
                                                b.object()->geometryId());
          if (aIsExternal == bIsExternal) {
            // If both are external or both are not external, sort by geometry
            // identifier
            return a.object()->geometryId() < b.object()->geometryId();
          }
          // If only one is external, it should come first
          return aIsExternal;
        });

    // Print surface information
    if (logger().doPrint(Logging::VERBOSE)) {
      std::ostringstream os;
      os << state.navSurfaces.size();
      os << " surface candidates found at path(s): ";
      for (auto& sfc : state.navSurfaces) {
        os << sfc.pathLength() << "  ";
      }
      logger().log(Logging::VERBOSE, os.str());
    }

    if (state.navSurfaces.empty()) {
      ACTS_VERBOSE(volInfo(state) << "No surface candidates found.");
    }
  }

  /// @brief Resolve compatible layers
  ///
  /// This function resolves the compatible layers for the navigation.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  void resolveLayers(State& state, const Vector3& position,
                     const Vector3& direction) const {
    ACTS_VERBOSE(volInfo(state) << "Searching for compatible layers.");

    NavigationOptions<Layer> navOpts;
    navOpts.resolveSensitive = m_cfg.resolveSensitive;
    navOpts.resolveMaterial = m_cfg.resolveMaterial;
    navOpts.resolvePassive = m_cfg.resolvePassive;
    navOpts.startObject = state.currentLayer;
    navOpts.nearLimit = state.options.nearLimit;
    navOpts.farLimit = state.options.farLimit;

    // Request the compatible layers
    state.navLayers = state.currentVolume->compatibleLayers(
        state.options.geoContext, position, direction, navOpts);
    std::ranges::sort(state.navLayers, [](const auto& a, const auto& b) {
      return SurfaceIntersection::pathLengthOrder(a.first, b.first);
    });

    // Print layer information
    if (logger().doPrint(Logging::VERBOSE)) {
      std::ostringstream os;
      os << state.navLayers.size();
      os << " layer candidates found at path(s): ";
      for (auto& lc : state.navLayers) {
        os << lc.first.pathLength() << "  ";
      }
      logger().log(Logging::VERBOSE, os.str());
    }

    if (state.navLayers.empty()) {
      ACTS_VERBOSE(volInfo(state) << "No layer candidates found.");
    }
  }

  /// @brief Resolve compatible boundaries
  ///
  /// This function resolves the compatible boundaries for the navigation.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  void resolveBoundaries(State& state, const Vector3& position,
                         const Vector3& direction) const {
    ACTS_VERBOSE(volInfo(state) << "Searching for compatible boundaries.");

    NavigationOptions<Surface> navOpts;
    navOpts.startObject = state.currentSurface;
    navOpts.nearLimit = state.options.nearLimit;
    navOpts.farLimit = state.options.farLimit;

    ACTS_VERBOSE(volInfo(state)
                 << "Try to find boundaries, we are at: " << toString(position)
                 << ", dir: " << toString(direction));

    // Request the compatible boundaries
    state.navBoundaries = state.currentVolume->compatibleBoundaries(
        state.options.geoContext, position, direction, navOpts, logger());
    std::ranges::sort(state.navBoundaries, [](const auto& a, const auto& b) {
      return SurfaceIntersection::pathLengthOrder(a.first, b.first);
    });

    // Print boundary information
    if (logger().doPrint(Logging::VERBOSE)) {
      std::ostringstream os;
      os << state.navBoundaries.size();
      os << " boundary candidates found at path(s): ";
      for (auto& bc : state.navBoundaries) {
        os << bc.first.pathLength() << "  ";
      }
      logger().log(Logging::VERBOSE, os.str());
    }

    if (state.navBoundaries.empty()) {
      ACTS_VERBOSE(volInfo(state) << "No boundary candidates found.");
    }
  }

  /// @brief Check if the navigator is inactive
  ///
  /// This function checks if the navigator is inactive.
  ///
  /// @param state The navigation state
  ///
  /// @return True if the navigator is inactive
  bool inactive(const State& state) const {
    // Void behavior in case no tracking geometry is present
    if (m_cfg.trackingGeometry == nullptr) {
      return true;
    }

    // Turn the navigator into void when you are instructed to do nothing
    if (!m_cfg.resolveSensitive && !m_cfg.resolveMaterial &&
        !m_cfg.resolvePassive) {
      return true;
    }

    if (state.navigationBreak) {
      return true;
    }

    return false;
  }

 private:
  template <typename propagator_state_t>
  std::string volInfo(const propagator_state_t& state) const {
    return (state.currentVolume != nullptr ? state.currentVolume->volumeName()
                                           : "No Volume") +
           " | ";
  }

  const Logger& logger() const { return *m_logger; }

  Config m_cfg;

  std::shared_ptr<const Logger> m_logger;
};

}  // namespace Acts
