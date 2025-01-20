// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Propagator/NavigatorStatistics.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {

class Surface;

/// @brief A navigator that does nothing
///
/// It does not provide any navigation action
///
class VoidNavigator {
 public:
  /// @brief Nested Config struct
  struct Config {};

  /// @brief Nested Options struct
  struct Options : public NavigatorPlainOptions {
    explicit Options(const GeometryContext& gctx)
        : NavigatorPlainOptions(gctx) {}

    void setPlainOptions(const NavigatorPlainOptions& options) {
      static_cast<NavigatorPlainOptions&>(*this) = options;
    }
  };

  /// @brief Nested State struct
  struct State {
    explicit State(const Options& options_) : options(options_) {}

    Options options;

    /// Navigation statistics
    NavigatorStatistics statistics;
  };

  State makeState(const Options& options) const {
    State state(options);
    return state;
  }

  const Surface* currentSurface(const State& /*state*/) const {
    return nullptr;
  }

  const Surface* startSurface(const State& /*state*/) const { return nullptr; }

  const Surface* targetSurface(const State& /*state*/) const { return nullptr; }

  bool navigationBreak(const State& /*state*/) const { return true; }

  Result<void> initialize(State& /*state*/, const Vector3& /*position*/,
                          const Vector3& /*direction*/,
                          Direction /*propagationDirection*/) const {
    return Result<void>::success();
  }

  NavigationTarget nextTarget(State& /*state*/, const Vector3& /*position*/,
                              const Vector3& /*direction*/) const {
    return NavigationTarget::None();
  }

  bool checkTargetValid(const State& /*state*/, const Vector3& /*position*/,
                        const Vector3& /*direction*/) const {
    return true;
  }

  void handleSurfaceReached(State& /*state*/, const Vector3& /*position*/,
                            const Vector3& /*direction*/,
                            const Surface& /*surface*/) const {
    return;
  }
};

}  // namespace Acts
