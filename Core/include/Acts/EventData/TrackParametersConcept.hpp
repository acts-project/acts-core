// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

#include <optional>

namespace Acts {

class Surface;

namespace Concepts {
template <typename Parameters>
concept BasicTrackParameters = requires {
  typename Parameters::ParametersVector;
  typename Parameters::CovarianceMatrix;

  requires requires(const Parameters &p) {
    { p.time() } -> std::floating_point;
    { p.direction() } -> std::same_as<Vector3>;
    { p.absoluteMomentum() } -> std::floating_point;
    { p.charge() } -> std::floating_point;
  };
};
}  // namespace Concepts

/// @brief Concept that asserts that a given type meets the requirements to be
/// considered free track parameters in Acts.
template <typename Parameters>
concept FreeTrackParametersConcept =
    Concepts::BasicTrackParameters<Parameters> &&
    requires(const Parameters &p) {
      { p.parameters() } -> std::convertible_to<FreeVector>;
      {
        p.covariance()
      } -> std::convertible_to<std::optional<FreeSquareMatrix>>;
      { p.fourPosition() } -> std::same_as<Vector4>;
      { p.position() } -> std::same_as<Vector3>;
    };

/// @brief Concept that asserts that a given type meets the requirements to be
/// considered bound track parameters in Acts.
template <typename Parameters>
concept BoundTrackParametersConcept =
    Concepts::BasicTrackParameters<Parameters> &&
    requires(const Parameters &p) {
      { p.parameters() } -> std::convertible_to<BoundVector>;
      {
        p.covariance()
      } -> std::convertible_to<std::optional<BoundSquareMatrix>>;
      { p.referenceSurface() } -> std::same_as<const Surface &>;

      requires requires(GeometryContext &c) {
        { p.position(c) } -> std::same_as<Vector3>;
        { p.fourPosition(c) } -> std::same_as<Vector4>;
        { p.position(c) } -> std::same_as<Vector3>;
      };
    };

namespace Concepts {
template <typename Parameters>
concept BoundConvertibleTrackParameters = requires(const Parameters &p) {
  { p.toBound() } -> BoundTrackParametersConcept;
};
}  // namespace Concepts

}  // namespace Acts
