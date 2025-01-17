// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <vector>

namespace Acts {

/// @class CylinderBounds
/// @image html CylinderBounds.gif
///
/// Bounds for a cylindrical Surface.
///
/// These bounds may be used for a CylinderSurface
/// In case of bounds for a StraightLineSurface the radius determines the radius
/// within a localPosition
/// is regarded as inside bounds.
///
/// CylinderBounds also enhance the possibility of a cylinder segment with an
/// opening angle @f$ 2\cdot\phi_{half}@f$
/// around an average @f$ \phi @f$ angle @f$ \phi_{ave} @f$.
///
/// CylinderBounds also supports beveled sides defined by an angle.
/// Different angles can be defined on both sides of the cylinder.
/// A positive angle is defined as "extruding" from the defined Zlength,
/// while a negative angle is "intruding" on the Zlength.
/// +    -            -   +
/// ________________________
/// \  |  /          \  |  /
///  \ | /            \ | /
///   \|/______________\|/
///     2 * ZhalfLength
///
class CylinderBounds : public SurfaceBounds {
 public:
  enum BoundValues : int {
    eR = 0,
    eHalfLengthZ = 1,
    eHalfPhiSector = 2,
    eAveragePhi = 3,
    eBevelMinZ = 4,
    eBevelMaxZ = 5,
    eSize = 6
  };

  /// Constructor - full cylinder
  ///
  /// @param r The radius of the cylinder
  /// @param halfZ The half length in z
  /// @param halfPhi The half opening angle
  /// @param avgPhi (optional) The phi value from which the opening angle spans
  /// @param bevelMinZ (optional) The bevel on the negative z side
  /// @param bevelMaxZ (optional) The bevel on the positive z sid The bevel on the positive z side
  CylinderBounds(double r, double halfZ, double halfPhi = std::numbers::pi,
                 double avgPhi = 0., double bevelMinZ = 0.,
                 double bevelMaxZ = 0.) noexcept(false)
      : m_values({r, halfZ, halfPhi, avgPhi, bevelMinZ, bevelMaxZ}),
        m_closed(std::abs(halfPhi - std::numbers::pi) < s_epsilon) {
    checkConsistency();
  }

  /// Constructor - from fixed size array
  ///
  /// @param values The parameter values
  CylinderBounds(const std::array<double, eSize>& values) noexcept(false)
      : m_values(values),
        m_closed(std::abs(values[eHalfPhiSector] - std::numbers::pi) <
                 s_epsilon) {
    checkConsistency();
  }

  /// @copydoc SurfaceBounds::type
  BoundsType type() const final { return eCylinder; }

  /// @copydoc SurfaceBounds::isCartesian
  bool isCartesian() const final { return true; }

  /// @copydoc SurfaceBounds::boundToCartesianJacobian
  SquareMatrix2 boundToCartesianJacobian(const Vector2& lposition) const final {
    (void)lposition;
    return SquareMatrix2::Identity();
  }

  /// @copydoc SurfaceBounds::boundToCartesianMetric
  SquareMatrix2 boundToCartesianMetric(const Vector2& lposition) const final {
    (void)lposition;
    return SquareMatrix2::Identity();
  }

  /// Return the bound values as dynamically sized vector
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

  /// @copydoc SurfaceBounds::inside
  bool inside(const Vector2& lposition) const final;

  /// @copydoc SurfaceBounds::closestPoint
  Vector2 closestPoint(const Vector2& lposition,
                       const std::optional<SquareMatrix2>& metric) const final;

  using SurfaceBounds::inside;

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  double get(BoundValues bValue) const { return m_values[bValue]; }

  /// Returns true for full phi coverage
  bool coversFullAzimuth() const { return m_closed; }

  /// Create the bow/circle vertices on either side of the cylinder
  ///
  /// @param transform is the global transform
  /// @param quarterSegments is the number of segments to approximate a quarter
  /// of a circle. In order to symmetrize fully closed and sectoral cylinders,
  /// also in the first case the two end points are given (albeit they overlap)
  /// in -pi / pi
  ///
  /// @return a singlevector containing the vertices from one side and then
  /// from the other side consecutively
  std::vector<Vector3> circleVertices(const Transform3 transform,
                                      unsigned int quarterSegments) const;

  /// Output Method for std::ostream
  std::ostream& toStream(std::ostream& sl) const final;

 private:
  /// The bound radius, half Z, half phi and average phi
  std::array<double, eSize> m_values;
  /// Indicator if the bounds are closed
  bool m_closed{false};

  /// Check the input values for consistency, will throw a logic_exception
  /// if consistency is not given
  void checkConsistency() noexcept(false);

  /// Helper method to shift into the phi-frame
  /// @param lposition the polar coordinates in the global frame
  Vector2 shifted(const Vector2& lposition) const;
};

}  // namespace Acts
