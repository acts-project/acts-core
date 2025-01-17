// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <array>
#include <iosfwd>
#include <vector>

namespace Acts {

/// @class TrapezoidBounds
///
/// Bounds for a trapezoidal, planar Surface.
///
/// @image html TrapezoidBounds.gif
///
/// @todo can be speed optimized by calculating kappa/delta and caching it
class TrapezoidBounds : public PlanarBounds {
 public:
  enum BoundValues {
    eHalfLengthXnegY = 0,
    eHalfLengthXposY = 1,
    eHalfLengthY = 2,
    eRotationAngle = 3,
    eSize = 4
  };

  /// Constructor for symmetric Trapezoid
  ///
  /// @param halfXnegY minimal half length X, definition at negative Y
  /// @param halfXposY maximal half length X, definition at positive Y
  /// @param halfY half length Y - defined at x=0
  /// @param rotAngle: rotation angle of the bounds w.r.t coordinate axes
  TrapezoidBounds(double halfXnegY, double halfXposY, double halfY,
                  double rotAngle = 0.) noexcept(false);

  /// Constructor for symmetric Trapezoid - from fixed size array
  ///
  /// @param values the values to be stream in
  TrapezoidBounds(const std::array<double, eSize>& values) noexcept(false);

  /// @copydoc SurfaceBounds::type
  BoundsType type() const final { return eTrapezoid; }

  /// @copydoc SurfaceBounds::values
  std::vector<double> values() const final;

  /// @copydoc SurfaceBounds::inside
  ///
  /// The orientation of the Trapezoid is according to the figure above,
  /// in words: the shorter of the two parallel sides of the trapezoid
  /// intersects
  /// with the negative @f$ y @f$ - axis of the local frame.
  ///
  /// <br>
  /// The cases are:<br>
  /// (0) @f$ y @f$ or @f$ x @f$ bounds are 0 || 0<br>
  /// (1) the local position is outside @f$ y @f$ bounds <br>
  /// (2) the local position is inside @f$ y @f$ bounds, but outside maximum @f$
  /// x
  /// @f$ bounds  <br>
  /// (3) the local position is inside @f$ y @f$ bounds AND inside minimum @f$ x
  /// @f$ bounds <br>
  /// (4) the local position is inside @f$ y @f$ bounds AND inside maximum @f$ x
  /// @f$ bounds, so that it depends on the @f$ eta @f$ coordinate
  /// (5) the local position fails test of (4) <br>
  ///
  /// The inside check is done using single equations of straight lines and one
  /// has
  /// to take care if a point
  /// lies on the positive @f$ x @f$ half area(I) or the negative one(II).
  /// Denoting
  /// @f$ |x_{min}| @f$ and
  /// @f$ | x_{max} | @f$ as \c minHalfX respectively \c maxHalfX, such as @f$ |
  /// y_{H} | @f$ as \c halfY,
  /// the equations for the straing lines in (I) and (II) can be written as:<br>
  ///  <br>
  /// - (I):  @f$ y = \kappa_{I} x + \delta_{I} @f$ <br>
  /// - (II): @f$ y = \kappa_{II} x + \delta_{II} @f$ ,<br>
  ///  <br>
  /// where @f$  \kappa_{I} = - \kappa_{II} = 2 \frac{y_{H}}{x_{max} - x_{min}}
  /// @f$
  /// <br>
  /// and   @f$  \delta_{I} = \delta_{II} = - \frac{1}{2}\kappa_{I}(x_{max} +
  /// x_{min}) @f$
  bool inside(const Vector2& lposition) const final;

  /// @copydoc SurfaceBounds::closestPoint
  Vector2 closestPoint(const Vector2& lposition,
                       const std::optional<SquareMatrix2>& metric) const final;

  using SurfaceBounds::inside;

  /// Return the vertices
  ///
  /// @param ignoredSegments is and ignored parameter used to describe
  /// the number of segments to approximate curved sectors.
  ///
  /// @note the number of segments is ignored in this representation
  ///
  /// @return vector for vertices in 2D
  std::vector<Vector2> vertices(unsigned int ignoredSegments = 0u) const final;

  // Bounding box representation
  const RectangleBounds& boundingBox() const final;

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostream to be dumped into
  std::ostream& toStream(std::ostream& sl) const final;

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  double get(BoundValues bValue) const { return m_values[bValue]; }

 private:
  std::array<double, eSize> m_values;
  RectangleBounds m_boundingBox;

  void rotateBoundingBox() noexcept(false);

  /// Check the input values for consistency, will throw a logic_exception
  /// if consistency is not given
  void checkConsistency() noexcept(false);
};

}  // namespace Acts
