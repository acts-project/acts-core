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
#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <array>
#include <iosfwd>
#include <vector>

namespace Acts {

/// @class RectangleBounds
///
/// @image html RectangleBounds.gif
///
/// Bounds for a rectangular, planar surface - it can be used to for
/// rectangles that are symmetrically centered around (0./0.) and for
/// generic shifted rectangles
class RectangleBounds : public PlanarBounds {
 public:
  enum BoundValues : int {
    eMinX = 0,
    eMinY = 1,
    eMaxX = 2,
    eMaxY = 3,
    eSize = 4
  };

  /// Constructor with halflength in x and y - symmetric
  ///
  /// @param halfX halflength in X
  /// @param halfY halflength in Y
  RectangleBounds(double halfX, double halfY) noexcept(false)
      : m_min({-halfX, -halfY}), m_max({halfX, halfY}) {
    checkConsistency();
  }

  /// Constructor - from fixed size array - generic
  ///
  /// @param values The parameter values
  RectangleBounds(const std::array<double, eSize>& values) noexcept(false)
      : m_min({values[eMinX], values[eMinY]}),
        m_max({values[eMaxX], values[eMaxY]}) {
    checkConsistency();
  }

  /// Constructor - from min/max - generic
  ///
  /// @param min The left bottom corner
  /// @param max The right top corning
  RectangleBounds(const Vector2& min, const Vector2& max) noexcept(false)
      : m_min(min), m_max(max) {
    checkConsistency();
  }

  /// @copydoc SurfaceBounds::type
  BoundsType type() const final { return eRectangle; }

  /// @copydoc SurfaceBounds::values
  std::vector<double> values() const final;

  /// @copydoc SurfaceBounds::inside
  bool inside(const Vector2& lposition) const final;

  /// @copydoc SurfaceBounds::closestPoint
  Vector2 closestPoint(const Vector2& lposition,
                       const std::optional<SquareMatrix2>& metric) const final;

  using SurfaceBounds::inside;

  /// Return the vertices
  ///
  /// @param quarterSegments is the number of segments used to describe curved
  /// segments in a quarter of the phi range.
  /// @note the number of segments is ignored in this representation
  ///
  /// @return vector for vertices in 2D
  std::vector<Vector2> vertices(unsigned int quarterSegments = 0u) const final;

  // Bounding box representation
  const RectangleBounds& boundingBox() const final;

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostream for the dump
  std::ostream& toStream(std::ostream& sl) const final;

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  double get(BoundValues bValue) const;

  /// Access to the half length in X
  double halfLengthX() const { return 0.5 * (m_max.x() - m_min.x()); }

  /// Access to the half length in Y
  double halfLengthY() const { return 0.5 * (m_max.y() - m_min.y()); }

  /// Get the min vertex defining the bounds
  /// @return The min vertex
  const Vector2& min() const { return m_min; }

  /// Get the max vertex defining the bounds
  /// @return The max vertex
  const Vector2& max() const { return m_max; }

 private:
  Vector2 m_min;
  Vector2 m_max;

  /// Check the input values for consistency, will throw a logic_exception
  /// if consistency is not given
  void checkConsistency() noexcept(false);
};

}  // namespace Acts
