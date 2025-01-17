// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <vector>

namespace Acts {

/// Forward declare rectangle bounds as boundary box
class RectangleBounds;

/// @class PlanarBounds
///
/// common base class for all bounds that are in a local x/y cartesian frame
///  - simply introduced to avoid wrong bound assignments to surfaces
///
class PlanarBounds : public SurfaceBounds {
 public:
  /// Return the vertices
  ///
  /// @param quarterSegments is the number of segments used to describe curved
  /// segments in a quarter of the phi range. If it is 1, then only the extrema
  /// points in phi are inserted next to the segment corners.
  ///
  /// @note for planar bounds without curved segments @c quarterSegments is ignored
  ///
  /// @return vector for vertices in 2D
  virtual std::vector<Vector2> vertices(
      unsigned int quarterSegments = 2u) const = 0;

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

  /// Bounding box parameters
  /// @return rectangle bounds for a bounding box
  virtual const RectangleBounds& boundingBox() const = 0;
};

}  // namespace Acts
