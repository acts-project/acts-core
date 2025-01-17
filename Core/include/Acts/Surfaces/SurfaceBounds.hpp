// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"

#include <cmath>
#include <optional>
#include <ostream>

namespace Acts {

/// @class SurfaceBounds
///
/// Interface for surface bounds.
///
/// Surface bounds provide:
/// - inside() checks
/// - distance to boundary calculations
/// - the BoundsType and a set of parameters to simplify persistency
///
class SurfaceBounds {
 public:
  /// @enum BoundsType
  /// This is nested to the SurfaceBounds, as also VolumeBounds will have
  /// Bounds Type.
  enum BoundsType : int {
    eCone = 0,
    eCylinder = 1,
    eDiamond = 2,
    eDisc = 3,
    eEllipse = 4,
    eLine = 5,
    eRectangle = 6,
    eTrapezoid = 7,
    eTriangle = 8,
    eDiscTrapezoid = 9,
    eConvexPolygon = 10,
    eAnnulus = 11,
    eBoundless = 12,
    eOther = 13
  };

  virtual ~SurfaceBounds() = default;

  /// Return the bounds type - for persistency optimization
  /// @return the bounds type
  virtual BoundsType type() const = 0;

  /// Check if the bound coordinates are cartesian
  /// @return true if the bound coordinates are cartesian
  virtual bool isCartesian() const = 0;

  /// Computes the bound to cartesian jacobian at a given local position
  /// @param lposition is the local position at which the jacobian is computed
  /// @return the bound to cartesian jacobian
  virtual SquareMatrix2 boundToCartesianJacobian(
      const Vector2& lposition) const = 0;

  /// Computes the bound to cartesian metric at a given local position
  /// @param lposition is the local position at which the metric is computed
  /// @return the bound to cartesian metric
  virtual SquareMatrix2 boundToCartesianMetric(
      const Vector2& lposition) const = 0;

  /// Access method for bound values, this is a dynamically sized
  /// vector containing the parameters needed to describe these bounds
  /// @return of the stored values for this SurfaceBounds object
  virtual std::vector<double> values() const = 0;

  /// Inside check for the bounds object
  /// @param lposition is the local position
  /// @return true if the local position is inside the bounds
  virtual bool inside(const Vector2& lposition) const = 0;

  /// Calculates the closest point on the bounds to a given local position
  /// @param lposition is the local position
  /// @param metric is the metric to be used for the distance calculation
  /// @return the closest point on the bounds
  virtual Vector2 closestPoint(
      const Vector2& lposition,
      const std::optional<SquareMatrix2>& metric) const = 0;

  /// Calculates the distance to the bounds from a given local position
  /// @param lposition is the local position
  /// @return the distance to the bounds
  virtual double distance(const Vector2& lposition) const {
    SquareMatrix2 metric = boundToCartesianMetric(lposition);

    Vector2 closest = closestPoint(lposition, metric);
    Vector2 diff = closest - lposition;
    return std::sqrt((diff.transpose() * metric * diff)(0, 0));
  }

  /// Inside check for the bounds object given a boundary tolerance.
  /// @param lposition is the local position
  /// @param boundaryTolerance is the boundary tolerance object
  /// @return true if the local position is inside the bounds and tolerance
  virtual bool inside(const Vector2& lposition,
                      const BoundaryTolerance& boundaryTolerance) const {
    using enum BoundaryTolerance::Mode;

    if (boundaryTolerance.isInfinite()) {
      return true;
    }

    BoundaryTolerance::Mode toleranceMode = boundaryTolerance.mode();
    bool strictlyInside = inside(lposition);

    if (toleranceMode == None) {
      return strictlyInside;
    }

    if (toleranceMode == Extend && strictlyInside) {
      return true;
    }

    std::optional<SquareMatrix2> jacobian;
    std::optional<SquareMatrix2> metric;
    if (boundaryTolerance.hasChi2Bound()) {
      SquareMatrix2 j = boundToCartesianJacobian(lposition);
      jacobian = j;
      metric = j.transpose() * boundaryTolerance.asChi2Bound().weight * j;
    } else if (!isCartesian()) {
      jacobian = boundToCartesianJacobian(lposition);
      metric = boundToCartesianMetric(lposition);
    }

    Vector2 closest = closestPoint(lposition, metric);
    Vector2 distance = closest - lposition;

    if (toleranceMode == Shrink) {
      return boundaryTolerance.isTolerated(distance, jacobian) &&
             strictlyInside;
    }
    return boundaryTolerance.isTolerated(distance, jacobian);
  }

  /// Output Method for std::ostream, to be overloaded by child classes
  /// @param os is the outstream in which the string dump is done
  virtual std::ostream& toStream(std::ostream& os) const = 0;

  friend bool operator==(const SurfaceBounds& lhs, const SurfaceBounds& rhs) {
    if (&lhs == &rhs) {
      return true;
    }
    return (lhs.type() == rhs.type()) && (lhs.values() == rhs.values());
  }

  friend std::ostream& operator<<(std::ostream& os, const SurfaceBounds& sb) {
    return sb.toStream(os);
  }
};

}  // namespace Acts
