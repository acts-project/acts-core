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

#include <algorithm>
#include <cmath>
#include <numbers>
#include <optional>
#include <span>
#include <utility>
#include <vector>

/// Helper methods for polyhedron vertices drawing and inside/outside checks.
namespace Acts::detail::VerticesHelper {

/// A method that inserts the cartesian extrema points and segments
/// a curved segment into sub segments
///
/// @param phiMin the minimum phi value
/// @param phiMax The second phi value
/// @param phiRef is a vector of reference phi values to be included as well
/// @param quarterSegments number of segments used to approximate a segment quarter
///
/// @return a vector of generated phi values
std::vector<double> phiSegments(double phiMin = -std::numbers::pi,
                                double phiMax = std::numbers::pi,
                                const std::vector<double>& phiRefs = {},
                                unsigned int quarterSegments = 2u);

/// Helper method to create a regular 2 or 3 D segment
/// between two phi values with a given number of segments
///
/// It will insert the phi at extrema points and reference points, it uses
/// a minimum approximation of a circle with 8 segments
///
/// @tparam vertex_t Type of vertex to be applied
/// @tparam transform_t Optional transform
///
/// @param rXY The radius description if first +/= second: ellipse
/// @param phiMin the minimum phi value
/// @param phiMax the second phi value
/// @param phiRef is a vector of reference phi values to be included as well
/// @param quarterSegments number of segments used to approximate a segment quarter
/// @param offset The out of plane offset position of the bow
/// @param transform The transform applied (optional)
///
/// @return a vector of vertices
template <typename vertex_t, typename transform_t>
std::vector<vertex_t> segmentVertices(
    std::pair<double, double> rXY, double phiMin, double phiMax,
    const std::vector<double>& phiRefs = {}, unsigned int quarterSegments = 2u,
    const vertex_t& offset = vertex_t::Zero(),
    const transform_t& transform = transform_t::Identity()) {
  std::vector<vertex_t> vertices;
  std::vector<double> phis =
      phiSegments(phiMin, phiMax, phiRefs, quarterSegments);
  for (double phi : phis) {
    vertex_t vertex = vertex_t::Zero();
    vertex(0) = rXY.first * std::cos(phi);
    vertex(1) = rXY.second * std::sin(phi);
    vertex = vertex + offset;
    vertices.push_back(transform * vertex);
  }
  return vertices;
}

/// Construct vertices on an ellipse-like bound object.
///
/// @param innerRx The radius of the inner ellipse (in x), 0 if sector
/// @param innerRy The radius of the inner ellipse (in y), 0 if sector
/// @param outerRx The radius of the outer ellipse (in x)
/// @param outerRy The radius of the outer ellipse (in y)
/// @param avgPhi The phi direction of the center if sector
/// @param halfPhi The half phi sector of the ellipse
/// @param quarterSegments number of segments used to approximate a segment quarter
///
/// @return a vector of 2d-vectors
std::vector<Vector2> ellipsoidVertices(double innerRx, double innerRy,
                                       double outerRx, double outerRy,
                                       double avgPhi = 0.,
                                       double halfPhi = std::numbers::pi,
                                       unsigned int quarterSegments = 2u);

/// Construct vertices on an disc/wheel-like bound object.
///
/// @param innerR The radius of the inner circle (sector)
/// @param outerR The radius of the outer circle (sector)
/// @param avgPhi The phi direction of the center if sector
/// @param halfPhi The half phi sector if sector
/// @param quarterSegments number of segments used to approximate a segment quarter
///
/// @return a vector of 2d-vectors
std::vector<Vector2> circularVertices(double innerR, double outerR,
                                      double avgPhi = 0.,
                                      double halfPhi = std::numbers::pi,
                                      unsigned int quarterSegments = 2u);

/// Check if the point is inside the polygon w/o any tolerances.
///
/// @tparam vertex_container_t is an iterable container
///
/// @param point is the Vector2Type to check
/// @param vertices Forward iterable container of convex polygon vertices.
///                 Calling `std::begin`/ `std::end` on the container must
///                 return an iterator where `*it` must be convertible to
///                 an `Vector2Type`.
/// @return bool for inside/outside
template <typename vertex_t, typename vertex_container_t>
bool isInsidePolygon(const vertex_t& point,
                     const vertex_container_t& vertices) {
  // when we move along the edges of a convex polygon, a point on the inside of
  // the polygon will always appear on the same side of each edge.
  // a point on the outside will switch sides at least once.

  // returns which side of the connecting line between `ll0` and `ll1` the point
  // `p` is on. computes the sign of the z-component of the cross-product
  // between the line normal vector and the vector from `ll0` to `p`.
  auto lineSide = [&](auto&& ll0, auto&& ll1) {
    auto normal = ll1 - ll0;
    auto delta = point - ll0;
    return std::signbit((normal[0] * delta[1]) - (normal[1] * delta[0]));
  };

  auto iv = std::begin(vertices);
  auto l0 = *iv;
  auto l1 = *(++iv);
  // use vertex0 to vertex1 to define reference sign and compare w/ all edges
  auto reference = lineSide(l0, l1);
  for (++iv; iv != std::end(vertices); ++iv) {
    l0 = l1;
    l1 = *iv;
    if (lineSide(l0, l1) != reference) {
      return false;
    }
  }
  // manual check for last edge from last vertex back to the first vertex
  if (lineSide(l1, *std::begin(vertices)) != reference) {
    return false;
  }
  // point was always on the same side. point must be inside.
  return true;
}

/// Check if the point is inside the rectangle.
///
/// @tparam vertex_t is vector with [0],[1] access
///
/// @param point is the Vector2Type to check
/// @param vertices Forward iterable container of convex polygon vertices.
///                 Calling `std::begin`/ `std::end` on the container must
///                 return an iterator where `*it` must be convertible to
///                 an `Vector2Type`.
/// @return bool for inside/outside
template <typename vertex_t>
bool isInsideRectangle(const vertex_t& point, const vertex_t& lowerLeft,
                       const vertex_t& upperRight) {
  return (lowerLeft[0] <= point[0]) && (point[0] < upperRight[0]) &&
         (lowerLeft[1] <= point[1]) && (point[1] < upperRight[1]);
}

/// This method checks if a cloud of points are on 2D hyper-plane in 3D space.
///
/// @param vertices The list of vertices to test
/// @param tolerance The allowed out of plane tolerance
/// @return boolean to indicate if all points are inside/outside
bool onHyperPlane(const std::vector<Vector3>& vertices,
                  double tolerance = s_onSurfaceTolerance);

/// Calculate the closest point on the polygon.
inline Vector2 computeClosestPointOnPolygon(const Vector2& point,
                                            std::span<const Vector2> vertices,
                                            const SquareMatrix2& metric) {
  auto squaredNorm = [&](const Vector2& x) {
    return (x.transpose() * metric * x).value();
  };

  // calculate the closest position on the segment between `ll0` and `ll1` to
  // the point as measured by the metric induced by the metric matrix
  auto closestOnSegment = [&](auto&& ll0, auto&& ll1) {
    // normal vector and position of the closest point along the normal
    auto n = ll1 - ll0;
    auto n_transformed = metric * n;
    auto f = n.dot(n_transformed);
    auto u = std::isnormal(f)
                 ? (point - ll0).dot(n_transformed) / f
                 : 0.5;  // ll0 and ll1 are so close it doesn't matter
    // u must be in [0, 1] to still be on the polygon segment
    return ll0 + std::clamp(u, 0.0, 1.0) * n;
  };

  auto iv = std::begin(vertices);
  Vector2 l0 = *iv;
  Vector2 l1 = *(++iv);
  Vector2 closest = closestOnSegment(l0, l1);
  auto closestDist = squaredNorm(closest - point);
  // Calculate the closest point on other connecting lines and compare distances
  for (++iv; iv != std::end(vertices); ++iv) {
    l0 = l1;
    l1 = *iv;
    Vector2 current = closestOnSegment(l0, l1);
    auto currentDist = squaredNorm(current - point);
    if (currentDist < closestDist) {
      closest = current;
      closestDist = currentDist;
    }
  }
  // final edge from last vertex back to the first vertex
  Vector2 last = closestOnSegment(l1, *std::begin(vertices));
  if (squaredNorm(last - point) < closestDist) {
    closest = last;
  }
  return closest;
}

/// Calculate the closest point on the box
inline Vector2 computeEuclideanClosestPointOnRectangle(
    const Vector2& point, const Vector2& lowerLeft, const Vector2& upperRight) {
  /*
   *
   *        |                 |
   *   IV   |       V         | I
   *        |                 |
   *  ------------------------------
   *        |                 |
   *        |                 |
   *   VIII |     INSIDE      | VI
   *        |                 |
   *        |                 |
   *  ------------------------------
   *        |                 |
   *   III  |      VII        | II
   *        |                 |
   *
   */

  double l0 = point[0], l1 = point[1];
  double loc0Min = lowerLeft[0], loc0Max = upperRight[0];
  double loc1Min = lowerLeft[1], loc1Max = upperRight[1];

  // check if inside
  if (loc0Min <= l0 && l0 < loc0Max && loc1Min <= l1 && l1 < loc1Max) {
    // INSIDE
    double dist = std::abs(loc0Max - l0);
    Vector2 cls(loc0Max, l1);

    double test = std::abs(loc0Min - l0);
    if (test <= dist) {
      dist = test;
      cls = {loc0Min, l1};
    }

    test = std::abs(loc1Max - l1);
    if (test <= dist) {
      dist = test;
      cls = {l0, loc1Max};
    }

    test = std::abs(loc1Min - l1);
    if (test <= dist) {
      return {l0, loc1Min};
    }
    return cls;
  } else {
    // OUTSIDE, check sectors
    if (l0 > loc0Max) {
      if (l1 > loc1Max) {  // I
        return {loc0Max, loc1Max};
      } else if (l1 <= loc1Min) {  // II
        return {loc0Max, loc1Min};
      } else {  // VI
        return {loc0Max, l1};
      }
    } else if (l0 < loc0Min) {
      if (l1 > loc1Max) {  // IV
        return {loc0Min, loc1Max};
      } else if (l1 <= loc1Min) {  // III
        return {loc0Min, loc1Min};
      } else {  // VIII
        return {loc0Min, l1};
      }
    } else {
      if (l1 > loc1Max) {  // V
        return {l0, loc1Max};
      } else {  // l1 <= loc1Min # VII
        return {l0, loc1Min};
      }
      // third case not necessary, see INSIDE above
    }
  }
}

inline Vector2 computeClosestPointOnAlignedBox(
    const Vector2& lowerLeft, const Vector2& upperRight, const Vector2& point,
    const std::optional<SquareMatrix2>& metricOpt) {
  Vector2 closestPoint;

  if (!metricOpt.has_value()) {
    closestPoint =
        detail::VerticesHelper::computeEuclideanClosestPointOnRectangle(
            point, lowerLeft, upperRight);
  } else {
    // TODO there might be a more optimal way to compute the closest point to a
    // box with metric

    std::array<Vector2, 4> vertices = {{lowerLeft,
                                        {upperRight[0], lowerLeft[1]},
                                        upperRight,
                                        {lowerLeft[0], upperRight[1]}}};

    closestPoint = detail::VerticesHelper::computeClosestPointOnPolygon(
        point, vertices, *metricOpt);
  }

  return closestPoint;
}

}  // namespace Acts::detail::VerticesHelper
