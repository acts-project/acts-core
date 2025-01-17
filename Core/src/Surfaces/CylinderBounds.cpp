// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/CylinderBounds.hpp"

#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <utility>

namespace Acts {

using VectorHelpers::perp;
using VectorHelpers::phi;

std::vector<double> CylinderBounds::values() const {
  return {m_values.begin(), m_values.end()};
}

Vector2 CylinderBounds::shifted(const Vector2& lposition) const {
  return {detail::radian_sym((lposition[0] / get(eR)) - get(eAveragePhi)),
          lposition[1]};
}

bool CylinderBounds::inside(const Vector2& lposition) const {
  double bevelMinZ = get(eBevelMinZ);
  double bevelMaxZ = get(eBevelMaxZ);

  double halfLengthZ = get(eHalfLengthZ);
  double halfPhi = get(eHalfPhiSector);

  if (bevelMinZ == 0. || bevelMaxZ == 0.) {
    return detail::VerticesHelper::isInsideRectangle(
        shifted(lposition), Vector2(-halfPhi, -halfLengthZ),
        Vector2(halfPhi, halfLengthZ));
  }

  double radius = get(eR);
  // Beleved sides will unwrap to a trapezoid
  ///////////////////////////////////
  //  ________
  // /| .  . |\ r/phi
  // \|______|/ r/phi
  // -Z   0  Z
  ///////////////////////////////////
  double localx =
      lposition[0] > radius ? 2 * radius - lposition[0] : lposition[0];
  Vector2 shiftedlposition = shifted(lposition);
  if ((std::abs(shiftedlposition[0]) <= halfPhi &&
       std::abs(shiftedlposition[1]) <= halfLengthZ)) {
    return true;
  }

  if ((lposition[1] >= -(localx * std::tan(bevelMinZ) + halfLengthZ)) &&
      (lposition[1] <= (localx * std::tan(bevelMaxZ) + halfLengthZ))) {
    return true;
  }

  Vector2 lowerLeft = {-radius, -halfLengthZ};
  Vector2 middleLeft = {0., -(halfLengthZ + radius * std::tan(bevelMinZ))};
  Vector2 upperLeft = {radius, -halfLengthZ};
  Vector2 upperRight = {radius, halfLengthZ};
  Vector2 middleRight = {0., (halfLengthZ + radius * std::tan(bevelMaxZ))};
  Vector2 lowerRight = {-radius, halfLengthZ};
  Vector2 vertices[] = {lowerLeft,  middleLeft,  upperLeft,
                        upperRight, middleRight, lowerRight};

  return detail::VerticesHelper::isInsidePolygon(lposition, vertices);
}

Vector2 CylinderBounds::closestPoint(
    const Vector2& lposition,
    const std::optional<SquareMatrix2>& metric) const {
  double bevelMinZ = get(eBevelMinZ);
  double bevelMaxZ = get(eBevelMaxZ);
  double halfLengthZ = get(eHalfLengthZ);
  double radius = get(eR);

  Vector2 lowerLeft = {-radius, -halfLengthZ};
  Vector2 middleLeft = {0., -(halfLengthZ + radius * std::tan(bevelMinZ))};
  Vector2 upperLeft = {radius, -halfLengthZ};
  Vector2 upperRight = {radius, halfLengthZ};
  Vector2 middleRight = {0., (halfLengthZ + radius * std::tan(bevelMaxZ))};
  Vector2 lowerRight = {-radius, halfLengthZ};
  Vector2 vertices[] = {lowerLeft,  middleLeft,  upperLeft,
                        upperRight, middleRight, lowerRight};

  return detail::VerticesHelper::computeClosestPointOnPolygon(
      lposition, vertices, metric.value_or(SquareMatrix2::Identity()));
}

std::ostream& CylinderBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::CylinderBounds: (radius, halfLengthZ, halfPhiSector, "
        "averagePhi, bevelMinZ, bevelMaxZ) = ";
  sl << "(" << get(eR) << ", " << get(eHalfLengthZ) << ", ";
  sl << get(eHalfPhiSector) << ", " << get(eAveragePhi) << ", ";
  sl << get(eBevelMinZ) << ", " << get(eBevelMaxZ) << ")";
  sl << std::setprecision(-1);
  return sl;
}

std::vector<Vector3> CylinderBounds::circleVertices(
    const Transform3 transform, unsigned int quarterSegments) const {
  std::vector<Vector3> vertices;

  double avgPhi = get(eAveragePhi);
  double halfPhi = get(eHalfPhiSector);

  std::vector<double> phiRef = {};
  if (bool fullCylinder = coversFullAzimuth(); fullCylinder) {
    phiRef = {avgPhi};
  }

  // Write the two bows/circles on either side
  std::vector<int> sides = {-1, 1};
  for (auto& side : sides) {
    // Helper method to create the segment
    auto svertices = detail::VerticesHelper::segmentVertices(
        {get(eR), get(eR)}, avgPhi - halfPhi, avgPhi + halfPhi, phiRef,
        quarterSegments, Vector3(0., 0., side * get(eHalfLengthZ)), transform);
    vertices.insert(vertices.end(), svertices.begin(), svertices.end());
  }

  double bevelMinZ = get(eBevelMinZ);
  double bevelMaxZ = get(eBevelMaxZ);

  // Modify the vertices position if bevel is defined
  if ((bevelMinZ != 0. || bevelMaxZ != 0.) && vertices.size() % 2 == 0) {
    auto halfWay = vertices.end() - vertices.size() / 2;
    double mult{1};
    auto invTransform = transform.inverse();
    auto func = [&mult, &transform, &invTransform](Vector3& v) {
      v = invTransform * v;
      v(2) += v(1) * mult;
      v = transform * v;
    };
    if (bevelMinZ != 0.) {
      mult = std::tan(-bevelMinZ);
      std::for_each(vertices.begin(), halfWay, func);
    }
    if (bevelMaxZ != 0.) {
      mult = std::tan(bevelMaxZ);
      std::for_each(halfWay, vertices.end(), func);
    }
  }
  return vertices;
}

void CylinderBounds::checkConsistency() noexcept(false) {
  if (get(eR) <= 0.) {
    throw std::invalid_argument(
        "CylinderBounds: invalid radial setup: radius is negative");
  }
  if (get(eHalfLengthZ) <= 0.) {
    throw std::invalid_argument(
        "CylinderBounds: invalid length setup: half length is negative");
  }
  if (get(eHalfPhiSector) <= 0. || get(eHalfPhiSector) > std::numbers::pi) {
    throw std::invalid_argument("CylinderBounds: invalid phi sector setup.");
  }
  if (get(eAveragePhi) != detail::radian_sym(get(eAveragePhi)) &&
      std::abs(std::abs(get(eAveragePhi)) - std::numbers::pi) > s_epsilon) {
    throw std::invalid_argument("CylinderBounds: invalid phi positioning.");
  }
  if (get(eBevelMinZ) != detail::radian_sym(get(eBevelMinZ))) {
    throw std::invalid_argument("CylinderBounds: invalid bevel at min Z.");
  }
  if (get(eBevelMaxZ) != detail::radian_sym(get(eBevelMaxZ))) {
    throw std::invalid_argument("CylinderBounds: invalid bevel at max Z.");
  }
}

}  // namespace Acts
