// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/RadialBounds.hpp"

#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace Acts {

std::vector<double> RadialBounds::values() const {
  return {m_values.begin(), m_values.end()};
}

void RadialBounds::checkConsistency() noexcept(false) {
  if (get(eMinR) < 0. || get(eMaxR) <= 0. || get(eMinR) > get(eMaxR)) {
    throw std::invalid_argument("RadialBounds: invalid radial setup");
  }
  if (get(eHalfPhiSector) < 0. || get(eHalfPhiSector) > std::numbers::pi) {
    throw std::invalid_argument("RadialBounds: invalid phi sector setup.");
  }
  if (get(eAveragePhi) != detail::radian_sym(get(eAveragePhi))) {
    throw std::invalid_argument("RadialBounds: invalid phi positioning.");
  }
}

SquareMatrix2 RadialBounds::boundToCartesianJacobian(
    const Vector2& lposition) const {
  SquareMatrix2 j;
  j(0, 0) = std::cos(lposition[1]);
  j(0, 1) = -lposition[0] * std::sin(lposition[1]);
  j(1, 0) = std::sin(lposition[1]);
  j(1, 1) = lposition[0] * std::cos(lposition[1]);
  return j;
}

SquareMatrix2 RadialBounds::boundToCartesianMetric(
    const Vector2& lposition) const {
  SquareMatrix2 m;
  m(0, 0) = 1;
  m(0, 1) = 0;
  m(1, 0) = 0;
  m(1, 1) = lposition[0] * lposition[0];
  return m;
}

Vector2 RadialBounds::shifted(const Vector2& lposition) const {
  Vector2 tmp;
  tmp[0] = lposition[0];
  tmp[1] = detail::radian_sym(lposition[1] - get(eAveragePhi));
  return tmp;
}

bool RadialBounds::inside(const Vector2& lposition) const {
  return detail::VerticesHelper::isInsideRectangle(
      shifted(lposition), Vector2(get(eMinR), -get(eHalfPhiSector)),
      Vector2(get(eMaxR), get(eHalfPhiSector)));
}

Vector2 RadialBounds::closestPoint(
    const Vector2& lposition,
    const std::optional<SquareMatrix2>& metric) const {
  return detail::VerticesHelper::computeClosestPointOnAlignedBox(
      Vector2(get(eMinR), -get(eHalfPhiSector)),
      Vector2(get(eMaxR), get(eHalfPhiSector)), shifted(lposition), metric);
}

std::vector<Vector2> RadialBounds::vertices(unsigned int lseg) const {
  return detail::VerticesHelper::circularVertices(
      get(eMinR), get(eMaxR), get(eAveragePhi), get(eHalfPhiSector), lseg);
}

std::ostream& RadialBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::RadialBounds:  (innerRadius, outerRadius, hPhiSector, "
        "averagePhi) = ";
  sl << "(" << get(eMinR) << ", " << get(eMaxR) << ", " << get(eHalfPhiSector)
     << ", " << get(eAveragePhi) << ")";
  sl << std::setprecision(-1);
  return sl;
}

}  // namespace Acts
