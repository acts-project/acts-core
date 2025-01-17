// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/DiamondBounds.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"

#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>

namespace Acts {

std::vector<double> DiamondBounds::values() const {
  return {m_values.begin(), m_values.end()};
}

void DiamondBounds::checkConsistency() noexcept(false) {
  if (std::ranges::any_of(m_values, [](auto v) { return v <= 0.; })) {
    throw std::invalid_argument("DiamondBounds: negative half length.");
  }
  if (get(eHalfLengthXnegY) > get(eHalfLengthXzeroY) ||
      get(eHalfLengthXposY) > get(eHalfLengthXzeroY)) {
    throw std::invalid_argument("DiamondBounds: not a diamond shape.");
  }
}

bool DiamondBounds::inside(const Vector2& lposition) const {
  // Vertices starting at lower left (min rel. phi)
  // counter-clockwise
  double x1 = get(DiamondBounds::eHalfLengthXnegY);
  double y1 = get(DiamondBounds::eHalfLengthYneg);
  double x2 = get(DiamondBounds::eHalfLengthXzeroY);
  double y2 = 0.;
  double x3 = get(DiamondBounds::eHalfLengthXposY);
  double y3 = get(DiamondBounds::eHalfLengthYpos);
  Vector2 vertices[] = {{-x1, -y1}, {x1, -y1}, {x2, y2},
                        {x3, y3},   {-x3, y3}, {-x2, y2}};
  return detail::VerticesHelper::isInsidePolygon(lposition, vertices);
}

Vector2 DiamondBounds::closestPoint(
    const Vector2& lposition,
    const std::optional<SquareMatrix2>& metric) const {
  // Vertices starting at lower left (min rel. phi)
  // counter-clockwise
  double x1 = get(DiamondBounds::eHalfLengthXnegY);
  double y1 = get(DiamondBounds::eHalfLengthYneg);
  double x2 = get(DiamondBounds::eHalfLengthXzeroY);
  double y2 = 0.;
  double x3 = get(DiamondBounds::eHalfLengthXposY);
  double y3 = get(DiamondBounds::eHalfLengthYpos);
  Vector2 vertices[] = {{-x1, -y1}, {x1, -y1}, {x2, y2},
                        {x3, y3},   {-x3, y3}, {-x2, y2}};
  return detail::VerticesHelper::computeClosestPointOnPolygon(
      lposition, vertices, metric.value_or(SquareMatrix2::Identity()));
}

std::vector<Vector2> DiamondBounds::vertices(
    unsigned int /*ignoredSegments*/) const {
  // Vertices starting at lower left (min rel. phi)
  // counter-clockwise
  double x1 = get(DiamondBounds::eHalfLengthXnegY);
  double y1 = get(DiamondBounds::eHalfLengthYneg);
  double x2 = get(DiamondBounds::eHalfLengthXzeroY);
  double y2 = 0.;
  double x3 = get(DiamondBounds::eHalfLengthXposY);
  double y3 = get(DiamondBounds::eHalfLengthYpos);
  return {{-x1, -y1}, {x1, -y1}, {x2, y2}, {x3, y3}, {-x3, y3}, {-x2, y2}};
}

const RectangleBounds& DiamondBounds::boundingBox() const {
  return m_boundingBox;
}

std::ostream& DiamondBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::DiamondBounds: (halfXatYneg, halfXatYzero, halfXatYpos, "
        "halfYneg, halfYpos) = ";
  sl << "(" << get(DiamondBounds::eHalfLengthXnegY) << ", "
     << get(DiamondBounds::eHalfLengthXzeroY) << ", "
     << get(DiamondBounds::eHalfLengthXposY) << ", "
     << get(DiamondBounds::eHalfLengthYneg) << ", "
     << get(DiamondBounds::eHalfLengthYpos) << ")";
  sl << std::setprecision(-1);
  return sl;
}

}  // namespace Acts
