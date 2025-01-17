// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace Acts {

DiscTrapezoidBounds::DiscTrapezoidBounds(double halfXminR, double halfXmaxR,
                                         double minR, double maxR,
                                         double avgPhi,
                                         double stereo) noexcept(false)
    : m_values({halfXminR, halfXmaxR, minR, maxR, avgPhi, stereo}) {
  checkConsistency();
  m_ymax = std::sqrt(get(eMaxR) * get(eMaxR) -
                     get(eHalfLengthXmaxR) * get(eHalfLengthXmaxR));
}

std::vector<double> DiscTrapezoidBounds::values() const {
  return {m_values.begin(), m_values.end()};
}

void DiscTrapezoidBounds::checkConsistency() noexcept(false) {
  if (get(eMinR) < 0. || get(eMaxR) <= 0. || get(eMinR) > get(eMaxR)) {
    throw std::invalid_argument("DiscTrapezoidBounds: invalid radial setup.");
  }
  if (get(eHalfLengthXminR) < 0. || get(eHalfLengthXmaxR) <= 0.) {
    throw std::invalid_argument("DiscTrapezoidBounds: negative length given.");
  }
  if (get(eAveragePhi) != detail::radian_sym(get(eAveragePhi))) {
    throw std::invalid_argument(
        "DiscTrapezoidBounds: invalid phi positioning.");
  }
}

Vector2 DiscTrapezoidBounds::toLocalCartesian(const Vector2& lposition) const {
  return {lposition[0] * std::sin(lposition[1] - get(eAveragePhi)),
          lposition[0] * std::cos(lposition[1] - get(eAveragePhi))};
}

SquareMatrix2 DiscTrapezoidBounds::boundToCartesianJacobian(
    const Vector2& lposition) const {
  SquareMatrix2 j;
  j(0, 0) = std::cos(lposition[1]);
  j(0, 1) = -lposition[0] * std::sin(lposition[1]);
  j(1, 0) = std::sin(lposition[1]);
  j(1, 1) = lposition[0] * std::cos(lposition[1]);
  return j;
}

SquareMatrix2 DiscTrapezoidBounds::boundToCartesianMetric(
    const Vector2& lposition) const {
  SquareMatrix2 m;
  m(0, 0) = 1;
  m(0, 1) = 0;
  m(1, 0) = 0;
  m(1, 1) = lposition[0] * lposition[0];
  return m;
}

bool DiscTrapezoidBounds::inside(const Vector2& lposition) const {
  Vector2 vertices[] = {{get(eHalfLengthXminR), get(eMinR)},
                        {get(eHalfLengthXmaxR), m_ymax},
                        {-get(eHalfLengthXmaxR), m_ymax},
                        {-get(eHalfLengthXminR), get(eMinR)}};
  return detail::VerticesHelper::isInsidePolygon(toLocalCartesian(lposition),
                                                 vertices);
}

Vector2 DiscTrapezoidBounds::closestPoint(
    const Vector2& lposition,
    const std::optional<SquareMatrix2>& metric) const {
  Vector2 vertices[] = {{get(eHalfLengthXminR), get(eMinR)},
                        {get(eHalfLengthXmaxR), m_ymax},
                        {-get(eHalfLengthXmaxR), m_ymax},
                        {-get(eHalfLengthXminR), get(eMinR)}};
  return detail::VerticesHelper::computeClosestPointOnPolygon(
      toLocalCartesian(lposition), vertices,
      metric.value_or(SquareMatrix2::Identity()));
}

std::vector<Vector2> DiscTrapezoidBounds::vertices(
    unsigned int /*ignoredSegments*/) const {
  Vector2 cAxis(std::cos(get(eAveragePhi)), std::sin(get(eAveragePhi)));
  Vector2 nAxis(cAxis.y(), -cAxis.x());
  auto ymin = std::sqrt(get(eMinR) * get(eMinR) -
                        get(eHalfLengthXminR) * get(eHalfLengthXminR));
  auto halfY = (m_ymax - ymin) / 2;
  return {-halfY * cAxis - get(eHalfLengthXminR) * nAxis,
          -halfY * cAxis + get(eHalfLengthXminR) * nAxis,
          halfY * cAxis + get(eHalfLengthXmaxR) * nAxis,
          halfY * cAxis - get(eHalfLengthXmaxR) * nAxis};
}

std::ostream& DiscTrapezoidBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::DiscTrapezoidBounds: (innerRadius, outerRadius, "
        "halfLengthXminR, "
        "halfLengthXmaxR, halfLengthY, halfPhiSector, averagePhi, rCenter, "
        "stereo) = ";
  sl << "(" << get(eMinR) << ", " << get(eMaxR) << ", " << get(eHalfLengthXminR)
     << ", " << get(eHalfLengthXmaxR) << ", " << halfLengthY() << ", "
     << halfPhiSector() << ", " << get(eAveragePhi) << ", " << rCenter() << ", "
     << stereo() << ")";
  sl << std::setprecision(-1);
  return sl;
}

}  // namespace Acts
