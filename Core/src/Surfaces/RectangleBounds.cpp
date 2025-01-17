// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/RectangleBounds.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace Acts {

double RectangleBounds::get(BoundValues bValue) const {
  switch (bValue) {
    case eMinX:
      return m_min.x();
    case eMinY:
      return m_min.y();
    case eMaxX:
      return m_max.x();
    case eMaxY:
      return m_max.y();
    default:
      assert(false && "Invalid BoundValue enum value");
      return std::numeric_limits<double>::quiet_NaN();
  }
}

std::vector<double> RectangleBounds::values() const {
  return {m_min.x(), m_min.y(), m_max.x(), m_max.y()};
}

void RectangleBounds::checkConsistency() noexcept(false) {
  if (get(eMinX) > get(eMaxX)) {
    throw std::invalid_argument("RectangleBounds: invalid local x setup");
  }
  if (get(eMinY) > get(eMaxY)) {
    throw std::invalid_argument("RectangleBounds: invalid local y setup");
  }
}

bool RectangleBounds::inside(const Vector2& lposition) const {
  return detail::VerticesHelper::isInsideRectangle(lposition, m_min, m_max);
}

Vector2 RectangleBounds::closestPoint(
    const Vector2& lposition,
    const std::optional<SquareMatrix2>& metric) const {
  // If no metric is provided we can use a shortcut for the rectangle
  if (!metric.has_value()) {
    return detail::VerticesHelper::computeEuclideanClosestPointOnRectangle(
        lposition, m_min, m_max);
  }

  // Otherwise we need to compute the closest point on the polygon
  std::array<Vector2, 4> vertices = {
      {m_min, {m_max[0], m_min[1]}, m_max, {m_min[0], m_max[1]}}};
  return detail::VerticesHelper::computeClosestPointOnPolygon(
      lposition, vertices, metric.value_or(SquareMatrix2::Identity()));
}

std::vector<Vector2> RectangleBounds::vertices(unsigned int /*lseg*/) const {
  // counter-clockwise starting from bottom-left corner
  return {m_min, {m_max.x(), m_min.y()}, m_max, {m_min.x(), m_max.y()}};
}

const RectangleBounds& RectangleBounds::boundingBox() const {
  return (*this);
}

std::ostream& RectangleBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::RectangleBounds:  (hlX, hlY) = "
     << "(" << 0.5 * (get(eMaxX) - get(eMinX)) << ", "
     << 0.5 * (get(eMaxY) - get(eMinY)) << ")";
  sl << "\n(lower left, upper right):\n";
  sl << min().transpose() << "\n" << max().transpose();
  sl << std::setprecision(-1);
  return sl;
}

}  // namespace Acts
