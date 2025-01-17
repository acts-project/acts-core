// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include "Acts/Surfaces/ConvexPolygonBounds.hpp"

#include <iomanip>
#include <iostream>

namespace Acts {

TrapezoidBounds::TrapezoidBounds(double halfXnegY, double halfXposY,
                                 double halfY, double rotAngle) noexcept(false)
    : m_values({halfXnegY, halfXposY, halfY, rotAngle}),
      m_boundingBox(std::max(halfXnegY, halfXposY), halfY) {
  rotateBoundingBox();
  checkConsistency();
}

TrapezoidBounds::TrapezoidBounds(
    const std::array<double, eSize>& values) noexcept(false)
    : m_values(values),
      m_boundingBox(
          std::max(values[eHalfLengthXnegY], values[eHalfLengthXposY]),
          values[eHalfLengthY]) {
  rotateBoundingBox();
  checkConsistency();
}

std::vector<double> TrapezoidBounds::values() const {
  return {m_values.begin(), m_values.end()};
}

bool TrapezoidBounds::inside(const Vector2& lposition) const {
  const double hlXnY = get(TrapezoidBounds::eHalfLengthXnegY);
  const double hlXpY = get(TrapezoidBounds::eHalfLengthXposY);
  const double hlY = get(TrapezoidBounds::eHalfLengthY);
  const double rotAngle = get(TrapezoidBounds::eRotationAngle);

  const Vector2 extPosition = Eigen::Rotation2Dd(rotAngle) * lposition;
  const double x = extPosition[0];
  const double y = extPosition[1];

  if (std::abs(y) - hlY > 0) {
    // outside y range
    return false;
  }

  if (std::abs(x) - std::max(hlXnY, hlXpY) > 0) {
    // outside x range
    return false;
  }

  if (std::abs(x) - std::min(hlXnY, hlXpY) <= 0) {
    // inside x range
    return true;
  }

  // at this stage, the point can only be in the triangles
  // run slow-ish polygon check
  Vector2 vertices[] = {
      {-hlXnY, -hlY}, {hlXnY, -hlY}, {hlXpY, hlY}, {-hlXpY, hlY}};
  return detail::VerticesHelper::isInsidePolygon(extPosition, vertices);
}

Vector2 TrapezoidBounds::closestPoint(
    const Vector2& lposition,
    const std::optional<SquareMatrix2>& metric) const {
  const double hlXnY = get(TrapezoidBounds::eHalfLengthXnegY);
  const double hlXpY = get(TrapezoidBounds::eHalfLengthXposY);
  const double hlY = get(TrapezoidBounds::eHalfLengthY);
  const double rotAngle = get(TrapezoidBounds::eRotationAngle);

  const Vector2 extPosition = Eigen::Rotation2Dd(rotAngle) * lposition;

  Vector2 vertices[] = {
      {-hlXnY, -hlY}, {hlXnY, -hlY}, {hlXpY, hlY}, {-hlXpY, hlY}};

  Vector2 extClosest = detail::VerticesHelper::computeClosestPointOnPolygon(
      extPosition, vertices, metric.value_or(SquareMatrix2::Identity()));

  return Eigen::Rotation2Dd(-rotAngle) * extClosest;
}

std::vector<Vector2> TrapezoidBounds::vertices(
    unsigned int /*ignoredSegments*/) const {
  const double hlXnY = get(TrapezoidBounds::eHalfLengthXnegY);
  const double hlXpY = get(TrapezoidBounds::eHalfLengthXposY);
  const double hlY = get(TrapezoidBounds::eHalfLengthY);
  const double rotAngle = get(TrapezoidBounds::eRotationAngle);

  std::vector<Vector2> vertices = {
      {-hlXnY, -hlY}, {hlXnY, -hlY}, {hlXpY, hlY}, {-hlXpY, hlY}};
  for (auto& v : vertices) {
    v = Eigen::Rotation2Dd(-rotAngle) * v;
  }
  return vertices;
}

const RectangleBounds& TrapezoidBounds::boundingBox() const {
  return m_boundingBox;
}

std::ostream& TrapezoidBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::TrapezoidBounds:  (halfXnegY, halfXposY, halfY, rotAngle) = "
     << "(" << get(eHalfLengthXnegY) << ", " << get(eHalfLengthXposY) << ", "
     << get(eHalfLengthY) << ", " << get(eRotationAngle) << ")";
  sl << std::setprecision(-1);
  return sl;
}

void TrapezoidBounds::rotateBoundingBox() noexcept(false) {
  const double rotAngle = get(eRotationAngle);

  if (rotAngle != 0.) {
    m_boundingBox = ConvexPolygonBounds<4>(vertices()).boundingBox();
  }
}

void TrapezoidBounds::checkConsistency() noexcept(false) {
  if (get(eHalfLengthXnegY) <= 0. || get(eHalfLengthXposY) <= 0.) {
    throw std::invalid_argument("TrapezoidBounds: invalid local x setup");
  }
  if (get(eHalfLengthY) <= 0.) {
    throw std::invalid_argument("TrapezoidBounds: invalid local y setup");
  }
}

}  // namespace Acts
