// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/LineBounds.hpp"

#include "Acts/Surfaces/detail/VerticesHelper.hpp"

#include <iomanip>
#include <iostream>

namespace Acts {

std::vector<double> LineBounds::values() const {
  return {m_values.begin(), m_values.end()};
}

void LineBounds::checkConsistency() noexcept(false) {
  if (get(eR) < 0.) {
    throw std::invalid_argument("LineBounds: zero radius.");
  }
  if (get(eHalfLengthZ) <= 0.) {
    throw std::invalid_argument("LineBounds: zero/negative length.");
  }
}

bool LineBounds::inside(const Vector2& lposition) const {
  double r = get(LineBounds::eR);
  double halfLengthZ = get(LineBounds::eHalfLengthZ);
  return detail::VerticesHelper::isInsideRectangle(
      lposition, Vector2(-r, -halfLengthZ), Vector2(r, halfLengthZ));
}

Vector2 LineBounds::closestPoint(
    const Vector2& lposition,
    const std::optional<SquareMatrix2>& metric) const {
  double r = get(LineBounds::eR);
  double halfLengthZ = get(LineBounds::eHalfLengthZ);
  return detail::VerticesHelper::computeClosestPointOnAlignedBox(
      Vector2(-r, -halfLengthZ), Vector2(r, halfLengthZ), lposition, metric);
}

std::ostream& LineBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::LineBounds: (radius, halflengthInZ) = ";
  sl << "(" << get(LineBounds::eR) << ", " << get(LineBounds::eHalfLengthZ)
     << ")";
  sl << std::setprecision(-1);
  return sl;
}

}  // namespace Acts
