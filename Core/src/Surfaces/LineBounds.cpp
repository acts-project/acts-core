// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/LineBounds.hpp"

#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/detail/BoundaryCheckHelper.hpp"

#include <iomanip>
#include <iostream>

namespace Acts {

SurfaceBounds::BoundsType LineBounds::type() const {
  return SurfaceBounds::eLine;
}

bool LineBounds::inside(const Vector2& lposition) const {
  double r = get(LineBounds::eR);
  double halfLengthZ = get(LineBounds::eHalfLengthZ);
  return detail::insideAlignedBox(
      Vector2(-r, -halfLengthZ), Vector2(r, halfLengthZ),
      BoundaryTolerance::None(), lposition, std::nullopt);
}

Vector2 LineBounds::closestPoint(
    const Vector2& lposition,
    const std::optional<SquareMatrix2>& metric) const {
  double r = get(LineBounds::eR);
  double halfLengthZ = get(LineBounds::eHalfLengthZ);
  return detail::computeClosestPointOnAlignedBox(
      Vector2(-r, -halfLengthZ), Vector2(r, halfLengthZ), lposition, metric);
}

bool LineBounds::inside(const Vector2& lposition,
                        const BoundaryTolerance& boundaryTolerance) const {
  double r = get(LineBounds::eR);
  double halfLengthZ = get(LineBounds::eHalfLengthZ);
  return detail::insideAlignedBox(Vector2(-r, -halfLengthZ),
                                  Vector2(r, halfLengthZ), boundaryTolerance,
                                  lposition, std::nullopt);
}

std::ostream& LineBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "LineBounds: (radius, halflengthInZ) = ";
  sl << "(" << get(LineBounds::eR) << ", " << get(LineBounds::eHalfLengthZ)
     << ")";
  sl << std::setprecision(-1);
  return sl;
}

}  // namespace Acts
