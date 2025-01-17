// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/EllipseBounds.hpp"

#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace Acts {

using VectorHelpers::perp;
using VectorHelpers::phi;

std::vector<double> EllipseBounds::values() const {
  return {m_values.begin(), m_values.end()};
}

void EllipseBounds::checkConsistency() noexcept(false) {
  if (get(eInnerRx) >= get(eOuterRx) || get(eInnerRx) < 0. ||
      get(eOuterRx) <= 0.) {
    throw std::invalid_argument("EllipseBounds: invalid along x axis");
  }
  if (get(eInnerRy) >= get(eOuterRy) || get(eInnerRy) < 0. ||
      get(eOuterRy) <= 0.) {
    throw std::invalid_argument("EllipseBounds: invalid along y axis.");
  }
  if (get(eHalfPhiSector) < 0. || get(eHalfPhiSector) > std::numbers::pi) {
    throw std::invalid_argument("EllipseBounds: invalid phi sector setup.");
  }
  if (get(eAveragePhi) != detail::radian_sym(get(eAveragePhi))) {
    throw std::invalid_argument("EllipseBounds: invalid phi positioning.");
  }
}

bool EllipseBounds::inside(const Vector2& lposition) const {
  double phi =
      detail::radian_sym(VectorHelpers::phi(lposition) - get(eAveragePhi));
  return (-get(eHalfPhiSector) <= phi) && (phi < get(eHalfPhiSector)) &&
         (square(lposition[eBoundLoc0] / get(eInnerRx)) +
          square(lposition[eBoundLoc1] / get(eInnerRy))) >= 1 &&
         (square(lposition[eBoundLoc0] / get(eOuterRx)) +
          square(lposition[eBoundLoc1] / get(eOuterRy))) < 1;
}

Vector2 EllipseBounds::closestPoint(
    const Vector2& /*lposition*/,
    const std::optional<SquareMatrix2>& /*metric*/) const {
  throw std::logic_error("Not implemented");
}

std::vector<Vector2> EllipseBounds::vertices(
    unsigned int quarterSegments) const {
  return detail::VerticesHelper::ellipsoidVertices(
      get(eInnerRx), get(eInnerRy), get(eOuterRx), get(eOuterRy),
      get(eAveragePhi), get(eHalfPhiSector), quarterSegments);
}

const RectangleBounds& EllipseBounds::boundingBox() const {
  return m_boundingBox;
}

std::ostream& EllipseBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::EllipseBounds:  (innerRadius0, outerRadius0, innerRadius1, "
        "outerRadius1, hPhiSector, averagePhi) = ";
  sl << "(" << get(eInnerRx) << ", " << get(eInnerRy) << ", " << get(eOuterRx)
     << ", " << get(eOuterRy) << ", " << get(eAveragePhi) << ", "
     << get(eHalfPhiSector) << ", " << get(eAveragePhi) << ")";
  sl << std::setprecision(-1);
  return sl;
}

}  // namespace Acts
