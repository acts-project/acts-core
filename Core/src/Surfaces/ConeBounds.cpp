// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/ConeBounds.hpp"

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>

namespace Acts {

ConeBounds::ConeBounds(double alpha, bool symm, double halfphi,
                       double avphi) noexcept(false)
    : m_values({alpha, symm ? -std::numeric_limits<double>::infinity() : 0,
                std::numeric_limits<double>::infinity(), halfphi, avphi}),
      m_tanAlpha(std::tan(alpha)) {
  checkConsistency();
}

ConeBounds::ConeBounds(double alpha, double minz, double maxz, double halfphi,
                       double avphi) noexcept(false)
    : m_values({alpha, minz, maxz, halfphi, avphi}),
      m_tanAlpha(std::tan(alpha)) {
  checkConsistency();
}

ConeBounds::ConeBounds(const std::array<double, eSize>& values) noexcept(false)
    : m_values(values), m_tanAlpha(std::tan(values[eAlpha])) {
  checkConsistency();
}

std::vector<double> ConeBounds::values() const {
  return {m_values.begin(), m_values.end()};
}

void ConeBounds::checkConsistency() noexcept(false) {
  if (get(eAlpha) < 0. || get(eAlpha) >= std::numbers::pi) {
    throw std::invalid_argument("ConeBounds: invalid open angle.");
  }
  if (get(eMinZ) > get(eMaxZ) ||
      std::abs(get(eMinZ) - get(eMaxZ)) < s_epsilon) {
    throw std::invalid_argument("ConeBounds: invalid z range setup.");
  }
  if (get(eHalfPhiSector) < 0. || std::abs(eHalfPhiSector) > std::numbers::pi) {
    throw std::invalid_argument("ConeBounds: invalid phi sector setup.");
  }
  if (get(eAveragePhi) != detail::radian_sym(get(eAveragePhi))) {
    throw std::invalid_argument("ConeBounds: invalid phi positioning.");
  }
}

Vector2 ConeBounds::shifted(const Vector2& lposition) const {
  using detail::radian_sym;

  auto x = r(lposition[1]);  // cone radius at the local position
  Vector2 shifted;
  shifted[1] = lposition[1];
  shifted[0] = std::isnormal(x)
                   ? (x * radian_sym((lposition[0] / x) - get(eAveragePhi)))
                   : lposition[0];
  return shifted;
}

bool ConeBounds::inside(const Vector2& lposition) const {
  auto rphiHalf = r(lposition[1]) * get(eHalfPhiSector);
  return detail::VerticesHelper::isInsideRectangle(
      shifted(lposition), Vector2(-rphiHalf, get(eMinZ)),
      Vector2(rphiHalf, get(eMaxZ)));
}

Vector2 ConeBounds::closestPoint(
    const Vector2& lposition,
    const std::optional<SquareMatrix2>& metric) const {
  auto rphiHalf = r(lposition[1]) * get(eHalfPhiSector);
  return detail::VerticesHelper::computeClosestPointOnAlignedBox(
      Vector2(-rphiHalf, get(eMinZ)), Vector2(rphiHalf, get(eMaxZ)),
      shifted(lposition), metric);
}

std::ostream& ConeBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::ConeBounds: (tanAlpha, minZ, maxZ, halfPhiSector, averagePhi) "
        "= ";
  sl << "(" << m_tanAlpha << ", " << get(eMinZ) << ", " << get(eMaxZ) << ", "
     << get(eHalfPhiSector) << ", " << get(eAveragePhi) << ")";
  sl << std::setprecision(-1);
  return sl;
}

}  // namespace Acts
