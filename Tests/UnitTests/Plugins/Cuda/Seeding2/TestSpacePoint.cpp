// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Local include(s).
#include "TestSpacePoint.hpp"

// System include(s).
#include <cmath>
#include <limits>

/// Difference allowed on floating point numbers to still be treated equal
static constexpr float allowedDiff = std::numeric_limits<float>::epsilon() * 4;

bool operator==(const TestSpacePoint& a, const TestSpacePoint& b) {
  return ((std::abs(a.m_x - b.m_x) < allowedDiff) &&
          (std::abs(a.m_y - b.m_y) < allowedDiff) &&
          (std::abs(a.m_z - b.m_z) < allowedDiff) &&
          (std::abs(a.m_r - b.m_r) < allowedDiff) &&
          (std::abs(a.m_x - b.m_x) < allowedDiff) &&
          (a.m_surface == b.m_surface) &&
          (std::abs(a.m_varianceR - b.m_varianceR) < allowedDiff) &&
          (std::abs(a.m_varianceZ - b.m_varianceZ) < allowedDiff));
}
