// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <memory>

namespace Acts {

struct SpacePointGridConfig {
  // magnetic field
  float bFieldInZ;
  // minimum pT to be found by seedfinder
  float minPt;
  // maximum extension of sensitive detector layer relevant for seeding as
  // distance from x=y=0 (i.e. in r)
  float rMax;
  // maximum extension of sensitive detector layer relevant for seeding in
  // positive direction in z
  float zMax;
  // maximum extension of sensitive detector layer relevant for seeding in
  // negative direction in z
  float zMin;
  // maximum distance in r from middle space point to bottom or top spacepoint
  float deltaRMax;
  // maximum forward direction expressed as cot(theta)
  float cotThetaMax;

  SpacePointGridConfig toInternalUnits() const {
    using namespace Acts::UnitLiterals;
    SpacePointGridConfig config = *this;

    config.bFieldInZ /= 1000_T;
    config.minPt /= 1_MeV;
    config.rMax /= 1_mm;
    config.zMax /= 1_mm;
    config.zMin /= 1_mm;
    config.deltaRMax /= 1_mm;

    return config;
  }
};
template <typename external_spacepoint_t>
using SpacePointGrid =
    detail::Grid<std::vector<std::unique_ptr<
                     const InternalSpacePoint<external_spacepoint_t>>>,
                 detail::Axis<detail::AxisType::Equidistant,
                              detail::AxisBoundaryType::Closed>,
                 detail::Axis<detail::AxisType::Equidistant,
                              detail::AxisBoundaryType::Bound>>;

class SpacePointGridCreator {
 public:
  template <typename external_spacepoint_t>
  static std::unique_ptr<SpacePointGrid<external_spacepoint_t>> createGrid(
      const Acts::SpacePointGridConfig& _config);
};
}  // namespace Acts
#include "Acts/Seeding/SpacePointGrid.ipp"
