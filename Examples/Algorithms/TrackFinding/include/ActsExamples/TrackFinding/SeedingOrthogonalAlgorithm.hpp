// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/SeedfinderConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include <string>
#include <vector>

namespace ActsExamples {

/// Construct track seeds from space points.
class SeedingOrthogonalAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input space point collections.
    ///
    /// We allow multiple space point collections to allow different parts of
    /// the detector to use different algorithms for space point construction,
    /// e.g. single-hit space points for pixel-like detectors or double-hit
    /// space points for strip-like detectors.
    std::vector<std::string> inputSpacePoints;
    /// Output track seed collection.
    std::string outputSeeds;
    /// Output proto track collection.
    std::string outputProtoTracks;
    float rMax = 200.;
    float deltaRMin = 1.;
    float deltaRMax = 60.;
    float collisionRegionMin = -250;
    float collisionRegionMax = 250.;
    float zMin = -2000.;
    float zMax = 2000.;
    float maxSeedsPerSpM = 1;
    float cotThetaMax = 7.40627;  // 2.7 eta
    float sigmaScattering = 50;
    float radLengthPerSeed = 0.1;
    float minPt = 500.;
    float bFieldInZ = 0.00199724;
    float beamPosX = 0;
    float beamPosY = 0;
    float impactMax = 3.;
  };

  /// Construct the seeding algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  SeedingOrthogonalAlgorithm(Config cfg, Acts::Logging::Level lvl);

  std::optional<std::pair<double, double>> validTriple(
    const Acts::InternalSpacePoint<ActsExamples::SimSpacePoint> & b,
    const Acts::InternalSpacePoint<ActsExamples::SimSpacePoint> & m,
    const Acts::InternalSpacePoint<ActsExamples::SimSpacePoint> & t
  ) const;

  std::vector<std::pair<float, std::unique_ptr<const Acts::InternalSeed<ActsExamples::SimSpacePoint>>>>
  helper(
    const Acts::InternalSpacePoint<ActsExamples::SimSpacePoint> & middle,
    std::vector<const Acts::InternalSpacePoint<ActsExamples::SimSpacePoint> *> &,
    std::vector<const Acts::InternalSpacePoint<ActsExamples::SimSpacePoint> *> &
  ) const;

  template<std::size_t NDims, typename tree_t>
  ActsExamples::SimSeedContainer strategy1(
    const std::vector<const Acts::InternalSpacePoint<ActsExamples::SimSpacePoint> *> spacePoints,
    const tree_t & tree
  ) const;

  template<std::size_t NDims, typename tree_t>
  ActsExamples::SimSeedContainer strategy2(
    const std::vector<const Acts::InternalSpacePoint<ActsExamples::SimSpacePoint> *> spacePoints,
    const tree_t & tree
  ) const;

  template<std::size_t NDims, typename tree_t>
  ActsExamples::SimSeedContainer strategy3(
    const std::vector<const Acts::InternalSpacePoint<ActsExamples::SimSpacePoint> *> spacePoints,
    const tree_t & tree
  ) const;


  /// Run the seeding algorithm.
  ///
  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final override;

 private:
  Config m_cfg;
  Acts::SpacePointGridConfig m_gridCfg;
  Acts::SeedfinderConfig<SimSpacePoint> m_finderCfg;
};

}  // namespace ActsExamples
