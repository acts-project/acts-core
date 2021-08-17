// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// System include(s)
#include <algorithm>
#include <cmath>
#include <utility>

// Acts include(s).
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"

// SYCL plugin include(s)
#include "Acts/Plugins/Sycl/Seeding/CreateSeedsForGroupSycl.hpp"
#include "Acts/Plugins/Sycl/Seeding/Seedfinder.hpp"

namespace Acts::Sycl {
template <typename external_spacepoint_t>
Seedfinder<external_spacepoint_t>::Seedfinder(
    Acts::SeedfinderConfig<external_spacepoint_t> config,
    const Acts::Sycl::DeviceExperimentCuts& cuts,
    Acts::Sycl::QueueWrapper wrappedQueue)
    : m_config(config.toInternalUnits()),
      m_deviceCuts(cuts),
      m_wrappedQueue(std::move(wrappedQueue)) {
  // init m_config
  m_config.highland = 13.6f * std::sqrt(m_config.radLengthPerSeed) *
                      (1 + 0.038f * std::log(m_config.radLengthPerSeed));
  float maxScatteringAngle = m_config.highland / m_config.minPt;
  m_config.maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;
  m_config.pTPerHelixRadius = 300.f * m_config.bFieldInZ;
  m_config.minHelixDiameter2 =
      std::pow(m_config.minPt * 2 / m_config.pTPerHelixRadius, 2);
  m_config.pT2perRadius =
      std::pow(m_config.highland / m_config.pTPerHelixRadius, 2);

  auto seedFilterConfig = m_config.seedFilter->getSeedFilterConfig();

  // init m_deviceConfig
  m_deviceConfig = Acts::Sycl::detail::DeviceSeedfinderConfig{
      m_config.deltaRMin,
      m_config.deltaRMax,
      m_config.cotThetaMax,
      m_config.collisionRegionMin,
      m_config.collisionRegionMax,
      m_config.maxScatteringAngle2,
      m_config.sigmaScattering,
      m_config.minHelixDiameter2,
      m_config.pT2perRadius,
      seedFilterConfig.deltaInvHelixDiameter,
      seedFilterConfig.impactWeightFactor,
      seedFilterConfig.deltaRMin,
      seedFilterConfig.compatSeedWeight,
      m_config.impactMax,
      seedFilterConfig.compatSeedLimit,
  };
}

template <typename external_spacepoint_t>
template <typename sp_range_t>
std::vector<Acts::Seed<external_spacepoint_t>>
Seedfinder<external_spacepoint_t>::createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const {
  std::vector<Seed<external_spacepoint_t>> outputVec;

  // As a first step, we create Arrays of Structures (AoS)
  // that are easily comprehensible by the GPU. This allows us
  // less memory access operations than with simple (float) arrays.

  std::vector<detail::DeviceSpacePoint> deviceBottomSPs;
  std::vector<detail::DeviceSpacePoint> deviceMiddleSPs;
  std::vector<detail::DeviceSpacePoint> deviceTopSPs;

  std::vector<const Acts::InternalSpacePoint<external_spacepoint_t>*>
      bottomSPvec;
  std::vector<const Acts::InternalSpacePoint<external_spacepoint_t>*>
      middleSPvec;
  std::vector<const Acts::InternalSpacePoint<external_spacepoint_t>*> topSPvec;

  for (const Acts::InternalSpacePoint<external_spacepoint_t>* SP : bottomSPs) {
    bottomSPvec.push_back(SP);
    deviceBottomSPs.push_back(
        detail::DeviceSpacePoint{SP->x(), SP->y(), SP->z(), SP->radius(),
                                 SP->varianceR(), SP->varianceZ()});
  }

  for (const Acts::InternalSpacePoint<external_spacepoint_t>* SP : middleSPs) {
    middleSPvec.push_back(SP);
    deviceMiddleSPs.push_back(
        detail::DeviceSpacePoint{SP->x(), SP->y(), SP->z(), SP->radius(),
                                 SP->varianceR(), SP->varianceZ()});
  }

  for (const Acts::InternalSpacePoint<external_spacepoint_t>* SP : topSPs) {
    topSPvec.push_back(SP);
    deviceTopSPs.push_back(
        detail::DeviceSpacePoint{SP->x(), SP->y(), SP->z(), SP->radius(),
                                 SP->varianceR(), SP->varianceZ()});
  }

  std::vector<std::vector<detail::SeedData>> seeds;

  // Call the SYCL seeding algorithm
  createSeedsForGroupSycl(m_wrappedQueue, m_deviceConfig, m_deviceCuts,
                          deviceBottomSPs, deviceMiddleSPs, deviceTopSPs,
                          seeds);

  // Iterate through seeds returned by the SYCL algorithm and perform the last
  // step of filtering for fixed middle SP.
  std::vector<std::pair<
      float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>
      seedsPerSPM;
  for (size_t mi = 0; mi < seeds.size(); ++mi) {
    seedsPerSPM.clear();
    for (size_t j = 0; j < seeds[mi].size(); ++j) {
      auto& bottomSP = *(bottomSPvec[seeds[mi][j].bottom]);
      auto& middleSP = *(middleSPvec[mi]);
      auto& topSP = *(topSPvec[seeds[mi][j].top]);
      float weight = seeds[mi][j].weight;

      seedsPerSPM.emplace_back(std::make_pair(
          weight, std::make_unique<const InternalSeed<external_spacepoint_t>>(
                      bottomSP, middleSP, topSP, 0)));
    }
    m_config.seedFilter->filterSeeds_1SpFixed(seedsPerSPM,
                                              std::back_inserter(outputVec));
  }
  return outputVec;
}
}  // namespace Acts::Sycl
