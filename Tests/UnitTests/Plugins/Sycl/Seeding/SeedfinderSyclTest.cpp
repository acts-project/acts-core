// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Sycl/Seeding/Seedfinder.hpp"
#include "Acts/Plugins/Sycl/Utilities/QueueWrapper.hpp"
#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SeedfinderConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>

#include <boost/type_erasure/any_cast.hpp>

#include "ATLASCuts.hpp"
#include "CommandLineArguments.h"
#include "SpacePoint.hpp"

using namespace Acts::UnitLiterals;

auto readFile(const std::string& filename) -> std::vector<const SpacePoint*> {
  std::string line;
  std::vector<const SpacePoint*> readSP;

  std::ifstream spFile(filename);
  if (spFile.is_open()) {
    int id = 0;
    while (std::getline(spFile, line)) {
      std::stringstream ss(line);
      std::string linetype;
      ss >> linetype;
      if (linetype == "lxyz") {
        float x;
        float y;
        float z;
        float r;
        float varianceR;
        float varianceZ;
        int layer;
        ss >> layer >> x >> y >> z >> varianceR >> varianceZ;
        r = std::sqrt(x * x + y * y);
        float f22 = varianceR;
        float wid = varianceZ;
        float cov = wid * wid * .08333;
        if (cov < f22)
          cov = f22;
        if (std::abs(z) > 450.) {
          varianceZ = 9. * cov;
          varianceR = .06;
        } else {
          varianceR = 9. * cov;
          varianceZ = .06;
        }
        readSP.emplace_back(
            new SpacePoint(x, y, z, r, layer, varianceR, varianceZ, id));
        ++id;
      }
    }
  }
  return readSP;
}

template <typename external_spacepoint_t>
auto setupSeedfinderConfiguration()
    -> Acts::SeedfinderConfig<external_spacepoint_t> {
  Acts::SeedfinderConfig<SpacePoint> config;
  // silicon detector max
  config.rMax = 160._mm;
  config.deltaRMin = 5._mm;
  config.deltaRMax = 160._mm;
  config.collisionRegionMin = -250._mm;
  config.collisionRegionMax = 250._mm;
  config.zMin = -2800._mm;
  config.zMax = 2800._mm;
  config.maxSeedsPerSpM = 5;
  // 2.7 eta
  config.cotThetaMax = 7.40627;
  config.sigmaScattering = 1.00000;
  config.minPt = 500._MeV;
  config.bFieldInZ = 1.99724_T;
  config.beamPos = {-.5_mm, -.5_mm};
  config.impactMax = 10._mm;

  // for sycl
  config.nTrplPerSpBLimit = 100;
  config.nAvgTrplPerSpBLimit = 6;
  return config;
}

template <typename external_spacepoint_t>
auto setupSpacePointGridConfig(
    const Acts::SeedfinderConfig<external_spacepoint_t>& config)
    -> Acts::SpacePointGridConfig {
  Acts::SpacePointGridConfig gridConf{};
  gridConf.bFieldInZ = config.bFieldInZ;
  gridConf.minPt = config.minPt;
  gridConf.rMax = config.rMax;
  gridConf.zMax = config.zMax;
  gridConf.zMin = config.zMin;
  gridConf.deltaRMax = config.deltaRMax;
  gridConf.cotThetaMax = config.cotThetaMax;
  return gridConf;
}

auto main(int argc, char** argv) -> int {
  auto start_prep = std::chrono::system_clock::now();

  CommandLineArguments cmdlTool;
  cmdlTool.parse(argc, argv);

  if (!cmdlTool.inpFileName.empty() && !cmdlTool.inpFileExists) {
    std::cerr << "Input file not found\n";
    return -1;
  }

  auto spVec = readFile(cmdlTool.inpFileName);

  auto bottomBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>());
  auto topBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>());
  auto config = setupSeedfinderConfiguration<SpacePoint>();

  Acts::ATLASCuts<SpacePoint> atlasCuts = Acts::ATLASCuts<SpacePoint>();
  Acts::Sycl::DeviceExperimentCuts deviceAtlasCuts;
  config.seedFilter = std::make_unique<Acts::SeedFilter<SpacePoint>>(
      Acts::SeedFilter<SpacePoint>(Acts::SeedFilterConfig(), &atlasCuts));

  const Acts::Logging::Level logLvl =
      cmdlTool.csvFormat ? Acts::Logging::WARNING : Acts::Logging::INFO;
  Acts::Sycl::Seedfinder<SpacePoint> syclSeedfinder(
      config, deviceAtlasCuts,
      Acts::Sycl::QueueWrapper(
          cmdlTool.deviceName,
          Acts::getDefaultLogger("Sycl::QueueWrapper", logLvl)));
  Acts::Seedfinder<SpacePoint> normalSeedfinder(config);
  auto globalTool =
      [=](const SpacePoint& sp, float /*unused*/, float /*unused*/,
          float_t /*unused*/) -> std::pair<Acts::Vector3, Acts::Vector2> {
    Acts::Vector3 position(sp.x(), sp.y(), sp.z());
    Acts::Vector2 covariance(sp.varianceR, sp.varianceZ);
    return std::make_pair(position, covariance);
  };

  std::unique_ptr<Acts::SpacePointGrid<SpacePoint>> grid =
      Acts::SpacePointGridCreator::createGrid<SpacePoint>(
          setupSpacePointGridConfig(config));

  auto spGroup = Acts::BinnedSPGroup<SpacePoint>(
      spVec.begin(), spVec.end(), globalTool, bottomBinFinder, topBinFinder,
      std::move(grid), config);

  auto end_prep = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsec_prep = end_prep - start_prep;
  double prepTime = elapsec_prep.count();

  if (!cmdlTool.csvFormat) {
    std::cout << "read " << spVec.size() << " SP from file "
              << cmdlTool.inpFileName << std::endl;
    std::cout << "Preparation time: " << std::to_string(prepTime) << std::endl;
  }

  // -------------------------------------- //
  // ----------- EXECUTE ON CPU ----------- //
  // -------------------------------------- //

  auto start_cpu = std::chrono::system_clock::now();
  uint group_count = 0;
  std::vector<std::vector<Acts::Seed<SpacePoint>>> seedVector_cpu;

  if (!cmdlTool.onlyGpu) {
    decltype(normalSeedfinder)::State state;
    for (auto groupIt = spGroup.begin(); !(groupIt == spGroup.end());
         ++groupIt) {
      normalSeedfinder.createSeedsForGroup(
          state, std::back_inserter(seedVector_cpu.emplace_back()),
          groupIt.bottom(), groupIt.middle(), groupIt.top());
      group_count++;
      if (!cmdlTool.allgroup && group_count >= cmdlTool.groups) {
        break;
      }
    }
  }

  auto end_cpu = std::chrono::system_clock::now();

  if (!cmdlTool.csvFormat) {
    std::cout << "Analyzed " << group_count << " groups for CPU" << std::endl;
  }

  // -------------------------------------- //
  // -------- EXECUTE ON GPU - SYCL ------- //
  // -------------------------------------- //

  auto start_sycl = std::chrono::system_clock::now();

  group_count = 0;
  std::vector<std::vector<Acts::Seed<SpacePoint>>> seedVector_sycl;

  for (auto groupIt = spGroup.begin(); !(groupIt == spGroup.end()); ++groupIt) {
    seedVector_sycl.push_back(syclSeedfinder.createSeedsForGroup(
        groupIt.bottom(), groupIt.middle(), groupIt.top()));
    group_count++;
    if (!cmdlTool.allgroup && group_count >= cmdlTool.groups) {
      break;
    }
  }
  auto end_sycl = std::chrono::system_clock::now();

  if (!cmdlTool.csvFormat) {
    std::cout << "Analyzed " << group_count << " groups for SYCL" << std::endl;
  }

  std::chrono::duration<double> elapsec_cpu = end_cpu - start_cpu;
  double cpuTime = elapsec_cpu.count();

  std::chrono::duration<double> elapsec_sycl = end_sycl - start_sycl;
  double syclTime = elapsec_sycl.count();

  auto textWidth = 20;
  auto numWidth = 11;

  int nSeed_cpu = 0;
  int nSeed_sycl = 0;
  int nMatch = 0;

  if (cmdlTool.matches && !cmdlTool.onlyGpu) {
    for (auto& outVec : seedVector_cpu) {
      nSeed_cpu += outVec.size();
    }

    for (auto& outVec : seedVector_sycl) {
      nSeed_sycl += outVec.size();
    }

    for (size_t i = 0; i < seedVector_cpu.size(); i++) {
      auto regionVec_cpu = seedVector_cpu[i];
      auto regionVec_sycl = seedVector_sycl[i];

      std::vector<std::vector<SpacePoint>> seeds_cpu;
      std::vector<std::vector<SpacePoint>> seeds_sycl;

      for (const auto& sd : regionVec_cpu) {
        std::vector<SpacePoint> seed_cpu;
        seed_cpu.push_back(*(sd.sp()[0]));
        seed_cpu.push_back(*(sd.sp()[1]));
        seed_cpu.push_back(*(sd.sp()[2]));
        seeds_cpu.push_back(seed_cpu);
      }
      for (const auto& sd : regionVec_sycl) {
        std::vector<SpacePoint> seed_sycl;
        seed_sycl.push_back(*(sd.sp()[0]));
        seed_sycl.push_back(*(sd.sp()[1]));
        seed_sycl.push_back(*(sd.sp()[2]));
        seeds_sycl.push_back(seed_sycl);
      }

      for (auto seed : seeds_cpu) {
        for (auto other : seeds_sycl) {
          if (seed[0] == other[0] && seed[1] == other[1] &&
              seed[2] == other[2]) {
            nMatch++;
            break;
          }
        }
      }
    }
  }

  if (!cmdlTool.csvFormat) {
    std::cout << std::endl;
    std::cout
        << "------------------------- Time Metric -------------------------"
        << std::endl;
    std::cout << std::setw(textWidth) << " Device:";
    std::cout << std::setw(numWidth) << "CPU";
    std::cout << std::setw(numWidth) << "SYCL";
    std::cout << std::setw(textWidth) << "Speedup/ Agreement" << std::endl;
    std::cout << std::setw(textWidth) << " Time (s):";
    std::cout << std::setw(numWidth) << std::to_string(cpuTime);
    std::cout << std::setw(numWidth) << std::to_string(syclTime);
    std::cout << std::setw(textWidth) << std::to_string(cpuTime / syclTime);
    std::cout << std::endl;

    if (cmdlTool.matches && !cmdlTool.onlyGpu) {
      std::cout << std::setw(textWidth) << " Seeds found:";
      std::cout << std::setw(numWidth) << std::to_string(nSeed_cpu);
      std::cout << std::setw(numWidth) << std::to_string(nSeed_sycl);
      std::cout << std::setw(textWidth)
                << std::to_string(float(nMatch) / float(nSeed_cpu) * 100);
      std::cout << std::endl;
    }

    std::cout
        << "---------------------------------------------------------------"
        << std::endl;
    std::cout << std::endl;
  } else {
    std::cout << cpuTime << ',' << syclTime << ',' << cpuTime / syclTime << ','
              << nSeed_cpu << ',' << nSeed_sycl << ',' << nMatch << '\n';
  }

  for (const auto* S : spVec) {
    delete[] S;
  }

  return 0;
}
