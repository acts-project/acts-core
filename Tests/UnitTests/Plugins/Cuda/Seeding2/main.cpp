// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Local include(s).
#include "CommandLineArguments.hpp"
#include "ReadSeedFile.hpp"
#include "TestDeviceCuts.hpp"
#include "TestHostCuts.hpp"
#include "TestSpacePoint.hpp"

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Seeding2/SeedFinder.hpp"
#include "Acts/Plugins/Cuda/Utilities/Info.hpp"
#include "Acts/Plugins/Cuda/Utilities/MemoryManager.hpp"

// Acts include(s).
#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/SeedFilterConfig.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SeedfinderConfig.hpp"

// System include(s).
#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>

using namespace Acts::UnitLiterals;

int main(int argc, char* argv[]) {
  // Interpret the command line arguments passed to the executable.
  CommandLineArguments cmdl;
  cmdl.interpret(argc, argv);

  // Read in the seeds from the input text file.
  auto spacepoints = readSeedFile(cmdl.spFile, cmdl.filterDuplicates);
  std::cout << "Read " << spacepoints.size()
            << " spacepoints from file: " << cmdl.spFile << std::endl;

  // Create a "view vector" on top of them. This is necessary to be able to pass
  // the objects to the Acts code. While the return type of readSeedFile(...) is
  // useful for simplified memory management...
  std::vector<const TestSpacePoint*> spView;
  spView.reserve(spacepoints.size());
  for (const auto& sp : spacepoints) {
    spView.push_back(sp.get());
  }

  // Create binned groups of these spacepoints.
  auto bottomBinFinder = std::make_shared<Acts::BinFinder<TestSpacePoint>>();
  auto topBinFinder = std::make_shared<Acts::BinFinder<TestSpacePoint>>();

  // Set up the seedfinder configuration.
  Acts::SeedfinderConfig<TestSpacePoint> sfConfig;
  // silicon detector max
  sfConfig.rMax = 160._mm;
  sfConfig.deltaRMin = 5._mm;
  sfConfig.deltaRMax = 160._mm;
  sfConfig.collisionRegionMin = -250._mm;
  sfConfig.collisionRegionMax = 250._mm;
  sfConfig.zMin = -2800._mm;
  sfConfig.zMax = 2800._mm;
  sfConfig.maxSeedsPerSpM = 5;
  // 2.7 eta
  sfConfig.cotThetaMax = 7.40627;
  sfConfig.sigmaScattering = 1.00000;
  sfConfig.minPt = 500._MeV;
  sfConfig.bFieldInZ = 1.99724_T;
  sfConfig.beamPos = {-.5_mm, -.5_mm};
  sfConfig.impactMax = 10._mm;

  // Use a size slightly smaller than what modern GPUs are capable of. This is
  // because for debugging we can't use all available threads in a block, and
  // because early testing shows that using this sort of block size results in
  // better performance than using the maximal one. (It probably works better
  // with the kind of branching that is present in the CUDA code.)
  sfConfig.maxBlockSize = 256;

  // Set up the spacepoint grid configuration.
  Acts::SpacePointGridConfig gridConfig;
  gridConfig.bFieldInZ = sfConfig.bFieldInZ;
  gridConfig.minPt = sfConfig.minPt;
  gridConfig.rMax = sfConfig.rMax;
  gridConfig.zMax = sfConfig.zMax;
  gridConfig.zMin = sfConfig.zMin;
  gridConfig.deltaRMax = sfConfig.deltaRMax;
  gridConfig.cotThetaMax = sfConfig.cotThetaMax;

  // Covariance tool, sets covariances per spacepoint as required.
  auto ct = [=](const TestSpacePoint& sp, float, float,
                float) -> std::pair<Acts::Vector3, Acts::Vector2> {
    Acts::Vector3 position(sp.x(), sp.y(), sp.z());
    Acts::Vector2 covariance(sp.m_varianceR, sp.m_varianceZ);
    return std::make_pair(position, covariance);
  };

  // Create a grid with bin sizes according to the configured geometry, and
  // split the spacepoints into groups according to that grid.
  auto grid =
      Acts::SpacePointGridCreator::createGrid<TestSpacePoint>(gridConfig);
  auto spGroup = Acts::BinnedSPGroup<TestSpacePoint>(
      spView.begin(), spView.end(), ct, bottomBinFinder, topBinFinder,
      std::move(grid), sfConfig);
  // Make a convenient iterator that will be used multiple times later on.
  auto spGroup_end = spGroup.end();

  // Allocate memory on the selected CUDA device.
  if (Acts::Cuda::Info::instance().devices().size() <=
      static_cast<std::size_t>(cmdl.cudaDevice)) {
    std::cerr << "Invalid CUDA device (" << cmdl.cudaDevice << ") requested"
              << std::endl;
    return 1;
  }
  static constexpr std::size_t MEGABYTES = 1024l * 1024l;
  std::size_t deviceMemoryAllocation = cmdl.cudaDeviceMemory * MEGABYTES;
  if (deviceMemoryAllocation == 0) {
    deviceMemoryAllocation =
        Acts::Cuda::Info::instance().devices()[cmdl.cudaDevice].totalMemory *
        0.8;
  }
  std::cout << "Allocating " << deviceMemoryAllocation / MEGABYTES
            << " MB memory on device:\n"
            << Acts::Cuda::Info::instance().devices()[cmdl.cudaDevice]
            << std::endl;
  Acts::Cuda::MemoryManager::instance().setMemorySize(deviceMemoryAllocation,
                                                      cmdl.cudaDevice);

  // Set up the seedfinder configuration objects.
  TestHostCuts hostCuts;
  Acts::SeedFilterConfig filterConfig;
  sfConfig.seedFilter = std::make_unique<Acts::SeedFilter<TestSpacePoint>>(
      filterConfig, &hostCuts);
  auto deviceCuts = testDeviceCuts();

  // Set up the seedfinder objects.
  Acts::Seedfinder<TestSpacePoint> seedfinder_host(sfConfig);
  Acts::Cuda::SeedFinder<TestSpacePoint> seedfinder_device(
      sfConfig, filterConfig, deviceCuts, cmdl.cudaDevice);

  //
  // Perform the seed finding on the host.
  //

  // Record the start time.
  auto start_host = std::chrono::system_clock::now();
  // Create the result object.
  std::vector<std::vector<Acts::Seed<TestSpacePoint>>> seeds_host;

  // Perform the seed finding.
  if (!cmdl.onlyGPU) {
    auto spGroup_itr = spGroup.begin();
    decltype(seedfinder_host)::State state;
    for (std::size_t i = 0;
         spGroup_itr != spGroup_end && i < cmdl.groupsToIterate;
         ++i, ++spGroup_itr) {
      auto& group = seeds_host.emplace_back();
      seedfinder_host.createSeedsForGroup(
          state, std::back_inserter(group), spGroup_itr.bottom(),
          spGroup_itr.middle(), spGroup_itr.top());
    }
  }

  // Record the finish time.
  auto end_host = std::chrono::system_clock::now();
  double time_host = std::chrono::duration_cast<std::chrono::milliseconds>(
                         end_host - start_host)
                         .count() *
                     0.001;
  if (!cmdl.onlyGPU) {
    std::cout << "Done with the seedfinding on the host" << std::endl;
  }

  //
  // Perform the seed finding on the accelerator.
  //

  // Record the start time.
  auto start_device = std::chrono::system_clock::now();
  // Create the result object.
  std::vector<std::vector<Acts::Seed<TestSpacePoint>>> seeds_device;

  // Perform the seed finding.
  auto spGroup_itr = spGroup.begin();
  for (std::size_t i = 0;
       spGroup_itr != spGroup_end && i < cmdl.groupsToIterate;
       ++i, ++spGroup_itr) {
    seeds_device.push_back(seedfinder_device.createSeedsForGroup(
        spGroup_itr.bottom(), spGroup_itr.middle(), spGroup_itr.top()));
  }

  // Record the finish time.
  auto end_device = std::chrono::system_clock::now();
  double time_device = std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_device - start_device)
                           .count() *
                       0.001;
  std::cout << "Done with the seedfinding on the device" << std::endl;

  //
  // Print some summary about the seed finding.
  //

  // Count the total number of reconstructed seeds.
  std::size_t nSeeds_host = 0, nSeeds_device = 0;
  for (const auto& seeds : seeds_host) {
    nSeeds_host += seeds.size();
  }
  for (const auto& seeds : seeds_device) {
    nSeeds_device += seeds.size();
  }

  // Count how many seeds, reconstructed on the host, can be matched with seeds
  // reconstructed on the accelerator.
  std::size_t nMatch = 0;
  double matchPercentage = 0.0;
  if (!cmdl.onlyGPU) {
    assert(seeds_host.size() == seeds_device.size());
    for (size_t i = 0; i < seeds_host.size(); i++) {
      // Access the seeds for this region.
      const auto& seeds_in_host_region = seeds_host[i];
      const auto& seeds_in_device_region = seeds_device[i];
      // Loop over all seeds found on the host.
      for (const auto& host_seed : seeds_in_host_region) {
        assert(host_seed.sp().size() == 3);
        // Try to find a matching seed that was found on the accelerator.
        for (const auto& device_seed : seeds_in_device_region) {
          assert(device_seed.sp().size() == 3);
          if ((*(host_seed.sp()[0]) == *(device_seed.sp()[0])) &&
              (*(host_seed.sp()[1]) == *(device_seed.sp()[1])) &&
              (*(host_seed.sp()[2]) == *(device_seed.sp()[2]))) {
            ++nMatch;
            break;
          }
        }
      }
    }
    matchPercentage = (100.0 * nMatch) / nSeeds_host;
  }

  // Print the summary results.
  std::cout << std::endl;
  std::cout << "-------------------------- Results ---------------------------"
            << std::endl;
  std::cout << "|          |     Host     |    Device    | Speedup/agreement |"
            << std::endl;
  std::cout << "--------------------------------------------------------------"
            << std::endl;
  std::cout << "| Time [s] |  " << std::setw(10)
            << (cmdl.onlyGPU ? "N/A   " : std::to_string(time_host)) << "  |  "
            << std::setw(10) << time_device << "  |     " << std::setw(10)
            << (cmdl.onlyGPU ? "N/A    "
                             : std::to_string(time_host / time_device))
            << "    |" << std::endl;
  std::cout << "|   Seeds  |  " << std::setw(10)
            << (cmdl.onlyGPU ? "N/A   " : std::to_string(nSeeds_host))
            << "  |  " << std::setw(10) << nSeeds_device << "  |     "
            << std::setw(10)
            << (cmdl.onlyGPU ? "N/A    " : std::to_string(matchPercentage))
            << "    |" << std::endl;
  std::cout << "--------------------------------------------------------------"
            << std::endl;
  std::cout << std::endl;

  // Return gracefully.
  return 0;
}
