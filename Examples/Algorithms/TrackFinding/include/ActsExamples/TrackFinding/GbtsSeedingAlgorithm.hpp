// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// basing off of SeedingOrtho

#pragma once

#include "Acts/Seeding/SeedFinderGbts.hpp"
#include "Acts/Seeding/SeedFinderGbtsConfig.hpp"
#include "ActsExamples/EventData/Cluster.hpp"

// in core
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <string>
#include <vector>

namespace ActsExamples {

class GbtsSeedingAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    std::vector<std::string> inputSpacePoints;

    std::string outputSeeds;

    Acts::SeedFinderGbtsConfig<SimSpacePoint> seedFinderConfig;
    Acts::SeedFinderOptions seedFinderOptions;

    std::string layerMappingFile;

    std::vector<Acts::GeometryIdentifier> geometrySelection;

    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

    std::map<std::pair<int, int>, std::pair<int, int>> ActsGbtsMap;

    bool fill_module_csv = false;

    std::string inputClusters;
  };

  // constructor:
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  GbtsSeedingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  // code to make the algorithm run
  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext &ctx) const override;

  // access to config
  const Config &config() const { return m_cfg; }

  // own class functions
  // make the map
  std::map<std::pair<int, int>, std::pair<int, int>> makeActsGbtsMap() const;
  // make the vector of space points with FTF Info
  std::vector<Acts::GbtsSP<SimSpacePoint>> MakeGbtsSpacePoints(
      const AlgorithmContext &ctx,
      std::map<std::pair<int, int>, std::pair<int, int>> map) const;
  // layer numbering
  std::vector<Acts::TrigInDetSiLayer> LayerNumbering() const;

 private:
  Config m_cfg;

  std::unique_ptr<Acts::GbtsGeometry<SimSpacePoint>> m_gbtsGeo;

  std::vector<std::unique_ptr<ReadDataHandle<SimSpacePointContainer>>>
      m_inputSpacePoints{};

  WriteDataHandle<SimSeedContainer> m_outputSeeds{this, "OutputSeeds"};

  ReadDataHandle<ClusterContainer> m_inputClusters{this, "InputClusters"};
};

}  // namespace ActsExamples
