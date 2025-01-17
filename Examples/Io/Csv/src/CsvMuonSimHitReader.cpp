// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvMuonSimHitReader.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/MuonSimHit.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Io/Csv/CsvInputOutput.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Hit.hpp"

#include <array>
#include <stdexcept>

#include "CsvOutputData.hpp"

ActsExamples::CsvMuonSimHitReader::CsvMuonSimHitReader(
    const ActsExamples::CsvMuonSimHitReader::Config& config,
    Acts::Logging::Level level)
    : m_cfg(config),
      // TODO check that all files (hits,cells,truth) exists
      m_eventsRange(
          determineEventFilesRange(m_cfg.inputDir, m_cfg.inputStem + ".csv")),
      m_logger(Acts::getDefaultLogger("CsvMuonSimHitReader", level)) {
  if (m_cfg.inputStem.empty()) {
    throw std::invalid_argument("Missing input filename stem");
  }
  if (m_cfg.outputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits output collection");
  }

  m_outputSimHits.initialize(m_cfg.outputSimHits);
}

std::string ActsExamples::CsvMuonSimHitReader::CsvMuonSimHitReader::name()
    const {
  return "CsvMuonSimHitReader";
}

std::pair<std::size_t, std::size_t>
ActsExamples::CsvMuonSimHitReader::availableEvents() const {
  return m_eventsRange;
}

ActsExamples::ProcessCode ActsExamples::CsvMuonSimHitReader::read(
    const ActsExamples::AlgorithmContext& ctx) {
  auto path = perEventFilepath(m_cfg.inputDir, m_cfg.inputStem + ".csv",
                               ctx.eventNumber);

  ActsExamples::NamedTupleCsvReader<MuonSimHitData> reader(path);

  MuonSimHitData data;

  SimHitContainer::sequence_type unordered;
  while (reader.read(data)) {
    Acts::Vector4 pos{data.LocalPositionExtrx * Acts::UnitConstants::mm,
                      data.LocalPositionExtry * Acts::UnitConstants::mm,
                      data.LocalPositionExtrz * Acts::UnitConstants::mm, 0};
    Acts::Vector4 mom{data.LocalDirectionx * Acts::UnitConstants::GeV,
                      data.LocalDirectiony * Acts::UnitConstants::GeV,
                      data.LocalDirectionz * Acts::UnitConstants::GeV,
                      std::sqrt(data.LocalDirectionx * data.LocalDirectionx +
                                data.LocalDirectiony * data.LocalDirectiony +
                                data.LocalDirectionz * data.LocalDirectionz) *
                          Acts::UnitConstants::GeV};
    MuonMdtIdentifierFields f;
    f.multilayer = 0;
    f.tube = 0;
    f.tubeLayer = 0;
    f.stationEta = data.StationEta;
    f.stationPhi = data.StationPhi;
    f.stationName = data.StationName;

    unordered.push_back(SimHit(Acts::GeometryIdentifier(compressId(f)),
                               data.pdgId, pos, mom, mom, -1));
  }
  SimHitContainer simHits;
  simHits.insert(unordered.begin(), unordered.end());

  // write the ordered data to the EventStore (according to geometry_id).
  m_outputSimHits(ctx, std::move(simHits));

  return ActsExamples::ProcessCode::SUCCESS;
}
