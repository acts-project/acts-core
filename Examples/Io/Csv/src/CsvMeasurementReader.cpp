// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvMeasurementReader.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Io/Csv/CsvInputOutput.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iterator>
#include <list>
#include <stdexcept>
#include <vector>

#include "CsvOutputData.hpp"

ActsExamples::CsvMeasurementReader::CsvMeasurementReader(
    const ActsExamples::CsvMeasurementReader::Config& config,
    Acts::Logging::Level level)
    : m_cfg(config),
      m_eventsRange(
          determineEventFilesRange(m_cfg.inputDir, "measurements.csv")),
      m_logger(Acts::getDefaultLogger("CsvMeasurementReader", level)) {
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurement output collection");
  }

  m_outputMeasurements.initialize(m_cfg.outputMeasurements);
  m_outputMeasurementSimHitsMap.initialize(m_cfg.outputMeasurementSimHitsMap);
  m_outputClusters.maybeInitialize(m_cfg.outputClusters);
  m_outputMeasurementParticlesMap.maybeInitialize(
      m_cfg.outputMeasurementParticlesMap);
  m_inputHits.maybeInitialize(m_cfg.inputSimHits);

  // Check if event ranges match (should also catch missing files)
  auto checkRange = [&](const std::string& fileStem) {
    const auto hitmapRange = determineEventFilesRange(m_cfg.inputDir, fileStem);
    if (hitmapRange.first > m_eventsRange.first ||
        hitmapRange.second < m_eventsRange.second) {
      throw std::runtime_error("event range mismatch for 'event**-" + fileStem +
                               "'");
    }
  };

  checkRange("measurement-simhit-map.csv");
  if (!m_cfg.outputClusters.empty()) {
    checkRange("cells.csv");
  }
}

std::string ActsExamples::CsvMeasurementReader::CsvMeasurementReader::name()
    const {
  return "CsvMeasurementReader";
}

std::pair<std::size_t, std::size_t>
ActsExamples::CsvMeasurementReader::availableEvents() const {
  return m_eventsRange;
}

namespace {
struct CompareHitId {
  // support transparent comparison between identifiers and full objects
  using is_transparent = void;
  template <typename T>
  constexpr bool operator()(const T& left, const T& right) const {
    return left.hit_id < right.hit_id;
  }
  template <typename T>
  constexpr bool operator()(std::uint64_t left_id, const T& right) const {
    return left_id < right.hit_id;
  }
  template <typename T>
  constexpr bool operator()(const T& left, std::uint64_t right_id) const {
    return left.hit_id < right_id;
  }
};

struct CompareGeometryId {
  bool operator()(const ActsExamples::MeasurementData& left,
                  const ActsExamples::MeasurementData& right) const {
    return left.geometry_id < right.geometry_id;
  }
};

template <typename Data>
inline std::vector<Data> readEverything(
    const std::string& inputDir, const std::string& filename,
    const std::vector<std::string>& optionalColumns, std::size_t event) {
  std::string path = ActsExamples::perEventFilepath(inputDir, filename, event);
  ActsExamples::NamedTupleCsvReader<Data> reader(path, optionalColumns);

  std::vector<Data> everything;
  Data one;
  while (reader.read(one)) {
    everything.push_back(one);
  }

  return everything;
}

std::vector<ActsExamples::MeasurementData> readMeasurementsByGeometryId(
    const std::string& inputDir, std::size_t event) {
  // geometry_id and t are optional columns
  auto measurements = readEverything<ActsExamples::MeasurementData>(
      inputDir, "measurements.csv", {"geometry_id", "t"}, event);
  // sort same way they will be sorted in the output container
  std::ranges::sort(measurements, CompareGeometryId{});
  return measurements;
}

ActsExamples::ClusterContainer makeClusters(
    const std::unordered_multimap<std::size_t, ActsExamples::CellData>&
        cellDataMap,
    std::size_t nMeasurements) {
  using namespace ActsExamples;
  ClusterContainer clusters;

  for (auto index = 0ul; index < nMeasurements; ++index) {
    auto [begin, end] = cellDataMap.equal_range(index);

    // Fill the channels with the iterators
    Cluster cluster;
    cluster.channels.reserve(std::distance(begin, end));

    for (auto it = begin; it != end; ++it) {
      const auto& cellData = it->second;
      ActsFatras::Segmentizer::Segment2D dummySegment = {Acts::Vector2::Zero(),
                                                         Acts::Vector2::Zero()};

      ActsFatras::Segmentizer::Bin2D bin{
          static_cast<unsigned int>(cellData.channel0),
          static_cast<unsigned int>(cellData.channel1)};

      cluster.channels.emplace_back(bin, dummySegment, cellData.value);
    }

    // update the iterator

    // Compute cluster size
    if (!cluster.channels.empty()) {
      auto compareX = [](const auto& a, const auto& b) {
        return a.bin[0] < b.bin[0];
      };
      auto compareY = [](const auto& a, const auto& b) {
        return a.bin[1] < b.bin[1];
      };

      auto [minX, maxX] = std::minmax_element(cluster.channels.begin(),
                                              cluster.channels.end(), compareX);
      auto [minY, maxY] = std::minmax_element(cluster.channels.begin(),
                                              cluster.channels.end(), compareY);
      cluster.sizeLoc0 = 1 + maxX->bin[0] - minX->bin[0];
      cluster.sizeLoc1 = 1 + maxY->bin[1] - minY->bin[1];
    }

    clusters.push_back(cluster);
  }
  return clusters;
}

}  // namespace

ActsExamples::ProcessCode ActsExamples::CsvMeasurementReader::read(
    const ActsExamples::AlgorithmContext& ctx) {
  // hit_id in the files is not required to be neither continuous nor
  // monotonic. internally, we want continuous indices within [0,#hits)
  // to simplify data handling. to be able to perform this mapping we first
  // read all data into memory before converting to the internal event data
  // types.
  //
  // Note: the cell data is optional
  auto measurementData =
      readMeasurementsByGeometryId(m_cfg.inputDir, ctx.eventNumber);

  // Prepare containers for the hit data using the framework event data types
  MeasurementContainer tmpMeasurements;
  GeometryIdMultimap<ConstVariableBoundMeasurementProxy> orderedMeasurements;
  IndexMultimap<Index> measurementSimHitsMap;

  tmpMeasurements.reserve(measurementData.size());
  orderedMeasurements.reserve(measurementData.size());
  // Safe long as we have single particle to sim hit association
  measurementSimHitsMap.reserve(measurementData.size());

  auto measurementSimHitLinkData =
      readEverything<ActsExamples::MeasurementSimHitLink>(
          m_cfg.inputDir, "measurement-simhit-map.csv", {}, ctx.eventNumber);
  for (auto mshLink : measurementSimHitLinkData) {
    measurementSimHitsMap.emplace_hint(measurementSimHitsMap.end(),
                                       mshLink.measurement_id, mshLink.hit_id);
  }

  for (const MeasurementData& m : measurementData) {
    Acts::GeometryIdentifier geoId{m.geometry_id};

    // Create the measurement
    DigitizedParameters dParameters;
    for (unsigned int ipar = 0;
         ipar < static_cast<unsigned int>(Acts::eBoundSize); ++ipar) {
      if (((m.local_key) & (1 << (ipar + 1))) != 0) {
        dParameters.indices.push_back(static_cast<Acts::BoundIndices>(ipar));
        switch (ipar) {
          case static_cast<unsigned int>(Acts::eBoundLoc0): {
            dParameters.values.push_back(m.local0);
            dParameters.variances.push_back(m.var_local0);
          }; break;
          case static_cast<unsigned int>(Acts::eBoundLoc1): {
            dParameters.values.push_back(m.local1);
            dParameters.variances.push_back(m.var_local1);
          }; break;
          case static_cast<unsigned int>(Acts::eBoundPhi): {
            dParameters.values.push_back(m.phi);
            dParameters.variances.push_back(m.var_phi);
          }; break;
          case static_cast<unsigned int>(Acts::eBoundTheta): {
            dParameters.values.push_back(m.theta);
            dParameters.variances.push_back(m.var_theta);
          }; break;
          case static_cast<unsigned int>(Acts::eBoundTime): {
            dParameters.values.push_back(m.time);
            dParameters.variances.push_back(m.var_time);
          }; break;
          default:
            break;
        }
      }
    }

    // The measurement container is unordered and the index under which
    // the measurement will be stored is known before adding it.
    auto measurement = createMeasurement(tmpMeasurements, geoId, dParameters);

    // Due to the previous sorting of the raw hit data by geometry id, new
    // measurements should always end up at the end of the container. previous
    // elements were not touched; cluster indices remain stable and can
    // be used to identify the m.
    auto inserted = orderedMeasurements.emplace_hint(orderedMeasurements.end(),
                                                     geoId, measurement);
    if (std::next(inserted) != orderedMeasurements.end()) {
      ACTS_FATAL("Something went horribly wrong with the hit sorting");
      return ProcessCode::ABORT;
    }
  }

  MeasurementContainer measurements;
  for (auto& [_, meas] : orderedMeasurements) {
    measurements.emplaceMeasurement(meas.size(), meas.geometryId(), meas);
  }

  // Generate measurement-particles-map
  if (m_inputHits.isInitialized() &&
      m_outputMeasurementParticlesMap.isInitialized()) {
    const auto hits = m_inputHits(ctx);

    IndexMultimap<ActsFatras::Barcode> outputMap;

    for (const auto& [measIdx, hitIdx] : measurementSimHitsMap) {
      const auto& hit = hits.nth(hitIdx);
      outputMap.emplace(measIdx, hit->particleId());
    }

    m_outputMeasurementParticlesMap(ctx, std::move(outputMap));
  }

  // Write the data to the EventStore
  m_outputMeasurements(ctx, std::move(measurements));
  m_outputMeasurementSimHitsMap(ctx, std::move(measurementSimHitsMap));

  /////////////////////////
  // Cluster information //
  /////////////////////////

  if (m_cfg.outputClusters.empty()) {
    return ActsExamples::ProcessCode::SUCCESS;
  }

  std::vector<ActsExamples::CellData> cellData;

  // This allows seamless import of files created with an older version where
  // the measurement_id-column is still named hit_id
  try {
    cellData = readEverything<ActsExamples::CellData>(
        m_cfg.inputDir, "cells.csv", {"timestamp"}, ctx.eventNumber);
  } catch (std::runtime_error& e) {
    // Rethrow exception if it is not about the measurement_id-column
    if (std::string(e.what()).find("Missing header column 'measurement_id'") ==
        std::string::npos) {
      throw;
    }

    const auto oldCellData = readEverything<ActsExamples::CellDataLegacy>(
        m_cfg.inputDir, "cells.csv", {"timestamp"}, ctx.eventNumber);

    auto fromLegacy = [](const CellDataLegacy& old) {
      return CellData{old.geometry_id, old.hit_id,    old.channel0,
                      old.channel1,    old.timestamp, old.value};
    };

    cellData.resize(oldCellData.size());
    std::transform(oldCellData.begin(), oldCellData.end(), cellData.begin(),
                   fromLegacy);
  }

  std::unordered_multimap<std::size_t, ActsExamples::CellData> cellDataMap;
  for (const auto& cd : cellData) {
    cellDataMap.emplace(cd.measurement_id, cd);
  }

  auto clusters = makeClusters(cellDataMap, orderedMeasurements.size());
  m_outputClusters(ctx, std::move(clusters));

  return ActsExamples::ProcessCode::SUCCESS;
}
