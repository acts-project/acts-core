// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/TrackFinding/SourceLinkAccessorConcept.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/MagneticField/MagneticField.hpp"

#include <functional>
#include <vector>

namespace ActsExamples {

class TrackFindingAlgorithm final : public BareAlgorithm {
 public:
  /// Track finder function that takes input measurements, initial trackstate
  /// and track finder options and returns some track-finder-specific result.
  using TrackFinderOptions =
      Acts::CombinatorialKalmanFilterOptions<IndexSourceLinkAccessor,
                                             MeasurementCalibrator,
                                             Acts::MeasurementSelector>;
  using TrackFinderResult = std::vector<
      Acts::Result<Acts::CombinatorialKalmanFilterResult<IndexSourceLink>>>;

  /// Find function that takes the above parameters
  /// @note This is separated into a virtual interface to keep compilation units
  /// small
  class TrackFinderFunction {
   public:
    virtual ~TrackFinderFunction() = default;
    virtual TrackFinderResult operator()(const IndexSourceLinkContainer&,
                                         const TrackParametersContainer&,
                                         const TrackFinderOptions&) const = 0;
  };

  /// Create the track finder function implementation.
  ///
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyways.
  static std::shared_ptr<TrackFinderFunction> makeTrackFinderFunction(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField);

  struct Config {
    /// Input measurements collection.
    std::string inputMeasurements;
    /// Input source links collection.
    std::string inputSourceLinks;
    /// Input initial track parameter estimates for for each proto track.
    std::string inputInitialTrackParameters;
    /// Output find trajectories collection.
    std::string outputTrajectories;
    /// Type erased track finder function.
    std::shared_ptr<TrackFinderFunction> findTracks;
    /// CKF measurement selector config
    Acts::MeasurementSelector::Config measurementSelectorCfg;
    /// Compute shared hit information
    bool computeSharedHits = false;
  };

  /// Constructor of the track finding algorithm
  ///
  /// @param config is the config struct to configure the algorithm
  /// @param level is the logging level
  TrackFindingAlgorithm(Config config, Acts::Logging::Level level);

  /// Framework execute method of the track finding algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algorithm flow
  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  template <typename source_link_accessor_container_t,
            typename source_link_accessor_value_t>
  void computeSharedHits(
      const source_link_accessor_container_t& sourcelinks,
      std::vector<Acts::Result<Acts::CombinatorialKalmanFilterResult<
          source_link_accessor_value_t>>>&) const;

 private:
  Config m_cfg;
};

template <typename source_link_accessor_container_t,
          typename source_link_accessor_value_t>
void TrackFindingAlgorithm::computeSharedHits(
    const source_link_accessor_container_t& sourceLinks,
    std::vector<Acts::Result<
        Acts::CombinatorialKalmanFilterResult<source_link_accessor_value_t>>>&
        results) const {
  // Compute shared hits from all the reconstructed tracks
  // Compute nSharedhits and Update ckf results
  // hit index -> list of multi traj indexes [traj, meas]
  static_assert(Acts::SourceLinkConcept<source_link_accessor_value_t>,
                "Source link does not fulfill SourceLinkConcept");

  std::vector<std::size_t> firstTrackOnTheHit(
      sourceLinks.size(), std::numeric_limits<std::size_t>::max());
  std::vector<std::size_t> firstStateOnTheHit(
      sourceLinks.size(), std::numeric_limits<std::size_t>::max());

  for (unsigned int iresult(0); iresult < results.size(); iresult++) {
    if (not results.at(iresult).ok()) {
      continue;
    }

    auto& ckfResult = results.at(iresult).value();
    auto& measIndexes = ckfResult.lastMeasurementIndices;

    for (auto measIndex : measIndexes) {
      ckfResult.fittedStates.visitBackwards(measIndex, [&](const auto& state) {
        if (not state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag))
          return;

        std::size_t hitIndex = state.uncalibrated().index();

        // Check if hit not already used
        if (firstTrackOnTheHit.at(hitIndex) ==
            std::numeric_limits<std::size_t>::max()) {
          firstTrackOnTheHit.at(hitIndex) = iresult;
          firstStateOnTheHit.at(hitIndex) = state.index();
          return;
        }

        // if already used, control if first track state has been marked
        // as shared
        int indexFirstTrack = firstTrackOnTheHit.at(hitIndex);
        int indexFirstState = firstStateOnTheHit.at(hitIndex);
        if (not results.at(indexFirstTrack)
                    .value()
                    .fittedStates.getTrackState(indexFirstState)
                    .typeFlags()
                    .test(Acts::TrackStateFlag::SharedHitFlag))
          results.at(indexFirstTrack)
              .value()
              .fittedStates.getTrackState(indexFirstState)
              .typeFlags()
              .set(Acts::TrackStateFlag::SharedHitFlag);

        // Decorate this track
        results.at(iresult)
            .value()
            .fittedStates.getTrackState(state.index())
            .typeFlags()
            .set(Acts::TrackStateFlag::SharedHitFlag);
      });
    }
  }
}

}  // namespace ActsExamples
