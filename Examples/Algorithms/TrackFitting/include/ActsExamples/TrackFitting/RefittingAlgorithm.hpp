// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/TrackFitting/TrackFitterFunction.hpp"

#include <memory>
#include <string>

namespace ActsExamples {
class TrackFitterFunction;
struct AlgorithmContext;

class RefittingAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// The input track collection
    std::string inputTracks;
    /// Output fitted tracks collection.
    std::string outputTracks;
    /// Type erased fitter function.
    std::shared_ptr<TrackFitterFunction> fit;
    /// Pick a single track for debugging (-1 process all tracks)
    int pickTrack = -1;
    /// Inflate initial covariance.
    std::array<double, 6> initialVarInflation = {1., 1., 1., 1., 1., 1.};
  };

  /// Constructor of the fitting algorithm
  ///
  /// @param config is the config struct to configure the algorithm
  /// @param level is the logging level
  RefittingAlgorithm(Config config, Acts::Logging::Level level);

  /// Framework execute method of the fitting algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algporithm flow
  ActsExamples::ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};
};

}  // namespace ActsExamples
