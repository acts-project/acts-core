// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <fstream>

namespace ActsExamples {

class CsvPropagationStepsWriter final
    : public ActsExamples::WriterT<std::vector<decltype(
          Acts::detail::SteppingLogger::this_result::steps)>> {
 public:
  using StepsVector =
      std::vector<decltype(Acts::detail::SteppingLogger::this_result::steps)>;

  struct Config {
    std::string collection;   ///< which collection to write
    std::string outputDir;    ///< where to place output files
    int outputPrecision = 6;  ///< floating point precision
  };

  /// Constructor with arguments
  ///
  /// @param cfg configuration struct
  /// @param level Output logging level
  CsvPropagationStepsWriter(const Config& cfg,
                            Acts::Logging::Level level = Acts::Logging::INFO)
      : ActsExamples::WriterT<StepsVector>(cfg.collection,
                                           "CSVPropgationStepsWriter", level),
        m_cfg(cfg) {
    if (m_cfg.collection.empty()) {
      throw std::invalid_argument("Missing input collection");
    }
  }

  /// Virtual destructor
  ~CsvPropagationStepsWriter() override = default;

  /// End-of-run hook
  ActsExamples::ProcessCode endRun() override {
    return ActsExamples::ProcessCode::SUCCESS;
  }

 private:
  Config m_cfg;  ///!< Internal configuration represenation

 protected:
  /// This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ActsExamples::ProcessCode writeT(
      const ActsExamples::AlgorithmContext& context,
      const StepsVector& stepCollection) override;
};

}  // namespace ActsExamples
