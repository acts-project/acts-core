// This file is part of the Acts project.
//
// Copyright (C) 2017-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>
#include <string>
#include <unordered_map>

namespace Acts {
class Surface;
}

namespace ActsExamples {

class PropagatorInterface;
struct AlgorithmContext;

struct PropagationSummary {
  explicit PropagationSummary(Acts::BoundTrackParameters startParameters_)
      : startParameters(std::move(startParameters_)) {}

  /// The start parameters
  Acts::BoundTrackParameters startParameters;

  /// Totoal number of successful steps
  std::size_t nSteps = 0;

  /// Totoal number of attempted steps
  std::size_t nStepTrials = 0;

  /// Path length
  double pathLength = 0;

  /// Steps
  std::vector<Acts::detail::Step> steps;
};

using PropagationSummaries = std::vector<PropagationSummary>;

/// Using some short hands for Recorded Material
using RecordedMaterial = Acts::MaterialInteractor::result_type;

/// And recorded material track
/// - this is start:  position, start momentum
///   and the Recorded material
using RecordedMaterialTrack =
    std::pair<std::pair<Acts::Vector3, Acts::Vector3>, Acts::RecordedMaterial>;

/// Finally the output of the propagation test
using PropagationOutput = std::pair<PropagationSummary, RecordedMaterial>;

/// @brief this test algorithm performs test propagation
/// within the Acts::Propagator
///
/// If the propagator is equipped appropriately, it can
/// also be used to test the Extrapolator within the geomtetry
class PropagationAlgorithm : public IAlgorithm {
 public:
  struct Config {
    /// The output step collection
    std::string outputSummaryCollection = "propagation_summary";
    /// The output material collection
    std::string outputMaterialCollection = "recorded_material";

    /// Instance of a propagator wrapper that performs the actual propagation
    std::shared_ptr<PropagatorInterface> propagatorImpl = nullptr;
    /// Switch the logger to sterile - for timing measurements
    bool sterileLogger = false;
    /// debug output
    bool debugOutput = false;
    /// Modify the behavior of the material interaction: energy loss
    bool energyLoss = true;
    /// Modify the behavior of the material interaction: scattering
    bool multipleScattering = true;
    /// Modify the behavior of the material interaction: record
    bool recordMaterialInteractions = true;
    /// looper protection
    double ptLoopers = 500 * Acts::UnitConstants::MeV;
    /// Max step size steering
    double maxStepSize = 5 * Acts::UnitConstants::m;
    /// Switch covariance transport on
    bool covarianceTransport = false;
    /// Input track parameters
    std::string inputTrackParameters = "InputTrackParameters";
    /// The step collection to be stored
    std::string outputSummarySollection = "PropagationSummary";
    /// The material collection to be stored
    std::string outputSummaryCollection = "RecordedMaterialTracks";
  };

  /// Constructor
  /// @param [in] config is the configuration struct
  /// @param [in] loglevel is the logging level
  PropagationAlgorithm(const Config& config, Acts::Logging::Level level);

  /// Framework execute method
  /// @param [in] the algorithm context for event consistency
  /// @return is a process code indicating success or not
  ActsExamples::ProcessCode execute(
      const AlgorithmContext& context) const override;

  /// Get const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<TrackParametersContainer> m_inputTrackParameters{
      this, "InputTrackParameters"};

  WriteDataHandle<PropagationSummaries> m_outputSummary{this, "OutputSummary"};


  WriteDataHandle<std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>>
      m_outputMaterialTracks{this, "RecordedMaterial"};
};

}  // namespace ActsExamples
