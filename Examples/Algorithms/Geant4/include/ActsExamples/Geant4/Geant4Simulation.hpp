// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>
#include <mutex>
#include <string>

#include "G4VUserDetectorConstruction.hh"

class G4RunManager;
class G4VUserPrimaryGeneratorAction;
class G4UserRunAction;
class G4UserEventAction;
class G4UserTrackingAction;
class G4UserSteppingAction;
class G4MagneticField;

namespace ActsExamples {

class Geant4SurfaceMapper;

/// Algorithm to run Geant4 simulation in the ActsExamples framework
///
/// This algorithm can be configured with the standardard Geant4
/// components:
/// - (a) the generator Action
/// - (b) detector construction and magnetic field
/// - (c) the user actions
///
/// In order to run within the ACTS framework, acces to the
/// EventData is provided by a EventStoreRegistry which provides
/// individual slots for the event containers and the store.
///
/// The Geant4Simulation algorithm clears those after processing.
class Geant4Simulation final : public BareAlgorithm {
 public:
  /// Nested configuration struct for the Geant4 simulation
  struct Config {
    // Name of the output collection : hits
    std::string outputSimHits = "";

    // Name of the output collection : initial particles
    std::string outputParticlesInitial = "";

    // Name of the output collection : final particles
    std::string outputParticlesFinal = "";

    /// The G4 run manager
    G4RunManager* runManager = nullptr;

    /// User Action: Primary generator action of the simulation
    G4VUserPrimaryGeneratorAction* primaryGeneratorAction = nullptr;

    /// User Actions: Run
    std::vector<G4UserRunAction*> runActions = {};

    /// User Actions: Event
    std::vector<G4UserEventAction*> eventActions = {};

    /// User Actions: Tracking
    std::vector<G4UserTrackingAction*> trackingActions = {};

    /// User Actions: Stepping
    std::vector<G4UserSteppingAction*> steppingActions = {};

    /// Detector construction object.
    G4VUserDetectorConstruction* detectorConstruction = nullptr;

    /// The ACTS Magnetic field provider
    G4MagneticField* magneticField = nullptr;
  };

  /// Constructor with arguments
  ///
  /// @param config is the configuration struct
  /// @param level is the logging level to be used
  Geant4Simulation(const Config& config,
                   Acts::Logging::Level level = Acts::Logging::INFO);
  ~Geant4Simulation();

  /// Algorithm execute method, called once per event with context
  ///
  /// @param ctx the AlgorithmContext for this event
  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& ctx) const final override;

  /// Readonly access to the configuration
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  // Has to be mutable; algorithm interface enforces object constness
  mutable std::mutex m_runManagerLock;
};

}  // namespace ActsExamples
