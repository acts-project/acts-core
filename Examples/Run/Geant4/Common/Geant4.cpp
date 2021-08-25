// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Geant4.hpp"

#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Geant4/EventStoreRegistry.hpp"
#include "ActsExamples/Geant4/G4DetectorConstructionFactory.hpp"
#include "ActsExamples/Geant4/Geant4Simulation.hpp"
#include "ActsExamples/Geant4/MaterialSteppingAction.hpp"
#include "ActsExamples/Geant4/SimParticleTranslation.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/ParticleGunOptions.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <memory>
#include <string>

#include <FTFP_BERT.hh>
#include <G4RunManager.hh>
#include <boost/program_options.hpp>

namespace ActsExamples {

void setupGeant4Simulation(const ActsExamples::Options::Variables& vars,
                           ActsExamples::Sequencer& sequencer,
                           G4VUserDetectorConstruction* detector,
                           G4VUserPrimaryGeneratorAction* generatorAction,
                           std::vector<G4UserRunAction*> runActions,
                           std::vector<G4UserEventAction*> eventActions,
                           std::vector<G4UserTrackingAction*> trackingActions,
                           std::vector<G4UserSteppingAction*> steppingActions) {
  // Create an event Store registry
  EventStoreRegistry esRegistry(vars["events"].as<size_t>());

  // The G4 Run Manager and the physics list go first
  auto g4RunManager = new G4RunManager();
  g4RunManager->SetUserInitialization(new FTFP_BERT);

  // Setup the main Geant4 algorithm
  Geant4Simulation::Config g4Cfg;
  g4Cfg.runManager = g4RunManager;

  g4Cfg.primaryGeneratorAction = generatorAction;
  g4Cfg.detectorConstruction = std::move(detector);

  // Set the user actions
  g4Cfg.runActions = runActions;
  g4Cfg.eventActions = eventActions;
  g4Cfg.trackingActions = trackingActions;
  g4Cfg.steppingActions = steppingActions;

  sequencer.addAlgorithm(std::make_shared<Geant4Simulation>(g4Cfg));

  return;
}

int runGeantinoRecording(
    const ActsExamples::Options::Variables& vars,
    std::shared_ptr<ActsExamples::G4DetectorConstructionFactory>
        g4DetectorFactory) {
  // Set up the sequencer
  Sequencer sequencer(Options::readSequencerConfig(vars));
  // Output level and log level
  auto logLevel = Options::readLogLevel(vars);
  auto outputDir =
      ensureWritableDirectory(vars["output-dir"].as<std::string>());

  // Basic services
  auto rnd =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vars));

  const std::string particles = "particles";
  const std::string outputMaterialTracks = "output_tracks";

  // Event generation w/ particle gun
  EventGenerator::Config evgen = Options::readParticleGunOptions(vars);
  evgen.outputParticles = "particles";
  evgen.randomNumbers = rnd;
  sequencer.addReader(std::make_shared<EventGenerator>(evgen, logLevel));

  // Read the particle from teh generator
  SimParticleTranslation::Config g4PrCfg;
  g4PrCfg.inputParticles = evgen.outputParticles;
  g4PrCfg.forceParticle = true;
  g4PrCfg.forcedMass = 0.;
  g4PrCfg.forcedPdgCode = 999;
  auto particleReader = new SimParticleTranslation(g4PrCfg);

  // Release the detector construction
  G4VUserDetectorConstruction* detector = (*g4DetectorFactory)().release();

  // The Geant4 actions needed
  std::vector<G4UserRunAction*> runActions = {};
  std::vector<G4UserEventAction*> eventActions = {};
  std::vector<G4UserTrackingAction*> trackingActions = {};

  MaterialSteppingAction::Config mStepCfg;
  mStepCfg.excludeMaterials = {"Air", "Vacuum"};
  std::vector<G4UserSteppingAction*> steppingActions = {
      new MaterialSteppingAction(mStepCfg)};

  // Set up the Geant4 Simulation
  setupGeant4Simulation(vars, sequencer, detector, particleReader, runActions,
                        eventActions, trackingActions, steppingActions);

  // setup the output writing
  if (vars["output-root"].as<bool>()) {
    // Write the propagation steps as ROOT TTree
    RootMaterialTrackWriter::Config materialTrackWriter;
    materialTrackWriter.prePostStep = true;
    materialTrackWriter.recalculateTotals = true;
    materialTrackWriter.collection = outputMaterialTracks;
    materialTrackWriter.filePath =
        joinPaths(outputDir, outputMaterialTracks + ".root");
    sequencer.addWriter(std::make_shared<RootMaterialTrackWriter>(
        materialTrackWriter, logLevel));
  }

  return sequencer.run();
}

}  // namespace ActsExamples
