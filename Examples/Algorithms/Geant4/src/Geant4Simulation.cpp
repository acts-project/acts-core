
// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/Geant4Simulation.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geant4/EventStoreRegistry.hpp"

#include <iostream>
#include <stdexcept>

#include <G4FieldManager.hh>
#include <G4MagneticField.hh>
#include <G4RunManager.hh>
#include <G4TransportationManager.hh>
#include <G4UniformMagField.hh>
#include <G4UserEventAction.hh>
#include <G4UserLimits.hh>
#include <G4UserRunAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4UserTrackingAction.hh>
#include <G4VUserDetectorConstruction.hh>

ActsExamples::Geant4Simulation::Geant4Simulation(
    const ActsExamples::Geant4Simulation::Config& config,
    Acts::Logging::Level level)
    : BareAlgorithm("Geant4Simulation", level), m_cfg(config) {
  if (m_cfg.runManager == nullptr) {
    throw std::invalid_argument("Missing G4 RunManager object");
  }
  if (m_cfg.detectorConstruction == nullptr) {
    throw std::invalid_argument("Missing G4 DetectorConstruction object");
  }
  if (m_cfg.primaryGeneratorAction == nullptr) {
    throw std::invalid_argument("Missing G4 PrimaryGeneratorAction object");
  }

  m_cfg.runManager->SetUserInitialization(m_cfg.detectorConstruction);
  // Set the primary generator action
  m_cfg.runManager->SetUserAction(m_cfg.primaryGeneratorAction);
  if (m_cfg.runAction != nullptr) {
    m_cfg.runManager->SetUserAction(m_cfg.runAction);
  }
  // Set the user actions
  if (m_cfg.runAction != nullptr) {
    m_cfg.runManager->SetUserAction(m_cfg.runAction);
  }
  if (m_cfg.eventAction != nullptr) {
    m_cfg.runManager->SetUserAction(m_cfg.eventAction);
  }
  if (m_cfg.trackingAction != nullptr) {
    m_cfg.runManager->SetUserAction(m_cfg.trackingAction);
  }
  if (m_cfg.steppingAction != nullptr) {
    m_cfg.runManager->SetUserAction(m_cfg.steppingAction);
  }

  // Initialize the Geant4 run manager
  m_cfg.runManager->Initialize();

  // Please note:
  // The following two blocks rely on the fact that the Acts
  // detector constructions cache the world volume

  // Set the magnetic field
  if (m_cfg.magneticField != nullptr) {
    ACTS_INFO("Setting ACTS configured field to Geant4.");
    // Get the g4World cache
    G4VPhysicalVolume* g4World = m_cfg.detectorConstruction->Construct();
    /// Set the field ot the G4Field manager
    G4FieldManager* fieldMgr = new G4FieldManager();
    fieldMgr->SetDetectorField(m_cfg.magneticField);
    fieldMgr->CreateChordFinder(m_cfg.magneticField);

    // Propagate down to all childrend
    g4World->GetLogicalVolume()->SetFieldManager(fieldMgr, true);
  }
}

ActsExamples::Geant4Simulation::~Geant4Simulation() {}

ActsExamples::ProcessCode ActsExamples::Geant4Simulation::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Ensure exclusive access to the geant run manager
  std::lock_guard<std::mutex> guard(m_runManagerLock);

  // Register the current event store to the registry
  // this will allow access from the User*Actions
  EventStoreRegistry::boards[ctx.eventNumber] = &(ctx.eventStore);

  // Start simulation. each track is simulated as a separate Geant4 event.
  m_cfg.runManager->BeamOn(1);

  // Output handling: Initial/Final particles
  if (not m_cfg.outputParticlesInitial.empty() and
      not m_cfg.outputParticlesFinal.empty()) {
    // Initial state of partciles
    SimParticleContainer outputParticlesInitial;
    outputParticlesInitial.insert(
        EventStoreRegistry::particlesInitial[ctx.eventNumber].begin(),
        EventStoreRegistry::particlesInitial[ctx.eventNumber].end());
    EventStoreRegistry::particlesInitial[ctx.eventNumber].clear();
    // Register to the event store
    ctx.eventStore.add(m_cfg.outputParticlesInitial,
                       std::move(outputParticlesInitial));
    // Final state of partciles
    SimParticleContainer outputParticlesFinal;
    outputParticlesFinal.insert(
        EventStoreRegistry::particlesFinal[ctx.eventNumber].begin(),
        EventStoreRegistry::particlesFinal[ctx.eventNumber].end());
    EventStoreRegistry::particlesFinal[ctx.eventNumber].clear();
    // Register to the event store
    ctx.eventStore.add(m_cfg.outputParticlesFinal,
                       std::move(outputParticlesFinal));
  }

  // Output handling: Simulated hits
  if (not m_cfg.outputSimHits.empty()) {
    SimHitContainer simHits;
    simHits.insert(EventStoreRegistry::hits[ctx.eventNumber].begin(),
                   EventStoreRegistry::hits[ctx.eventNumber].end());
    EventStoreRegistry::hits[ctx.eventNumber].clear();
    // Register to the event store
    ctx.eventStore.add(m_cfg.outputSimHits, std::move(simHits));
  }

  return ActsExamples::ProcessCode::SUCCESS;
}
