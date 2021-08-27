
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

namespace {
/// Helper method to add the user actions
/// @tparam manager_t the run manager type
/// @tparam actions_t the actions iterable list
template <typename manager_t, typename actions_t>
void setUserActions(manager_t& manager, actions_t& actions) {
  for (const auto& action : actions) {
    if (action != nullptr) {
      manager.SetUserAction(action);
    }
  }
}
}  // namespace

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

  // Set the detector construction
  m_cfg.runManager->SetUserInitialization(m_cfg.detectorConstruction);

  // Set the primary generator action
  m_cfg.runManager->SetUserAction(m_cfg.primaryGeneratorAction);

  // Set the configured user actions
  setUserActions(*m_cfg.runManager, m_cfg.runActions);
  setUserActions(*m_cfg.runManager, m_cfg.eventActions);
  setUserActions(*m_cfg.runManager, m_cfg.trackingActions);
  setUserActions(*m_cfg.runManager, m_cfg.steppingActions);

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

  ACTS_DEBUG("Sending Geant RunManager the BeamOn() command.");
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

  // Output handling: Material tracks
  if (not m_cfg.outputMaterialTracks.empty()) {
    ctx.eventStore.add(
        m_cfg.outputMaterialTracks,
        std::move(EventStoreRegistry::recordedMaterial[ctx.eventNumber]));
  }

  // Clear the EventStoreRegistry
  EventStoreRegistry::clearEvent(ctx.eventNumber);

  return ActsExamples::ProcessCode::SUCCESS;
}
