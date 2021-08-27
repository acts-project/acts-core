// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/SimParticleTranslation.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/PdgParticle.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geant4/EventStoreRegistry.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <G4ChargedGeantino.hh>
#include <G4Event.hh>
#include <G4Geantino.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>
#include <G4UnitsTable.hh>
#include <globals.hh>

ActsExamples::SimParticleTranslation::SimParticleTranslation(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : G4VUserPrimaryGeneratorAction(), m_cfg(cfg), m_logger(std::move(logger)) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particle collection");
  }
}

ActsExamples::SimParticleTranslation::~SimParticleTranslation() {}

void ActsExamples::SimParticleTranslation::GeneratePrimaries(G4Event* anEvent) {
  anEvent->SetEventID(m_eventNr++);
  unsigned int eventID = anEvent->GetEventID();

  ACTS_DEBUG("Primary Generator Action for Event: " << eventID);

  WhiteBoard* eventStore = EventStoreRegistry::boards[eventID];
  if (eventStore == nullptr) {
    ACTS_WARNING("No EventStore instance could be found for this event!");
    return;
  }

  // Get the number of input particles
  const auto inputParticles =
      eventStore->get<ActsExamples::SimParticleContainer>(m_cfg.inputParticles);

  // Feserve appropriate resources for initial/final particles
  EventStoreRegistry::particlesInitial[eventID].reserve(inputParticles.size());
  EventStoreRegistry::particlesFinal[eventID].reserve(inputParticles.size());

  // Reserve hopefully enough hit space
  EventStoreRegistry::hits[eventID].reserve(inputParticles.size() *
                                            m_cfg.reserveHitsPerParticle);

  // Default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4PrimaryVertex* pVertex = nullptr;

  // We are looping through the particles and flush per vertex
  std::unique_ptr<SimParticle::Vector4> lastVertex = nullptr;

  constexpr double convertLength = CLHEP::mm / Acts::UnitConstants::mm;
  constexpr double convertEnergy = CLHEP::GeV / Acts::UnitConstants::GeV;

  unsigned int pCounter = 0;
  // Loop over the input partilces and run
  for (const auto& part : inputParticles) {
    auto currentVertex = part.fourPosition();
    if (lastVertex == nullptr or not currentVertex.isApprox(*lastVertex)) {
      // Add the vertex to the event
      if (pVertex != nullptr) {
        anEvent->AddPrimaryVertex(pVertex);
        ACTS_DEBUG("Flushing " << pCounter
                               << " particles associated with vertex "
                               << Acts::toString(*lastVertex));
        pCounter = 0;
      }
      lastVertex = std::make_unique<SimParticle::Vector4>(currentVertex);
      pVertex = new G4PrimaryVertex(
          currentVertex[0] * convertLength, currentVertex[1] * convertLength,
          currentVertex[2] * convertLength, currentVertex[3]);
    }
    // Add a new primary to the vertex
    auto mom4 = part.fourMomentum() * convertEnergy;

    // Particle properties, may be forced to specific value
    G4int particlePdgCode =
        m_cfg.forceParticle ? m_cfg.forcedPdgCode : part.pdg();
    G4double particleMass =
        m_cfg.forceParticle ? part.mass() * convertEnergy : 0.;

    // Check if it is a Geantino / ChargedGeantino
    G4ParticleDefinition* particleDefinition =
        particleTable->FindParticle(particlePdgCode);
    if (particleDefinition == nullptr) {
      switch (particlePdgCode) {
        case 999: {
          particleDefinition = G4Geantino::Definition();
        } break;
        case 998: {
          particleDefinition = G4ChargedGeantino::Definition();
        } break;
        default:
          break;
      }
    }

    // Skip if tranlation failed
    if (particleDefinition == nullptr) {
      ACTS_VERBOSE(
          "Could not translate particle with PDG code : " << particlePdgCode);
      continue;
    }

    ACTS_VERBOSE("Adding particle with name '"
                 << particleDefinition->GetParticleName() << "'.");

    G4PrimaryParticle* particle = new G4PrimaryParticle(particleDefinition);

    particle->SetMass(particleMass);
    particle->Set4Momentum(mom4[0], mom4[1], mom4[2], mom4[3]);
    particle->SetCharge(part.charge());
    particle->SetTrackID(pCounter);
    // Add the primary to the vertex
    pVertex->SetPrimary(particle);
    ++pCounter;
  }
  // Final vertex to be added
  if (pVertex != nullptr) {
    anEvent->AddPrimaryVertex(pVertex);
    ACTS_DEBUG("Flushing " << pCounter << " particles associated with vertex "
                           << Acts::toString(*lastVertex));
  }
}
