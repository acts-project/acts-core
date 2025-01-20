// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4HepMC/EventRecording.hpp"

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>

#include <FTFP_BERT.hh>
#include <G4RunManager.hh>
#include <G4VUserDetectorConstruction.hh>
#include <HepMC3/GenParticle.h>

#include "EventAction.hpp"
#include "PrimaryGeneratorAction.hpp"
#include "RunAction.hpp"
#include "SteppingAction.hpp"

namespace ActsExamples {

EventRecording::~EventRecording() {
  m_runManager = nullptr;
}

EventRecording::EventRecording(const EventRecording::Config& config,
                               Acts::Logging::Level level)
    : IAlgorithm("EventRecording", level),
      m_cfg(config),
      m_runManager(std::make_unique<G4RunManager>()) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particle collection");
  }
  if (m_cfg.outputHepMcTracks.empty()) {
    throw std::invalid_argument("Missing output event collection");
  }
  if (m_cfg.detector == nullptr) {
    throw std::invalid_argument("Missing detector construction object");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_outputEvents.initialize(m_cfg.outputHepMcTracks);

  // Now set up the Geant4 simulation

  // G4RunManager deals with the lifetime of these objects
  m_runManager->SetUserInitialization(
      m_cfg.detector->buildGeant4DetectorConstruction(m_cfg.constructionOptions)
          .release());
  m_runManager->SetUserInitialization(new FTFP_BERT);
  m_runManager->SetUserAction(new Geant4::HepMC3::RunAction());
  m_runManager->SetUserAction(
      new Geant4::HepMC3::EventAction(m_cfg.processesCombine));
  m_runManager->SetUserAction(
      new Geant4::HepMC3::PrimaryGeneratorAction(m_cfg.seed1, m_cfg.seed2));
  m_runManager->SetUserAction(
      new Geant4::HepMC3::SteppingAction(m_cfg.processesReject));
  m_runManager->Initialize();
}

ProcessCode EventRecording::execute(const AlgorithmContext& context) const {
  // ensure exclusive access to the geant run manager
  std::lock_guard<std::mutex> guard(m_runManagerLock);

  // Retrieve the initial particles
  const auto initialParticles = m_inputParticles(context);

  // Storage of events that will be produced
  std::vector<HepMC3::GenEvent> events;
  events.reserve(initialParticles.size());

  for (const auto& part : initialParticles) {
    // Prepare the particle gun
    Geant4::HepMC3::PrimaryGeneratorAction::instance()->prepareParticleGun(
        part);

    // Begin with the simulation
    // MARK: fpeMaskBegin(FLTINV, 1, #4021)
    m_runManager->BeamOn(1);  // this is where the actual FPE occurs

    // Test if the event was aborted
    if (Geant4::HepMC3::SteppingAction::instance()
            ->eventAborted()) {  // this is where the FPE is attributed to due
                                 // to inlining
      continue;
    }
    // MARK: fpeMaskEnd(FLTINV)

    // Set event start time
    HepMC3::GenEvent event = Geant4::HepMC3::EventAction::instance()->event();
    HepMC3::FourVector shift(0., 0., 0., part.time() / Acts::UnitConstants::mm);
    event.shift_position_by(shift);

    // Set beam particle properties
    const Acts::Vector4 momentum4 =
        part.fourMomentum() / Acts::UnitConstants::GeV;
    HepMC3::FourVector beamMom4(momentum4[0], momentum4[1], momentum4[2],
                                momentum4[3]);
    auto beamParticle = event.particles()[0];
    beamParticle->set_momentum(beamMom4);
    beamParticle->set_pid(part.pdg());

    if (m_cfg.processSelect.empty()) {
      // Store the result
      events.push_back(std::move(event));
    } else {
      bool storeEvent = false;
      // Test if the event has a process of interest in it
      for (const auto& vertex : event.vertices()) {
        if (vertex->id() == -1) {
          vertex->add_particle_in(beamParticle);
        }
        const std::vector<std::string> vertexAttributes =
            vertex->attribute_names();
        for (const auto& att : vertexAttributes) {
          if ((vertex->attribute_as_string(att).find(m_cfg.processSelect) !=
               std::string::npos) &&
              !vertex->particles_in().empty() &&
              vertex->particles_in()[0]->attribute<HepMC3::IntAttribute>(
                  "TrackID") &&
              vertex->particles_in()[0]
                      ->attribute<HepMC3::IntAttribute>("TrackID")
                      ->value() == 1) {
            storeEvent = true;
            break;
          }
        }
        if (storeEvent) {
          break;
        }
      }
      // Store the result
      if (storeEvent) {
        // Remove vertices w/o outgoing particles and particles w/o production
        // vertices
        while (true) {
          bool sane = true;
          for (const auto& v : event.vertices()) {
            if (!v) {
              continue;
            }
            if (v->particles_out().empty()) {
              event.remove_vertex(v);
              sane = false;
            }
          }
          for (const auto& p : event.particles()) {
            if (!p) {
              continue;
            }
            if (!p->production_vertex()) {
              event.remove_particle(p);
              sane = false;
            }
          }
          if (sane) {
            break;
          }
        }
        events.push_back(std::move(event));
      }
    }
  }

  ACTS_INFO(initialParticles.size() << " initial particles provided");
  ACTS_INFO(events.size() << " tracks generated");

  // Write the recorded material to the event store
  m_outputEvents(context, std::move(events));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
