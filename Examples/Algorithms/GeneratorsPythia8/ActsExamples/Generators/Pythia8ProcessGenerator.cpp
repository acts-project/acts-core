// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Generators/Pythia8ProcessGenerator.hpp"

#include <algorithm>
#include <iterator>
#include <random>

#include <Pythia8/Pythia.h>

namespace {
struct FrameworkRndmEngine : public Pythia8::RndmEngine {
  ActsExamples::RandomEngine& rng;

  FrameworkRndmEngine(ActsExamples::RandomEngine& rng_) : rng(rng_) {}
  double flat() {
    return std::uniform_real_distribution<double>(0.0, 1.0)(rng);
  }
};
}  // namespace

ActsExamples::Pythia8Generator::Pythia8Generator(const Config& cfg,
                                                 Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("Pythia8Generator", lvl)),
      m_pythia8(std::make_unique<Pythia8::Pythia>("", false)) {
  // disable all output by default but allow reenable via config
  m_pythia8->settings.flag("Print:quiet", true);
  for (const auto& setting : m_cfg.settings) {
    ACTS_VERBOSE("use Pythia8 setting '" << setting << "'");
    m_pythia8->readString(setting.c_str());
  }
  m_pythia8->settings.mode("Beams:idA", m_cfg.pdgBeam0);
  m_pythia8->settings.mode("Beams:idB", m_cfg.pdgBeam1);
  m_pythia8->settings.mode("Beams:frameType", 1);
  m_pythia8->settings.parm("Beams:eCM",
                           m_cfg.cmsEnergy / Acts::UnitConstants::GeV);
  m_pythia8->init();
}

// needed to allow unique_ptr of forward-declared Pythia class
ActsExamples::Pythia8Generator::~Pythia8Generator() {}

ActsExamples::SimParticleContainer ActsExamples::Pythia8Generator::operator()(
    RandomEngine& rng) {
  using namespace Acts::UnitLiterals;

  SimParticleContainer::sequence_type generated;
  std::vector<SimParticle::Vector4> vertexPositions;

  // pythia8 is not thread safe and generation needs to be protected
  std::lock_guard<std::mutex> lock(m_pythia8Mutex);
  // use per-thread random engine also in pythia
  FrameworkRndmEngine rndmEngine(rng);
  m_pythia8->rndm.rndmEnginePtr(&rndmEngine);
  m_pythia8->next();

  // convert generated final state particles into internal format
  for (int ip = 0; ip < m_pythia8->event.size(); ++ip) {
    const auto& genParticle = m_pythia8->event[ip];

    // ignore beam particles
    if (genParticle.statusHepMC() == 4) {
      continue;
    }
    // only interested in final, visible particles
    if (not genParticle.isFinal()) {
      continue;
    }
    if (not genParticle.isVisible()) {
      continue;
    }

    // production vertex. Pythia8 time uses units mm/c, and we use c=1
    SimParticle::Vector4 pos4(
        genParticle.xProd() * 1_mm, genParticle.yProd() * 1_mm,
        genParticle.zProd() * 1_mm, genParticle.tProd() * 1_mm);

    // define the particle identifier including possible secondary vertices
    ActsFatras::Barcode particleId(0u);
    // ensure particle identifier component is non-zero
    particleId.setParticle(1u + generated.size());
    // only secondaries have a defined vertex position
    if (genParticle.hasVertex()) {
      // either add to existing secondary vertex if exists or create new one
      // TODO can we do this w/o the manual search and position check?
      auto it = std::find_if(
          vertexPositions.begin(), vertexPositions.end(),
          [=](const SimParticle::Vector4& pos) { return (pos == pos4); });
      if (it == vertexPositions.end()) {
        // no matching secondary vertex exists -> create new one
        vertexPositions.emplace_back(pos4);
        particleId.setVertexSecondary(vertexPositions.size());
        ACTS_VERBOSE("created new secondary vertex " << pos4.transpose());
      } else {
        particleId.setVertexSecondary(
            1u + std::distance(vertexPositions.begin(), it));
      }
    }

    // construct internal particle
    const auto pdg = static_cast<Acts::PdgParticle>(genParticle.id());
    const auto charge = genParticle.charge() * 1_e;
    const auto mass = genParticle.m0() * 1_GeV;
    ActsFatras::Particle particle(particleId, pdg, charge, mass);
    particle.setPosition4(pos4);
    // normalization/ units are not import for the direction
    particle.setDirection(genParticle.px(), genParticle.py(), genParticle.pz());
    particle.setAbsoluteMomentum(
        std::hypot(genParticle.px(), genParticle.py(), genParticle.pz()) *
        1_GeV);

    generated.push_back(std::move(particle));
  }

  SimParticleContainer out;
  out.insert(generated.begin(), generated.end());
  return out;
}
