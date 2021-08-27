// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Generators/EventGenerator.hpp"

#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <algorithm>
#include <cstdint>
#include <stdexcept>

ActsExamples::EventGenerator::EventGenerator(const Config& cfg,
                                             Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("EventGenerator", lvl)) {
  if (m_cfg.outputParticles.empty()) {
    throw std::invalid_argument("Missing output particles collection");
  }
  if (m_cfg.generators.empty()) {
    throw std::invalid_argument("No generators are configured");
  }
  if (!m_cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers service");
  }
}

std::string ActsExamples::EventGenerator::name() const {
  return "EventGenerator";
}

std::pair<size_t, size_t> ActsExamples::EventGenerator::availableEvents()
    const {
  return {0u, SIZE_MAX};
}

ActsExamples::ProcessCode ActsExamples::EventGenerator::read(
    const AlgorithmContext& ctx) {
  SimParticleContainer particles;

  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);

  size_t nPrimaryVertices = 0;
  for (size_t iGenerate = 0; iGenerate < m_cfg.generators.size(); ++iGenerate) {
    auto& generate = m_cfg.generators[iGenerate];

    // generate the primary vertices from this generator
    for (size_t n = (*generate.multiplicity)(rng); 0 < n; --n) {
      nPrimaryVertices += 1;

      // generate primary vertex position
      auto vertexPosition = (*generate.vertex)(rng);
      // generate particles associated to this vertex
      auto vertexParticles = (*generate.particles)(rng);

      ACTS_VERBOSE("Generate vertex at " << vertexPosition.transpose());

      auto updateParticleInPlace = [&](ActsFatras::Particle& particle) {
        // only set the primary vertex, leave everything else as-is
        // using the number of primary vertices as the index ensures
        // that barcode=0 is not used, since it is used elsewhere
        // to signify elements w/o an associated particle.
        const auto pid = ActsFatras::Barcode(particle.particleId())
                             .setVertexPrimary(nPrimaryVertices);
        // move particle to the vertex
        const auto pos4 = (vertexPosition + particle.fourPosition()).eval();
        ACTS_VERBOSE(" - particle at " << pos4.transpose());
        // `withParticleId` returns a copy because it changes the identity
        particle = particle.withParticleId(pid).setPosition4(pos4);
      };
      for (auto& vertexParticle : vertexParticles) {
        updateParticleInPlace(vertexParticle);
      }

      ACTS_VERBOSE("event=" << ctx.eventNumber << " generator=" << iGenerate
                            << " primary_vertex=" << nPrimaryVertices
                            << " n_particles=" << vertexParticles.size());

      particles.merge(std::move(vertexParticles));
    }
  }

  ACTS_DEBUG("event=" << ctx.eventNumber
                      << " n_primary_vertices=" << nPrimaryVertices
                      << " n_particles=" << particles.size());

  // move generated event to the store
  ctx.eventStore.add(m_cfg.outputParticles, std::move(particles));
  return ProcessCode::SUCCESS;
}
