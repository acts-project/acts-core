// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Generators/ParametricParticleGenerator.hpp"

#include "ActsFatras/Utilities/ParticleData.hpp"

#include <limits>
#include <random>

ActsExamples::ParametricParticleGenerator::ParametricParticleGenerator(
    const Config& cfg)
    : m_cfg(cfg),
      m_charge(ActsFatras::findCharge(m_cfg.pdg)),
      m_mass(ActsFatras::findMass(m_cfg.pdg)),
      // since we want to draw the direction uniform on the unit sphere, we must
      // draw from cos(theta) instead of theta. see e.g.
      // https://mathworld.wolfram.com/SpherePointPicking.html
      m_cosThetaMin(std::cos(m_cfg.thetaMin)),
      // ensure upper bound is included. see e.g.
      // https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
      m_cosThetaMax(std::nextafter(std::cos(m_cfg.thetaMax),
                                   std::numeric_limits<double>::max())) {}

ActsExamples::SimParticleContainer
ActsExamples::ParametricParticleGenerator::operator()(RandomEngine& rng) {
  using UniformIndex = std::uniform_int_distribution<unsigned int>;
  using UniformReal = std::uniform_real_distribution<double>;

  // choose between particle/anti-particle if requested
  // the upper limit of the distribution is inclusive
  UniformIndex particleTypeChoice(0u, m_cfg.randomizeCharge ? 1u : 0u);
  // (anti-)particle choice is one random draw but defines two properties
  const Acts::PdgParticle pdgChoices[] = {
      m_cfg.pdg,
      static_cast<Acts::PdgParticle>(-m_cfg.pdg),
  };
  const double qChoices[] = {
      m_charge,
      -m_charge,
  };
  UniformReal phiDist(m_cfg.phiMin, m_cfg.phiMax);
  UniformReal cosThetaDist(m_cosThetaMin, m_cosThetaMax);
  UniformReal pDist(m_cfg.pMin, m_cfg.pMax);

  SimParticleContainer particles;
  particles.reserve(m_cfg.numParticles);

  // counter will be reused as barcode particle number which must be non-zero.
  for (size_t ip = 1; ip <= m_cfg.numParticles; ++ip) {
    // all particles are treated as originating from the same primary vertex
    const auto pid = ActsFatras::Barcode(0u).setParticle(ip);

    // draw parameters
    const unsigned int type = particleTypeChoice(rng);
    const Acts::PdgParticle pdg = pdgChoices[type];
    const double q = qChoices[type];
    const double phi = phiDist(rng);
    const double cosTheta = cosThetaDist(rng);
    const double sinTheta = std::sqrt(1 - cosTheta * cosTheta);
    double p = pDist(rng);

    // we already have sin/cos theta. they can be used directly to
    Acts::Vector3 dir;
    dir[Acts::eMom0] = sinTheta * std::cos(phi);
    dir[Acts::eMom1] = sinTheta * std::sin(phi);
    dir[Acts::eMom2] = cosTheta;

    // construct the particle;
    ActsFatras::Particle particle(pid, pdg, q, m_mass);
    particle.setDirection(dir);
    p *= m_cfg.pTransverse ? 1. / sinTheta : 1.;
    particle.setAbsoluteMomentum(p);

    // generated particle ids are already ordered and should end up at the end
    particles.insert(particles.end(), std::move(particle));
  }

  return particles;
}
