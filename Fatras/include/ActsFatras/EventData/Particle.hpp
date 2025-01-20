// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/ParticleOutcome.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <cmath>
#include <iosfwd>
#include <optional>

namespace ActsFatras {

/// Particle identity information and kinematic state.
///
/// Also stores some simulation-specific properties.
class Particle {
 public:
  /// Construct a default particle with invalid identity.
  Particle() = default;
  /// Construct a particle at rest with explicit mass and charge.
  ///
  /// @param particleId Particle identifier within an event
  /// @param pdg PDG id
  /// @param charge Particle charge in native units
  /// @param mass Particle mass in native units
  ///
  /// @warning It is the users responsibility that charge and mass match
  ///          the PDG particle number.
  Particle(Barcode particleId, Acts::PdgParticle pdg, double charge,
           double mass)
      : m_particleId(particleId), m_pdg(pdg), m_charge(charge), m_mass(mass) {}
  /// Construct a particle at rest from a PDG particle number.
  ///
  /// @param particleId Particle identifier within an event
  /// @param pdg PDG particle number
  ///
  /// Charge and mass are retrieved from the particle data table.
  Particle(Barcode particleId, Acts::PdgParticle pdg);
  Particle(const Particle &) = default;
  Particle(Particle &&) = default;
  Particle &operator=(const Particle &) = default;
  Particle &operator=(Particle &&) = default;

  /// Construct a new particle with a new identifier but same kinematics.
  ///
  /// @note This is intentionally not a regular setter. The particle id
  ///       is used to identify the whole particle. Setting it on an existing
  ///       particle is usually a mistake.
  Particle withParticleId(Barcode particleId) const {
    Particle p = *this;
    p.m_particleId = particleId;
    return p;
  }

  /// Set the process type that generated this particle.
  Particle &setProcess(ProcessType proc) {
    m_process = proc;
    return *this;
  }
  /// Set the pdg.
  Particle setPdg(Acts::PdgParticle pdg) {
    m_pdg = pdg;
    return *this;
  }
  /// Set the charge.
  Particle setCharge(double charge) {
    m_charge = charge;
    return *this;
  }
  /// Set the mass.
  Particle setMass(double mass) {
    m_mass = mass;
    return *this;
  }
  /// Set the particle ID.
  Particle &setParticleId(Barcode barcode) {
    m_particleId = barcode;
    return *this;
  }
  /// Set the space-time position four-vector.
  Particle &setPosition4(const Acts::Vector4 &pos4) {
    m_position4 = pos4;
    return *this;
  }
  /// Set the space-time position four-vector from three-position and time.
  Particle &setPosition4(const Acts::Vector3 &position, double time) {
    m_position4.segment<3>(Acts::ePos0) = position;
    m_position4[Acts::eTime] = time;
    return *this;
  }
  /// Set the space-time position four-vector from scalar components.
  Particle &setPosition4(double x, double y, double z, double time) {
    m_position4[Acts::ePos0] = x;
    m_position4[Acts::ePos1] = y;
    m_position4[Acts::ePos2] = z;
    m_position4[Acts::eTime] = time;
    return *this;
  }
  /// Set the direction three-vector
  Particle &setDirection(const Acts::Vector3 &direction) {
    m_direction = direction;
    m_direction.normalize();
    return *this;
  }
  /// Set the direction three-vector from scalar components.
  Particle &setDirection(double dx, double dy, double dz) {
    m_direction[Acts::ePos0] = dx;
    m_direction[Acts::ePos1] = dy;
    m_direction[Acts::ePos2] = dz;
    m_direction.normalize();
    return *this;
  }
  /// Set the absolute momentum.
  Particle &setAbsoluteMomentum(double absMomentum) {
    m_absMomentum = absMomentum;
    return *this;
  }

  /// Change the energy by the given amount.
  ///
  /// Energy loss corresponds to a negative change. If the updated energy
  /// would result in an unphysical value, the particle is put to rest, i.e.
  /// its absolute momentum is set to zero.
  Particle &correctEnergy(double delta) {
    const auto newEnergy = std::hypot(m_mass, m_absMomentum) + delta;
    if (newEnergy <= m_mass) {
      m_absMomentum = 0.;
    } else {
      m_absMomentum = std::sqrt(newEnergy * newEnergy - m_mass * m_mass);
    }
    return *this;
  }

  /// Particle identifier within an event.
  Barcode particleId() const { return m_particleId; }
  /// Which type of process generated this particle.
  ProcessType process() const { return m_process; }
  /// PDG particle number that identifies the type.
  Acts::PdgParticle pdg() const { return m_pdg; }
  /// Absolute PDG particle number that identifies the type.
  Acts::PdgParticle absolutePdg() const {
    return Acts::makeAbsolutePdgParticle(pdg());
  }
  /// Particle charge.
  double charge() const { return m_charge; }
  /// Particle absolute charge.
  double absoluteCharge() const { return std::abs(m_charge); }
  /// Particle mass.
  double mass() const { return m_mass; }

  /// Particle hypothesis.
  Acts::ParticleHypothesis hypothesis() const {
    return Acts::ParticleHypothesis(
        absolutePdg(), mass(),
        Acts::AnyCharge{static_cast<float>(absoluteCharge())});
  }
  /// Particl qOverP.
  double qOverP() const {
    return hypothesis().qOverP(absoluteMomentum(), charge());
  }

  /// Space-time position four-vector.
  const Acts::Vector4 &fourPosition() const { return m_position4; }
  /// Three-position, i.e. spatial coordinates without the time.
  auto position() const { return m_position4.segment<3>(Acts::ePos0); }
  /// Time coordinate.
  double time() const { return m_position4[Acts::eTime]; }
  /// Energy-momentum four-vector.
  Acts::Vector4 fourMomentum() const {
    Acts::Vector4 mom4;
    // stored direction is always normalized
    mom4[Acts::eMom0] = m_absMomentum * m_direction[Acts::ePos0];
    mom4[Acts::eMom1] = m_absMomentum * m_direction[Acts::ePos1];
    mom4[Acts::eMom2] = m_absMomentum * m_direction[Acts::ePos2];
    mom4[Acts::eEnergy] = energy();
    return mom4;
  }
  /// Unit three-direction, i.e. the normalized momentum three-vector.
  const Acts::Vector3 &direction() const { return m_direction; }
  /// Polar angle.
  double theta() const { return Acts::VectorHelpers::theta(direction()); }
  /// Azimuthal angle.
  double phi() const { return Acts::VectorHelpers::phi(direction()); }
  /// Absolute momentum in the x-y plane.
  double transverseMomentum() const {
    return m_absMomentum * m_direction.segment<2>(Acts::eMom0).norm();
  }
  /// Absolute momentum.
  double absoluteMomentum() const { return m_absMomentum; }
  /// Absolute momentum.
  Acts::Vector3 momentum() const { return absoluteMomentum() * direction(); }
  /// Total energy, i.e. norm of the four-momentum.
  double energy() const { return std::hypot(m_mass, m_absMomentum); }

  /// Check if the particle is alive, i.e. is not at rest.
  bool isAlive() const { return 0. < m_absMomentum; }

  /// Check if this is a secondary particle.
  bool isSecondary() const {
    return particleId().vertexSecondary() != 0 ||
           particleId().generation() != 0 || particleId().subParticle() != 0;
  }

  // simulation specific properties

  /// Set the proper time in the particle rest frame.
  ///
  /// @param properTime passed proper time in the rest frame
  Particle &setProperTime(double properTime) {
    m_properTime = properTime;
    return *this;
  }
  /// Proper time in the particle rest frame.
  double properTime() const { return m_properTime; }

  /// Set the accumulated material measured in radiation/interaction lengths.
  ///
  /// @param pathInX0 accumulated material measured in radiation lengths
  /// @param pathInL0 accumulated material measured in interaction lengths
  Particle &setMaterialPassed(double pathInX0, double pathInL0) {
    m_pathInX0 = pathInX0;
    m_pathInL0 = pathInL0;
    return *this;
  }
  /// Accumulated path within material measured in radiation lengths.
  double pathInX0() const { return m_pathInX0; }
  /// Accumulated path within material measured in interaction lengths.
  double pathInL0() const { return m_pathInL0; }

  /// Set the reference surface.
  ///
  /// @param surface reference surface
  Particle &setReferenceSurface(const Acts::Surface *surface) {
    m_referenceSurface = surface;
    return *this;
  }

  /// Reference surface.
  const Acts::Surface *referenceSurface() const { return m_referenceSurface; }

  /// Check if the particle has a reference surface.
  bool hasReferenceSurface() const { return m_referenceSurface != nullptr; }

  /// Bound track parameters.
  Acts::Result<Acts::BoundTrackParameters> boundParameters(
      const Acts::GeometryContext &gctx) const {
    if (!hasReferenceSurface()) {
      return Acts::Result<Acts::BoundTrackParameters>::failure(
          std::error_code());
    }
    Acts::Result<Acts::Vector2> localResult =
        m_referenceSurface->globalToLocal(gctx, position(), direction());
    if (!localResult.ok()) {
      return localResult.error();
    }
    Acts::BoundVector params;
    params << localResult.value(), phi(), theta(), qOverP(), time();
    return Acts::BoundTrackParameters(referenceSurface()->getSharedPtr(),
                                      params, std::nullopt, hypothesis());
  }

  Acts::CurvilinearTrackParameters curvilinearParameters() const {
    return Acts::CurvilinearTrackParameters(
        fourPosition(), direction(), qOverP(), std::nullopt, hypothesis());
  }

  /// Set the number of hits.
  ///
  /// @param nHits number of hits
  Particle &setNumberOfHits(std::uint32_t nHits) {
    m_numberOfHits = nHits;
    return *this;
  }

  /// Number of hits.
  std::uint32_t numberOfHits() const { return m_numberOfHits; }

  /// Set the outcome of particle.
  ///
  /// @param outcome outcome code
  Particle &setOutcome(ParticleOutcome outcome) {
    m_outcome = outcome;
    return *this;
  }

  /// Particle outcome.
  ParticleOutcome outcome() const { return m_outcome; }

 private:
  // identity, i.e. things that do not change over the particle lifetime.
  /// Particle identifier within the event.
  Barcode m_particleId;
  /// Process type specifier.
  ProcessType m_process = ProcessType::eUndefined;
  /// PDG particle number.
  Acts::PdgParticle m_pdg = Acts::PdgParticle::eInvalid;
  // Particle charge and mass.
  double m_charge = 0.;
  double m_mass = 0.;
  // kinematics, i.e. things that change over the particle lifetime.
  Acts::Vector3 m_direction = Acts::Vector3::UnitZ();
  double m_absMomentum = 0.;
  Acts::Vector4 m_position4 = Acts::Vector4::Zero();
  /// proper time in the particle rest frame
  double m_properTime = 0.;
  // accumulated material
  double m_pathInX0 = 0.;
  double m_pathInL0 = 0.;
  /// number of hits
  std::uint32_t m_numberOfHits = 0;
  /// reference surface
  const Acts::Surface *m_referenceSurface{nullptr};
  /// outcome
  ParticleOutcome m_outcome = ParticleOutcome::Alive;
};

std::ostream &operator<<(std::ostream &os, const Particle &particle);

}  // namespace ActsFatras
