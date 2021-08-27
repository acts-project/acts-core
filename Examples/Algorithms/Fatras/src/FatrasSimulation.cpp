// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Fatras/FatrasSimulation.hpp"

#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsFatras/Kernel/InteractionList.hpp"
#include "ActsFatras/Kernel/Simulation.hpp"
#include "ActsFatras/Physics/Decay/NoDecay.hpp"
#include "ActsFatras/Physics/ElectroMagnetic/PhotonConversion.hpp"
#include "ActsFatras/Physics/StandardInteractions.hpp"
#include "ActsFatras/Selectors/KinematicCasts.hpp"
#include "ActsFatras/Selectors/SelectorHelpers.hpp"
#include "ActsFatras/Selectors/SurfaceSelectors.hpp"

namespace {

/// Simple struct to select surfaces where hits should be generated.
struct HitSurfaceSelector {
  bool sensitive = false;
  bool material = false;
  bool passive = false;

  /// Check if the surface should be used.
  bool operator()(const Acts::Surface &surface) const {
    // sensitive/material are not mutually exclusive
    bool isSensitive = surface.associatedDetectorElement();
    bool isMaterial = surface.surfaceMaterial();
    // passive should be an orthogonal category
    bool isPassive = not(isSensitive or isMaterial);
    return (isSensitive and sensitive) or (isMaterial and material) or
           (isPassive and passive);
  }
};

}  // namespace

// Same interface as `ActsFatras::Simulation` but with concrete types.
struct ActsExamples::detail::FatrasSimulation {
  virtual ~FatrasSimulation() = default;
  virtual Acts::Result<std::vector<ActsFatras::FailedParticle>> simulate(
      const Acts::GeometryContext &, const Acts::MagneticFieldContext &,
      ActsExamples::RandomEngine &, const ActsExamples::SimParticleContainer &,
      ActsExamples::SimParticleContainer::sequence_type &,
      ActsExamples::SimParticleContainer::sequence_type &,
      ActsExamples::SimHitContainer::sequence_type &) const = 0;
};

namespace {

// Magnetic-field specific PIMPL implementation.
//
// This always uses the EigenStepper with default extensions for charged
// particle propagation and is thus limited to propagation in vacuum at the
// moment.
// @TODO: Remove this, unneeded after #675
struct FatrasSimulationT final : ActsExamples::detail::FatrasSimulation {
  using CutPMin = ActsFatras::Min<ActsFatras::Casts::P>;

  // typedefs for charge particle simulation
  // propagate charged particles numerically in the given magnetic field
  using ChargedStepper = Acts::EigenStepper<>;
  using ChargedPropagator = Acts::Propagator<ChargedStepper, Acts::Navigator>;
  // charged particles w/ standard em physics list and selectable hits
  using ChargedSelector = CutPMin;
  using ChargedSimulation = ActsFatras::SingleParticleSimulation<
      ChargedPropagator, ActsFatras::StandardChargedElectroMagneticInteractions,
      HitSurfaceSelector, ActsFatras::NoDecay>;

  // typedefs for neutral particle simulation
  // propagate neutral particles with just straight lines
  using NeutralStepper = Acts::StraightLineStepper;
  using NeutralPropagator = Acts::Propagator<NeutralStepper, Acts::Navigator>;
  // neutral particles w/ photon conversion and no hits
  using NeutralSelector = CutPMin;
  using NeutralInteractions =
      ActsFatras::InteractionList<ActsFatras::PhotonConversion>;
  using NeutralSimulation = ActsFatras::SingleParticleSimulation<
      NeutralPropagator, NeutralInteractions, ActsFatras::NoSurface,
      ActsFatras::NoDecay>;

  // combined simulation type
  using Simulation = ActsFatras::Simulation<ChargedSelector, ChargedSimulation,
                                            NeutralSelector, NeutralSimulation>;

  Simulation simulation;

  FatrasSimulationT(const ActsExamples::FatrasSimulation::Config &cfg,
                    Acts::Logging::Level lvl)
      : simulation(
            ChargedSimulation(
                ChargedPropagator(ChargedStepper(cfg.magneticField),
                                  Acts::Navigator{{cfg.trackingGeometry}}),
                Acts::getDefaultLogger("Simulation", lvl)),
            NeutralSimulation(
                NeutralPropagator(NeutralStepper(),
                                  Acts::Navigator{{cfg.trackingGeometry}}),
                Acts::getDefaultLogger("Simulation", lvl))) {
    using namespace ActsFatras;
    using namespace ActsFatras::detail;
    // apply the configuration

    // minimal p cut on input particles and as is-alive check for interactions
    simulation.selectCharged.valMin = cfg.pMin;
    simulation.selectNeutral.valMin = cfg.pMin;
    simulation.charged.interactions =
        makeStandardChargedElectroMagneticInteractions(cfg.pMin);

    // processes are enabled by default
    if (not cfg.emScattering) {
      simulation.charged.interactions.template disable<StandardScattering>();
    }
    if (not cfg.emEnergyLossIonisation) {
      simulation.charged.interactions.template disable<StandardBetheBloch>();
    }
    if (not cfg.emEnergyLossRadiation) {
      simulation.charged.interactions.template disable<StandardBetheHeitler>();
    }
    if (not cfg.emPhotonConversion) {
      simulation.neutral.interactions.template disable<PhotonConversion>();
    }

    // configure hit surfaces for charged particles
    simulation.charged.selectHitSurface.sensitive = cfg.generateHitsOnSensitive;
    simulation.charged.selectHitSurface.material = cfg.generateHitsOnMaterial;
    simulation.charged.selectHitSurface.passive = cfg.generateHitsOnPassive;
  }
  ~FatrasSimulationT() final override = default;

  Acts::Result<std::vector<ActsFatras::FailedParticle>> simulate(
      const Acts::GeometryContext &geoCtx,
      const Acts::MagneticFieldContext &magCtx, ActsExamples::RandomEngine &rng,
      const ActsExamples::SimParticleContainer &inputParticles,
      ActsExamples::SimParticleContainer::sequence_type
          &simulatedParticlesInitial,
      ActsExamples::SimParticleContainer::sequence_type
          &simulatedParticlesFinal,
      ActsExamples::SimHitContainer::sequence_type &simHits)
      const final override {
    return simulation.simulate(geoCtx, magCtx, rng, inputParticles,
                               simulatedParticlesInitial,
                               simulatedParticlesFinal, simHits);
  }
};

}  // namespace

ActsExamples::FatrasSimulation::FatrasSimulation(Config cfg,
                                                 Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("FatrasSimulation", lvl),
      m_cfg(std::move(cfg)) {
  ACTS_DEBUG("hits on sensitive surfaces: " << m_cfg.generateHitsOnSensitive);
  ACTS_DEBUG("hits on material surfaces: " << m_cfg.generateHitsOnMaterial);
  ACTS_DEBUG("hits on passive surfaces: " << m_cfg.generateHitsOnPassive);

  if (!m_cfg.generateHitsOnSensitive && !m_cfg.generateHitsOnMaterial &&
      !m_cfg.generateHitsOnPassive) {
    ACTS_WARNING("FatrasSimulation not configured to generate any hits!");
  }

  if (!m_cfg.trackingGeometry) {
    throw std::invalid_argument{"Missing tracking geometry"};
  }
  if (!m_cfg.magneticField) {
    throw std::invalid_argument{"Missing magnetic field"};
  }
  if (!m_cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers tool");
  }

  // construct the simulation for the specific magnetic field
  m_sim = std::make_unique<FatrasSimulationT>(m_cfg, lvl);
}

// explicit destructor needed for the PIMPL implementation to work
ActsExamples::FatrasSimulation::~FatrasSimulation() {}

ActsExamples::ProcessCode ActsExamples::FatrasSimulation::execute(
    const AlgorithmContext &ctx) const {
  // read input containers
  const auto &inputParticles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);

  ACTS_DEBUG(inputParticles.size() << " input particles");

  // prepare output containers
  SimParticleContainer::sequence_type particlesInitialUnordered;
  SimParticleContainer::sequence_type particlesFinalUnordered;
  SimHitContainer::sequence_type simHitsUnordered;
  // reserve appropriate resources
  particlesInitialUnordered.reserve(inputParticles.size());
  particlesFinalUnordered.reserve(inputParticles.size());
  simHitsUnordered.reserve(inputParticles.size() *
                           m_cfg.averageHitsPerParticle);

  // run the simulation w/ a local random generator
  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);
  auto ret = m_sim->simulate(ctx.geoContext, ctx.magFieldContext, rng,
                             inputParticles, particlesInitialUnordered,
                             particlesFinalUnordered, simHitsUnordered);
  // fatal error leads to panic
  if (not ret.ok()) {
    ACTS_FATAL("event " << ctx.eventNumber << " simulation failed with error "
                        << ret.error());
    return ProcessCode::ABORT;
  }
  // failed particles are just logged. assumes that failed particles are due
  // to edge-cases representing a tiny fraction of the event; not due to a
  // fundamental issue.
  for (const auto &failed : ret.value()) {
    ACTS_ERROR("event " << ctx.eventNumber << " particle " << failed.particle
                        << " failed to simulate with error " << failed.error
                        << ": " << failed.error.message());
  }

  ACTS_DEBUG(particlesInitialUnordered.size()
             << " simulated particles (initial state)");
  ACTS_DEBUG(particlesFinalUnordered.size()
             << " simulated particles (final state)");
  ACTS_DEBUG(simHitsUnordered.size() << " simulated hits");

  // order output containers
  SimParticleContainer particlesInitial;
  SimParticleContainer particlesFinal;
  SimHitContainer simHits;
  particlesInitial.insert(particlesInitialUnordered.begin(),
                          particlesInitialUnordered.end());
  particlesFinal.insert(particlesFinalUnordered.begin(),
                        particlesFinalUnordered.end());
  simHits.insert(simHitsUnordered.begin(), simHitsUnordered.end());
  // store ordered output containers
  ctx.eventStore.add(m_cfg.outputParticlesInitial, std::move(particlesInitial));
  ctx.eventStore.add(m_cfg.outputParticlesFinal, std::move(particlesFinal));
  ctx.eventStore.add(m_cfg.outputSimHits, std::move(simHits));

  return ActsExamples::ProcessCode::SUCCESS;
}
