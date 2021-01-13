// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Charge.hpp"
#include "Acts/EventData/SingleCurvilinearTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsFatras/EventData/Hit.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Kernel/SimulationResult.hpp"
#include "ActsFatras/Kernel/detail/Interactor.hpp"
#include "ActsFatras/Kernel/detail/SimulatorError.hpp"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <memory>
#include <vector>

namespace ActsFatras {

/// Single particle simulator with a fixed propagator and physics list.
///
/// @tparam generator_t random number generator
/// @tparam continuous_physics_t physics lists for continuous interactions
/// @tparam pointlike_physics_t physics lists for point-like interactions
/// @tparam hit_surface_selector_t sensitive hit surfaces selector
/// @tparam decay_t decay module
template <typename propagator_t, typename continuous_physics_t,
          typename pointlike_physics_t, typename hit_surface_selector_t,
          typename decay_t>
struct ParticleSimulator {
  /// How and within which geometry to propagate the particle.
  propagator_t propagator;
  /// Physics list detailing the simulated continuous interactions.
  /// Will be copied to the per-call interactor.
  continuous_physics_t continuous;
  /// Physics list detailing the simulated point-like interactions.
  /// Will be copied to the per-call interactor.
  pointlike_physics_t pointlike;
  /// Where hits are registiered. Will be copied to the per-call interactor.
  hit_surface_selector_t selectHitSurface;
  /// Decay module.
  decay_t decay;
  /// Local logger for debug output.
  std::shared_ptr<const Acts::Logger> localLogger = nullptr;

  /// Construct the simulator with the underlying propagator.
  ParticleSimulator(propagator_t &&propagator_, Acts::Logging::Level lvl)
      : propagator(propagator_),
        localLogger(Acts::getDefaultLogger("Simulator", lvl)) {}

  /// Provide access to the local logger instance, e.g. for logging macros.
  const Acts::Logger &logger() const { return *localLogger; }

  /// Simulate a single particle without secondaries.
  ///
  /// @param geoCtx is the geometry context to access surface geometries
  /// @param magCtx is the magnetic field context to access field values
  /// @param generator is the random number generator
  /// @param particle is the initial particle state
  /// @returns the result of the corresponding Interactor propagator action.
  ///
  /// @tparam generator_t is the type of the random number generator
  template <typename generator_t>
  Acts::Result<SimulationResult> simulate(
      const Acts::GeometryContext &geoCtx,
      const Acts::MagneticFieldContext &magCtx, generator_t &generator,
      const Particle &particle) const {
    assert(localLogger and "Missing local logger");

    // propagator-related additional types
    using Interactor =
        detail::Interactor<generator_t, continuous_physics_t,
                           pointlike_physics_t, hit_surface_selector_t>;
    using InteractorResult = typename Interactor::result_type;
    using Actions = Acts::ActionList<Interactor>;
    using Abort = Acts::AbortList<typename Interactor::ParticleNotAlive,
                                  Acts::EndOfWorldReached>;
    using PropagatorOptions = Acts::PropagatorOptions<Actions, Abort>;

    // Construct per-call options.
    PropagatorOptions options(geoCtx, magCtx,
                              Acts::LoggerWrapper{*localLogger});
    options.absPdgCode = Acts::makeAbsolutePdgParticle(particle.pdg());
    options.mass = particle.mass();
    // setup the interactor as part of the propagator options
    auto &interactor = options.actionList.template get<Interactor>();
    interactor.generator = &generator;
    interactor.continuous = continuous;
    interactor.pointlike = pointlike;
    interactor.selectHitSurface = selectHitSurface;
    interactor.initialParticle = particle;
    interactor.properTimeLimit =
        decay.generateProperTimeLimit(generator, particle);
    // use AnyCharge to be able to handle neutral and charged parameters
    Acts::SingleCurvilinearTrackParameters<Acts::AnyCharge> start(
        particle.fourPosition(), particle.unitDirection(),
        particle.absoluteMomentum(), particle.charge());
    auto result = propagator.propagate(start, options);
    if (not result.ok()) {
      return result.error();
    }
    auto &value = result.value().template get<InteractorResult>();

    // add decay products to the generated particles if the particle has decayed
    if (value.particleStatus == SimulationParticleStatus::eDecayed) {
      auto descendants = decay.run(generator, value.particle);
      for (auto &&descendant : descendants) {
        value.generatedParticles.emplace_back(std::move(descendant));
      }
    }

    return std::move(value);
  }
};

/// Multi-particle simulator.
///
/// @tparam charged_selector_t Callable selector type for charged particles
/// @tparam charged_simulator_t Single particle simulator for charged particles
/// @tparam neutral_selector_t Callable selector type for neutral particles
/// @tparam neutral_simulator_t Single particle simulator for neutral particles
template <typename charged_selector_t, typename charged_simulator_t,
          typename neutral_selector_t, typename neutral_simulator_t>
struct Simulator {
  /// A particle that failed to simulate.
  struct FailedParticle {
    /// Initial particle state of the failed particle.
    ///
    /// This must store the full particle state to be able to handle secondaries
    /// that are not in the input particle list. Otherwise they could not be
    /// referenced.
    Particle particle;
    /// The associated error code for this particular failure case.
    std::error_code error;
  };

  charged_selector_t selectCharged;
  neutral_selector_t selectNeutral;
  charged_simulator_t charged;
  neutral_simulator_t neutral;

  /// Construct from the single charged/neutral particle simulators.
  Simulator(charged_simulator_t &&charged_, neutral_simulator_t &&neutral_)
      : charged(std::move(charged_)), neutral(std::move(neutral_)) {}

  /// Simulate multiple particles and generated secondaries.
  ///
  /// @param geoCtx is the geometry context to access surface geometries
  /// @param magCtx is the magnetic field context to access field values
  /// @param inputParticles contains all particles that should be simulated
  /// @param simulatedParticlesInitial contains initial particle states
  /// @param simulatedParticlesFinal contains final particle states
  /// @param hits contains all generated hits
  /// @retval Acts::Result::Error if there is a fundamental issue
  /// @retval Acts::Result::Success with all particles that failed to simulate
  ///
  /// @warning Particle-hit association is based on particle ids generated
  ///          during the simulation. This requires that all input particles
  ///          **must** have generation and sub-particle number set to zero.
  /// @note Parameter edge-cases can lead to errors in the underlying propagator
  ///       and thus to particles that fail to simulate. Here, full events are
  ///       simulated and the failure to simulate one particle should not be
  ///       considered a general failure of the simulator. Instead, a list of
  ///       particles that fail to simulate is provided to the user. It is the
  ///       users responsibility to handle them.
  /// @note Failed particles are removed from the regular output, i.e. they do
  ///       not appear in the simulated particles containers nor do they
  ///       generate hits.
  ///
  /// This takes all input particles and simulates those passing the selection
  /// using the appropriate simulator. All selected particle states including
  /// additional ones generated from interactions are stored in separate output
  /// containers; both the initial state at the production vertex and the final
  /// state after propagation are stored. Hits generated from selected input and
  /// generated particles are stored in the hit container.
  ///
  /// @tparam generator_t is the type of the random number generator
  /// @tparam input_particles_t is a Container for particles
  /// @tparam output_particles_t is a SequenceContainer for particles
  /// @tparam hits_t is a SequenceContainer for hits
  template <typename generator_t, typename input_particles_t,
            typename output_particles_t, typename hits_t>
  Acts::Result<std::vector<FailedParticle>> simulate(
      const Acts::GeometryContext &geoCtx,
      const Acts::MagneticFieldContext &magCtx, generator_t &generator,
      const input_particles_t &inputParticles,
      output_particles_t &simulatedParticlesInitial,
      output_particles_t &simulatedParticlesFinal, hits_t &hits) const {
    assert(
        (simulatedParticlesInitial.size() == simulatedParticlesFinal.size()) and
        "Inconsistent initial sizes of the simulated particle containers");

    using ParticleSimulatorResult = Acts::Result<SimulationResult>;

    std::vector<FailedParticle> failedParticles;

    for (const Particle &inputParticle : inputParticles) {
      // only consider simulatable particles
      if (not selectParticle(inputParticle)) {
        continue;
      }
      // required to allow correct particle id numbering for secondaries later
      if ((inputParticle.particleId().generation() != 0u) or
          (inputParticle.particleId().subParticle() != 0u)) {
        return detail::SimulatorError::eInvalidInputParticleId;
      }

      // Do a *depth-first* simulation of the particle and its secondaries,
      // i.e. we simulate all secondaries, tertiaries, ... before simulating
      // the next primary particle. Use the end of the output container as
      // a queue to store particles that should be simulated.
      //
      // WARNING the initial particle state output container will be modified
      //         during iteration. New secondaries are added to and failed
      //         particles might be removed. To avoid issues, access must always
      //         occur via indices.
      auto iinitial = simulatedParticlesInitial.size();
      simulatedParticlesInitial.push_back(inputParticle);
      for (; iinitial < simulatedParticlesInitial.size(); ++iinitial) {
        const auto &initialParticle = simulatedParticlesInitial[iinitial];

        // only simulatable particles are pushed to the container.
        // they must therefore be either charged or neutral.
        ParticleSimulatorResult result = ParticleSimulatorResult::success({});
        if (selectCharged(initialParticle)) {
          result = charged.simulate(geoCtx, magCtx, generator, initialParticle);
        } else {
          result = neutral.simulate(geoCtx, magCtx, generator, initialParticle);
        }

        if (not result.ok()) {
          // remove particle from output container since it was not simulated.
          simulatedParticlesInitial.erase(
              std::next(simulatedParticlesInitial.begin(), iinitial));
          // record the particle as failed
          failedParticles.push_back({initialParticle, result.error()});
          continue;
        }

        copyOutputs(result.value(), simulatedParticlesInitial,
                    simulatedParticlesFinal, hits);
        // since physics processes are independent, there can be particle id
        // collisions within the generated secondaries. they can be resolved by
        // renumbering within each sub-particle generation. this must happen
        // before the particle is simulated since the particle id is used to
        // associate generated hits back to the particle.
        renumberTailParticleIds(simulatedParticlesInitial, iinitial);
      }
    }

    // the overall function call succeeded, i.e. no fatal errors occured.
    // yet, there might have been some particle for which the propagation
    // failed. thus, the successful result contains a list of failed particles.
    // sounds a bit weird, but that is the way it is.
    return failedParticles;
  }

 private:
  /// Select if the particle should be simulated at all.
  ///
  /// This also enforces mutual-exclusivity of the two charge selections. If
  /// both charge selections evaluate true, they are probably not setup
  /// correctly and not simulating them at all is a reasonable fall-back.
  bool selectParticle(const Particle &particle) const {
    const bool isValidCharged = selectCharged(particle);
    const bool isValidNeutral = selectNeutral(particle);
    assert(not(isValidCharged and isValidNeutral) and
           "Charge selectors are not mutually exclusive");
    return isValidCharged xor isValidNeutral;
  }

  /// Copy Interactor results to output containers.
  ///
  /// @tparam particles_t is a SequenceContainer for particles
  /// @tparam hits_t is a SequenceContainer for hits
  template <typename particles_t, typename hits_t>
  void copyOutputs(const SimulationResult &result,
                   particles_t &particlesInitial, particles_t &particlesFinal,
                   hits_t &hits) const {
    // initial particle state was already pushed to the container before
    // store final particle state at the end of the simulation
    particlesFinal.push_back(result.particle);
    // move generated secondaries that should be simulated to the output
    std::copy_if(
        result.generatedParticles.begin(), result.generatedParticles.end(),
        std::back_inserter(particlesInitial),
        [this](const Particle &particle) { return selectParticle(particle); });
    std::copy(result.hits.begin(), result.hits.end(), std::back_inserter(hits));
  }

  /// Renumber particle ids in the tail of the container.
  ///
  /// Ensures particle ids are unique by modifying the sub-particle number
  /// within each generation.
  ///
  /// @param particles particle container in which particles are renumbered
  /// @param lastValid index of the last particle with a valid particle id
  ///
  /// @tparam particles_t is a SequenceContainer for particles
  ///
  /// @note This function assumes that all particles in the tail have the same
  ///       vertex numbers (primary/secondary) and particle number and are
  ///       ordered according to their generation number.
  ///
  template <typename particles_t>
  static void renumberTailParticleIds(particles_t &particles,
                                      std::size_t lastValid) {
    // iterate over adjacent pairs; potentially modify the second element.
    // assume e.g. a primary particle 2 with generation=subparticle=0 that
    // generates two secondaries during simulation. we have the following
    // ids (decoded as vertex|particle|generation|subparticle)
    //
    //     v|2|0|0, v|2|1|0, v|2|1|0
    //
    // where v represents the vertex numbers. this will be renumbered to
    //
    //     v|2|0|0, v|2|1|0, v|2|1|1
    //
    // if each secondary then generates two tertiaries we could have e.g.
    //
    //     v|2|0|0, v|2|1|0, v|2|1|1, v|2|2|0, v|2|2|1, v|2|2|0, v|2|2|1
    //
    // which would then be renumbered to
    //
    //     v|2|0|0, v|2|1|0, v|2|1|1, v|2|2|0, v|2|2|1, v|2|2|2, v|2|2|3
    //
    for (auto j = lastValid; (j + 1u) < particles.size(); ++j) {
      const auto prevId = particles[j].particleId();
      auto currId = particles[j + 1u].particleId();
      // NOTE primary/secondary vertex and particle are assumed to be equal
      // only renumber within one generation
      if (prevId.generation() != currId.generation()) {
        continue;
      }
      // ensure the sub-particle is strictly monotonic within a generation
      if (prevId.subParticle() < currId.subParticle()) {
        continue;
      }
      // sub-particle numbering must be non-zero
      currId.setSubParticle(prevId.subParticle() + 1u);
      particles[j + 1u] = particles[j + 1u].withParticleId(currId);
    }
  }
};

}  // namespace ActsFatras
