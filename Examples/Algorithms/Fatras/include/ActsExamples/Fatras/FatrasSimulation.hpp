// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/MagneticField/MagneticField.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"
#include "ActsFatras/Physics/NuclearInteraction/NuclearInteraction.hpp"

#include <memory>
#include <string>

namespace ActsExamples {
namespace detail {
struct FatrasSimulation;
}

/// Fast track simulation using the Acts propagation and navigation.
class FatrasSimulation final : public BareAlgorithm {
 public:
  struct Config {
    /// The particles input collection.
    std::string inputParticles;
    /// The simulated particles initial state collection.
    std::string outputParticlesInitial;
    /// The simulated particles final state collection.
    std::string outputParticlesFinal;
    /// The simulated hits output collection.
    std::string outputSimHits;
    /// Parametrisation of nuclear interaction
    std::string imputParametrisationNuclearInteraction =
        "nuclearInteractionParameters";
    /// Random number service.
    std::shared_ptr<const RandomNumbers> randomNumbers;
    /// The tracking geometry that should be used.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// The magnetic field that should be used.
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;

    // tuning parameters
    /// Minimal absolute momentum for particles to be simulated.
    double pMin = 0.5 * Acts::UnitConstants::GeV;
    /// Simulate (multiple) scattering for charged particles.
    bool emScattering = false;
    /// Simulate ionisiation/excitation energy loss of charged particles.
    bool emEnergyLossIonisation = false;
    /// Simulate radiative energy loss of charged particles.
    bool emEnergyLossRadiation = false;
    /// Simulate electron-positron pair production by photon conversion.
    bool emPhotonConversion = false;
    /// Generate simulation hits on sensitive surfaces.
    bool generateHitsOnSensitive = false;
    /// Generate simulation hits on surfaces with associated material.
    bool generateHitsOnMaterial = false;
    /// Generate simulation hits on passive surfaces, i.e neither sensitive nor
    /// have associated material.
    bool generateHitsOnPassive = false;

    /// Expected average number of hits generated per particle.
    ///
    /// This is just a performance optimization hint and has no impact on the
    /// algorithm function. It is used to guess the amount of memory to
    /// pre-allocate to avoid allocation during event simulation.
    size_t averageHitsPerParticle = 16u;
  };

  /// Add options for the particle selector.
  static void addOptions(Options::Description& desc);
  /// Construct particle selector config from user variables.
  static Config readConfig(const Options::Variables& vars);

  /// Construct the algorithm from a config.
  ///
  /// @param cfg is the configuration struct
  /// @param lvl is the logging level
  FatrasSimulation(Config cfg, Acts::Logging::Level lvl);
  ~FatrasSimulation() final override;

  /// Run the simulation for a single event.
  ///
  /// @param ctx the algorithm context containing all event information
  ActsExamples::ProcessCode execute(
      const AlgorithmContext& ctx) const final override;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::unique_ptr<detail::FatrasSimulation> m_sim;
};

}  // namespace ActsExamples
