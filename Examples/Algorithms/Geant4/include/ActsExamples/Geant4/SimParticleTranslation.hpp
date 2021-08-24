// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>

#include <G4VUserPrimaryGeneratorAction.hh>

class G4Event;

namespace ActsExamples {

/// @class SimParticleTranslation
///
/// Translates the Generated SimParticles from ACTS into
/// Geant4 particles and acts as the associated Geant4 primary
/// generator action.
///
/// This also ensures that the event numbers correspond to the
/// ACTS framework event numbers and hence harmonizes the EventStoreRegistry.
class SimParticleTranslation final : public G4VUserPrimaryGeneratorAction {
 public:
  /// Nested configuration struct that contains the
  /// input particle collection name,
  struct Config {
    /// The input particle collection
    std::string inputParticles = "";

    /// The number of hits per particle to be expected
    /// @note best to include secondaries for that
    unsigned int reserveHitsPerParticle = 20;
  };

  /// Construct the generator action
  ///
  /// @param cfg the configuration struct
  /// @param logger the ACTS logging instance
  SimParticleTranslation(const Config& cfg,
                         std::unique_ptr<const Acts::Logger> logger =
                             Acts::getDefaultLogger("SimParticleTranslation",
                                                    Acts::Logging::INFO));

  ~SimParticleTranslation() final override;

  /// Interface method to generate the primary
  ///
  /// @param anEvent is the event that will be run
  void GeneratePrimaries(G4Event* anEvent) final override;

 protected:
  Config m_cfg;

 private:
  /// Event number cache for EventStoreRegistry harmonization
  unsigned int m_eventNr = 0;

  /// Private access method to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// The looging instance
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace ActsExamples
