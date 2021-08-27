// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>
#include <vector>

#include <G4VUserPhysicsList.hh>

namespace ActsExamples {

/// @class MaterialPhysicsList
///
/// @brief Stripped down physics list to Geantinos and only Transport
/// This speeds up the initialization of the MaterialRecording job
///
class MaterialPhysicsList final : public G4VUserPhysicsList {
 public:
  /// Construct the action
  ///
  /// @param cfg the configuration struct for this Stepping action
  /// @param logger is an Acts::Logger for unique logging
  MaterialPhysicsList(std::unique_ptr<const Acts::Logger> logger =
                          Acts::getDefaultLogger("MaterialPhysicsList",
                                                 Acts::Logging::INFO));
  ~MaterialPhysicsList() final override = default;

  /// @brief Interface particle construction method
  ///
  /// This will add only Geantino and ChargedGeantino
  void ConstructParticle() final override;

  /// @brief Interface process construction method
  ///
  /// This will add only Transport
  void ConstructProcess() final override;

  /// @brief Interface process construction method
  ///
  /// This will add only Transport
  void SetCuts() final override;

 private:
  /// Private access method to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// The looging instance
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace ActsExamples
