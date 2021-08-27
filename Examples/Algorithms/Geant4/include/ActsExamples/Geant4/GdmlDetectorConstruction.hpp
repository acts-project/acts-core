// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Geant4/G4DetectorConstructionFactory.hpp"

#include <G4VUserDetectorConstruction.hh>

namespace ActsExamples {

/// Construct the Geant4 detector from a Gdml file.
class GdmlDetectorConstruction final : public G4VUserDetectorConstruction {
 public:
  /// @param path is the path to the Gdml file
  GdmlDetectorConstruction(std::string path);

  /// Read the file and parse it to construct the Geant4 description
  ///
  /// @note to simplify further setup withiin the ACTS framework
  /// the world is cached when first constructed
  G4VPhysicalVolume* Construct() final override;

 private:
  /// Path to the Gdml file
  std::string m_path;
  /// Cached worled volume
  G4VPhysicalVolume* m_world;
};

class GdmlDetectorConstructionFactory : public G4DetectorConstructionFactory {
 public:
  /// @brief Construct a new GDML Detector Factory
  ///
  /// @param path The input GDML file path
  GdmlDetectorConstructionFactory(const std::string& path);

  /// @brief Main factory method following the interface
  ///
  /// @return Detector construction based on GDML.
  std::unique_ptr<G4VUserDetectorConstruction> operator()() const override;

 private:
  std::string m_path;
};

}  // namespace ActsExamples
