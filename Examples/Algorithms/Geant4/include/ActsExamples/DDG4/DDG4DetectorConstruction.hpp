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

namespace dd4hep {
class Detector;
}
namespace ActsExamples {

/// Construct the Geant4 detector from a DD4hep description.
class DDG4DetectorConstruction final : public G4VUserDetectorConstruction {
 public:
  DDG4DetectorConstruction(dd4hep::Detector& detector);

  /// Convert the stored DD4hep detector to a Geant4 description.
  ///
  /// Transfers ownership of the created object as all volumes (including world)
  /// are deleted in ~G4PhysicalVolumeStore().
  ///
  /// @note for facilitating configuration within the ACTS framework the world
  /// volume is cached
  G4VPhysicalVolume* Construct() final override;

 private:
  /// The DD4hep detector instrance
  dd4hep::Detector& m_detector;
  /// The world volume
  G4VPhysicalVolume* m_world = nullptr;
};

class DDG4DetectorConstructionFactory : public G4DetectorConstructionFactory {
 public:
  /// @brief Construct a new DD4hep-Geant4 detector factory
  ///
  /// @param detector DD4hep detector instance to construct G4 geometry from
  DDG4DetectorConstructionFactory(dd4hep::Detector& detector);

  /// @brief Main factory method
  ///
  /// @return Detector construction based on a DD4hep geometry
  std::unique_ptr<G4VUserDetectorConstruction> operator()() const override;

 private:
  dd4hep::Detector& m_detector;
};

}  // namespace ActsExamples
