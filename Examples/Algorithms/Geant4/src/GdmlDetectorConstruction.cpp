// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"

#include <G4GDMLParser.hh>

using namespace ActsExamples;

GdmlDetectorConstruction::GdmlDetectorConstruction(std::string path)
    : G4VUserDetectorConstruction(), m_path(std::move(path)) {}

G4VPhysicalVolume* GdmlDetectorConstruction::Construct() {
  if (m_world == nullptr) {
    G4GDMLParser parser;
    // TODO how to handle errors
    parser.Read(m_path);
    m_world = parser.GetWorldVolume();
  }
  return m_world;
}

GdmlDetectorConstructionFactory::GdmlDetectorConstructionFactory(
    const std::string& path)
    : m_path{path} {}

std::unique_ptr<G4VUserDetectorConstruction>
GdmlDetectorConstructionFactory::operator()() const {
  return std::make_unique<GdmlDetectorConstruction>(m_path);
}
