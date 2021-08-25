// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

#include <memory>
#include <vector>

class G4VUserPrimaryGeneratorAction;
class G4VUserDetectorConstruction;
class G4UserRunAction;
class G4UserEventAction;
class G4UserTrackingAction;
class G4UserSteppingAction;

namespace ActsExamples {

class G4DetectorConstructionFactory;

/// Main function for running Geant4 with a specific detector.
///
/// @param vars the parsed variables
/// @param sequencer the event sequencer
/// @param detector is the detector to be used
/// @param generatorAction the Geant4 user generator action
/// @param runActions the list of Geant4 user run action
/// @param eventActions the list of Geant4 user event action
/// @param trackingActions the list of Geant4 user tracking action
/// @param steppingActions the list of Geant4 user stepping action
///
void setupGeant4Simulation(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer, G4VUserDetectorConstruction* detector,
    G4VUserPrimaryGeneratorAction* generatorAction,
    std::vector<G4UserRunAction*> runActions = {},
    std::vector<G4UserEventAction*> eventActions = {},
    std::vector<G4UserTrackingAction*> trackingActions = {},
    std::vector<G4UserSteppingAction*> steppingActions = {});

/// Specific setup: Material Recording
///
/// @param vars the parsed variables
/// @param sequencer the event sequencer
/// @param g4DetectorFactory is the detector to be used
int runGeantinoRecording(
    const ActsExamples::Options::Variables& vars,
    std::shared_ptr<ActsExamples::G4DetectorConstructionFactory>
        g4DetectorFactory);

}  // namespace ActsExamples
