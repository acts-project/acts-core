// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Utilities/OptionsFwd.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"

namespace ActsExamples {
namespace Simulation {

    /// Collection names
    static constexpr const char* kParticlesInput = "particles_input";
    static constexpr const char* kParticlesSelection = "particles_selection";
    static constexpr const char* kParticlesInitial = "particles_initial";
    static constexpr const char* kParticlesFinal = "particles_final";
    static constexpr const char* kSimHits = "simhits";

    /// Add input options 
    ///
    /// @param desc is the boost options descr format
    void addInputOptions(ActsExamples::Options::Description& desc);

    /// Setup the input from the provided options
    ///
    /// @param vars the parsed variables from the boost options
    /// @param sequencer the non-const sequencer
    /// @param randomNumbers the randomNumbers shared pointer
    /// 
    void setupInput(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers);

    /// Setup the output from the provided options
    ///
    /// @param vars the parsed variables from the boost options
    /// @param sequencer the non-const sequencer
    void setupOutput(const ActsExamples::Options::Variables& vars,
                 ActsExamples::Sequencer& sequencer);


} // namespace Simulation
} // namespace ActsExamples
