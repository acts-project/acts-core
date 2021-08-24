// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/EventStoreRegistry.hpp"

std::vector<ActsExamples::WhiteBoard*>
    ActsExamples::EventStoreRegistry::boards =
        std::vector<ActsExamples::WhiteBoard*>();

std::vector<ActsExamples::SimHitContainer::sequence_type>
    ActsExamples::EventStoreRegistry::hits =
        std::vector<ActsExamples::SimHitContainer::sequence_type>();

std::vector<ActsExamples::SimParticleContainer::sequence_type>
    ActsExamples::EventStoreRegistry::particlesInitial =
        std::vector<ActsExamples::SimParticleContainer::sequence_type>();

std::vector<ActsExamples::SimParticleContainer::sequence_type>
    ActsExamples::EventStoreRegistry::particlesFinal =
        std::vector<ActsExamples::SimParticleContainer::sequence_type>();

ActsExamples::EventStoreRegistry::EventStoreRegistry(size_t nevents) {
  boards = std::vector<WhiteBoard*>(nevents, nullptr);
  hits = std::vector<SimHitContainer::sequence_type>(
      nevents, SimHitContainer::sequence_type());
  particlesInitial =
      std::vector<ActsExamples::SimParticleContainer::sequence_type>(
          nevents, ActsExamples::SimParticleContainer::sequence_type());
  particlesFinal =
      std::vector<ActsExamples::SimParticleContainer::sequence_type>(
          nevents, ActsExamples::SimParticleContainer::sequence_type());
}
