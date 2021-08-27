// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/Geant4Options.hpp"

#include "ActsExamples/Utilities/Options.hpp"

#include <string>

#include <boost/program_options.hpp>

void ActsExamples::Options::addGeant4Options(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::bool_switch;
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("g4-rnd-seed1", value<unsigned int>()->default_value(287362910),
      "The first seed of the Geant4 random number generation");
  opt("g4-rnd-seed2", value<unsigned int>()->default_value(730284537),
      "The second seed of the Geant4 random number generation");
  opt("g4-loglevel", value<unsigned int>()->default_value(2),
      "Screen output log level for Geant4 actions.");
}
