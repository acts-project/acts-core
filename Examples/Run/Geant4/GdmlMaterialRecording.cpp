// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"
#include "ActsExamples/Geant4/Geant4Options.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/ParticleGunOptions.hpp"

#include <boost/program_options.hpp>

#include "Geant4.hpp"

using namespace ActsExamples;

int main(int argc, char* argv[]) {
  // Setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::Root);
  Options::addGeant4Options(desc);
  Options::addParticleGunOptions(desc);
  Options::addRandomNumbersOptions(desc);
  desc.add_options()(
      "gdml-file",
      boost::program_options::value<std::string>()->default_value(""),
      "GDML detector file.");

  auto vm = Options::parse(desc, argc, argv);

  if (vm.empty()) {
    return EXIT_FAILURE;
  }
  auto gdmlFile = vm["gdml-file"].as<std::string>();

  return runMaterialRecording(
      vm, std::make_unique<GdmlDetectorConstructionFactory>(gdmlFile));
}
