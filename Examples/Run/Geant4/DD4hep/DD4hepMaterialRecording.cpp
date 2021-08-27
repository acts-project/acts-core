// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/DD4hepDetectorOptions.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"
#include "ActsExamples/DDG4/DDG4DetectorConstruction.hpp"
#include "ActsExamples/Geant4/Geant4Options.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/ParticleGunOptions.hpp"

#include <boost/program_options.hpp>

#include "../Common/Geant4.hpp"

using namespace ActsExamples;

int main(int argc, char* argv[]) {
  // Setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::Root);
  Options::addDD4hepOptions(desc);
  Options::addGeant4Options(desc);
  Options::addParticleGunOptions(desc);
  Options::addRandomNumbersOptions(desc);
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  // Setup the DD4hep detector
  auto dd4hepCfg = Options::readDD4hepConfig<po::variables_map>(vm);
  auto geometrySvc = std::make_shared<DD4hep::DD4hepGeometryService>(dd4hepCfg);

  return runMaterialRecording(
      vm,
      std::make_unique<DDG4DetectorConstructionFactory>(*geometrySvc->lcdd()));
}
