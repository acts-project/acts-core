// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Fatras.hpp"

#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Fatras/FatrasSimulation.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Generators/EventGenerator.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleWriter.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitWriter.hpp"
#include "ActsExamples/Io/Root/RootParticleWriter.hpp"
#include "ActsExamples/Io/Root/RootSimHitWriter.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/ParticleGunOptions.hpp"
#include "ActsExamples/Simulation/CommonSimulation.hpp"
#include "ActsExamples/TruthTracking/ParticleSelector.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <memory>
#include <string>

#include <boost/program_options.hpp>

using namespace ActsExamples;

namespace {

// simulation handling

void setupFatrasSimulation(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry) {
  auto logLevel = Options::readLogLevel(vars);
  auto fatrasCfg = FatrasSimulation::readConfig(vars);
  fatrasCfg.inputParticles = Simulation::kParticlesSelection;
  fatrasCfg.outputParticlesInitial = Simulation::kParticlesInitial;
  fatrasCfg.outputParticlesFinal = Simulation::kParticlesFinal;
  fatrasCfg.outputSimHits = Simulation::kSimHits;
  fatrasCfg.randomNumbers = randomNumbers;
  fatrasCfg.trackingGeometry = trackingGeometry;
  fatrasCfg.magneticField = ActsExamples::Options::readMagneticField(vars);

  sequencer.addAlgorithm(
      std::make_shared<FatrasSimulation>(std::move(fatrasCfg), logLevel));
}

}  // namespace

/// Fatras main function
///
/// Standard arguments @param argc and @param argv[] are forwarded
/// @param detector abstracts the used detector input
int runFatras(int argc, char* argv[],
              std::shared_ptr<ActsExamples::IBaseDetector> detector) {
  using namespace ActsExamples;

  // Setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Simulation::addInputOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::Root | OutputFormat::Csv);
  // Add general and detector-specific geometry options
  Options::addGeometryOptions(desc);
  detector->addOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addMagneticFieldOptions(desc);
  // Algorithm-specific options
  FatrasSimulation::addOptions(desc);

  auto vars = Options::parse(desc, argc, argv);
  if (vars.empty()) {
    return EXIT_FAILURE;
  }

  // Basic services
  auto randomNumbers =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vars));

  // Setup sequencer
  Sequencer sequencer(Options::readSequencerConfig(vars));
  // Setup detector geometry and material and the magnetic field
  auto [trackingGeometry, contextDecorators] = Geometry::build(vars, *detector);
  for (auto cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }
  // Setup algorithm chain
  Simulation::setupInput(vars, sequencer, randomNumbers);
  setupFatrasSimulation(vars, sequencer, randomNumbers, trackingGeometry);
  Simulation::setupOutput(vars, sequencer);

  // Run the simulation
  return sequencer.run();
}
