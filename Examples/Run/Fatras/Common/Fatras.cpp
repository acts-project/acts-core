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
#include "ActsExamples/TruthTracking/ParticleSelector.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <memory>
#include <string>

#include <boost/program_options.hpp>

using namespace ActsExamples;

namespace {

// collection names
static constexpr const char* kParticlesInput = "particles_input";
static constexpr const char* kParticlesSelection = "particles_selection";
static constexpr const char* kParticlesInitial = "particles_initial";
static constexpr const char* kParticlesFinal = "particles_final";
static constexpr const char* kSimHits = "simhits";

// input handling

void addInputOptions(ActsExamples::Options::Description& desc) {
  using boost::program_options::value;

  auto opt = desc.add_options();
  // do not use the common input options: they have the wrong doc text and add a
  // lot of unused options.
  opt("input-dir", value<std::string>(),
      "Read particle input from CSV files in the given directory. If not "
      "given, particles are generated with the particle gun");
  ActsExamples::Options::addParticleGunOptions(desc);
  ActsExamples::ParticleSelector::addOptions(desc);
}

void setupInput(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers) {
  auto logLevel = Options::readLogLevel(vars);

  if (not vars["input-dir"].empty()) {
    // read particle input from csv file
    CsvParticleReader::Config readParticles;
    readParticles.inputDir = vars["input-dir"].as<std::string>();
    readParticles.inputStem = "particles";
    readParticles.outputParticles = kParticlesInput;
    sequencer.addReader(
        std::make_shared<CsvParticleReader>(readParticles, logLevel));

  } else {
    // generate particle input from a particle gun
    auto gen = Options::readParticleGunOptions(vars);
    gen.outputParticles = kParticlesInput;
    gen.randomNumbers = randomNumbers;
    sequencer.addReader(std::make_shared<EventGenerator>(gen, logLevel));
  }

  // add additional particle selection
  auto select = ActsExamples::ParticleSelector::readConfig(vars);
  select.inputParticles = kParticlesInput;
  select.outputParticles = kParticlesSelection;
  sequencer.addAlgorithm(
      std::make_shared<ActsExamples::ParticleSelector>(select, logLevel));
}

// output handling

// output options are just the common output options

void setupOutput(const ActsExamples::Options::Variables& vars,
                 ActsExamples::Sequencer& sequencer) {
  auto logLevel = Options::readLogLevel(vars);
  auto outputDir =
      ensureWritableDirectory(vars["output-dir"].template as<std::string>());

  // Write simulation information as CSV files
  if (vars["output-csv"].template as<bool>()) {
    // NOTE: the collection names are an internal implementation detail but the
    //   output file stems are a public user interface. the former can be
    //   changed at-will, but the later should always be the same.

    // write simulated particle initial states
    CsvParticleWriter::Config writeInitial;
    writeInitial.inputParticles = kParticlesInitial;
    writeInitial.outputDir = outputDir;
    writeInitial.outputStem = "particles_initial";
    sequencer.addWriter(
        std::make_shared<CsvParticleWriter>(writeInitial, logLevel));

    // write simulated particle final states
    CsvParticleWriter::Config writeFinal;
    writeFinal.inputParticles = kParticlesFinal;
    writeFinal.outputDir = outputDir;
    writeFinal.outputStem = "particles_final";
    sequencer.addWriter(
        std::make_shared<CsvParticleWriter>(writeFinal, logLevel));

    // write simulated hits
    CsvSimHitWriter::Config writeSimHits;
    writeSimHits.inputSimHits = kSimHits;
    writeSimHits.outputDir = outputDir;
    writeSimHits.outputStem = "hits";
    sequencer.addWriter(
        std::make_shared<CsvSimHitWriter>(writeSimHits, logLevel));
  }

  // Write simulation information as ROOT files
  if (vars["output-root"].template as<bool>()) {
    // NOTE: the collection names are an internal implementation detail but the
    //   output file names are a public user interface. the former can be
    //   changed at-will, but the later should always be the same.

    // write simulated particle initial states
    RootParticleWriter::Config writeInitial;
    writeInitial.inputParticles = kParticlesInitial;
    writeInitial.filePath = joinPaths(outputDir, "particles_initial.root");
    sequencer.addWriter(
        std::make_shared<RootParticleWriter>(writeInitial, logLevel));

    // write simulated particle final states
    RootParticleWriter::Config writeFinal;
    writeFinal.inputParticles = kParticlesFinal;
    writeFinal.filePath = joinPaths(outputDir, "particles_final.root");
    sequencer.addWriter(
        std::make_shared<RootParticleWriter>(writeFinal, logLevel));

    // write simulated hits
    RootSimHitWriter::Config writeHits;
    writeHits.inputSimHits = kSimHits;
    writeHits.filePath = joinPaths(outputDir, "hits.root");
    sequencer.addWriter(
        std::make_shared<RootSimHitWriter>(writeHits, logLevel));
  }
}

// simulation handling

void setupSimulation(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry) {
  auto logLevel = Options::readLogLevel(vars);
  auto fatrasCfg = FatrasSimulation::readConfig(vars);
  fatrasCfg.inputParticles = kParticlesSelection;
  fatrasCfg.outputParticlesInitial = kParticlesInitial;
  fatrasCfg.outputParticlesFinal = kParticlesFinal;
  fatrasCfg.outputSimHits = kSimHits;
  fatrasCfg.randomNumbers = randomNumbers;
  fatrasCfg.trackingGeometry = trackingGeometry;
  fatrasCfg.magneticField = ActsExamples::Options::readMagneticField(vars);

  sequencer.addAlgorithm(
      std::make_shared<FatrasSimulation>(std::move(fatrasCfg), logLevel));
}

}  // namespace

// fatras main function

int runFatras(int argc, char* argv[],
              std::shared_ptr<ActsExamples::IBaseDetector> detector) {
  using namespace ActsExamples;

  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  addInputOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::Root | OutputFormat::Csv);
  // add general and detector-specific geometry options
  Options::addGeometryOptions(desc);
  detector->addOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addMagneticFieldOptions(desc);
  // algorithm-specific options
  FatrasSimulation::addOptions(desc);

  auto vars = Options::parse(desc, argc, argv);
  if (vars.empty()) {
    return EXIT_FAILURE;
  }

  // basic services
  auto randomNumbers =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vars));

  // setup sequencer
  Sequencer sequencer(Options::readSequencerConfig(vars));
  // setup detector geometry and material and the magnetic field
  auto [trackingGeometry, contextDecorators] = Geometry::build(vars, *detector);
  for (auto cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }
  // setup algorithm chain
  setupInput(vars, sequencer, randomNumbers);
  setupSimulation(vars, sequencer, randomNumbers, trackingGeometry);
  setupOutput(vars, sequencer);

  // run the simulation
  return sequencer.run();
}
