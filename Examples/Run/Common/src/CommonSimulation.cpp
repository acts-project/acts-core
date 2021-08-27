// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Simulation/CommonSimulation.hpp"

#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleWriter.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitWriter.hpp"
#include "ActsExamples/Io/Root/RootParticleWriter.hpp"
#include "ActsExamples/Io/Root/RootSimHitWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/ParticleGunOptions.hpp"
#include "ActsExamples/TruthTracking/ParticleSelector.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <boost/program_options.hpp>

namespace ActsExamples {
namespace Simulation {

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

}  // namespace Simulation
}  // namespace ActsExamples
