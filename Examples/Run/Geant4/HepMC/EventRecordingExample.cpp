// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepDetectorOptions.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"
#include "ActsExamples/DDG4/DDG4DetectorConstruction.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Geant4/Geant4Options.hpp"
#include "ActsExamples/Geant4HepMC/EventRecording.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Options.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Writer.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <fstream>
#include <string>

#include <boost/program_options.hpp>

int main(int argc, char* argv[]) {
  // Declare the supported program options.
  // Setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  ActsExamples::Options::addSequencerOptions(desc);
  ActsExamples::Options::addOutputOptions(
      desc, ActsExamples::OutputFormat::DirectoryOnly);
  ActsExamples::Options::addInputOptions(desc);
  ActsExamples::Options::addDD4hepOptions(desc);
  ActsExamples::Options::addGeant4Options(desc);
  ActsExamples::Options::addHepMC3WriterOptions(desc);

  auto vm = ActsExamples::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  auto logLevel = ActsExamples::Options::readLogLevel(vm);

  ActsExamples::Sequencer sequencer(
      ActsExamples::Options::readSequencerConfig(vm));

  // Read particles (initial states) and clusters from CSV files
  auto particleReader = ActsExamples::Options::readCsvParticleReaderConfig(vm);
  particleReader.inputStem = "particles";
  particleReader.outputParticles = "particles";
  sequencer.addReader(std::make_shared<ActsExamples::CsvParticleReader>(
      particleReader, logLevel));

  // Prepare the detector
  auto dd4hepCfg = ActsExamples::Options::readDD4hepConfig(vm);
  auto geometrySvc =
      std::make_shared<ActsExamples::DD4hep::DD4hepGeometryService>(dd4hepCfg);

  // Prepare the recording
  ActsExamples::EventRecording::Config erConfig;
  erConfig.inputParticles = particleReader.outputParticles;
  erConfig.outputHepMcTracks = "geant-event";
  erConfig.detectorConstructionFactory =
      std::make_unique<ActsExamples::DDG4DetectorConstructionFactory>(
          *geometrySvc->lcdd());
  erConfig.seed1 = vm["g4-rnd-seed1"].as<unsigned int>();
  erConfig.seed2 = vm["g4-rnd-seed2"].as<unsigned int>();

  // Create the writer
  auto hepMC3WriterConfig = ActsExamples::Options::readHepMC3WriterOptions(vm);
  hepMC3WriterConfig.inputEvents = erConfig.outputHepMcTracks;

  // Add to the sequencer
  sequencer.addAlgorithm(std::make_shared<ActsExamples::EventRecording>(
      std::move(erConfig), logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::HepMC3AsciiWriter>(
      std::move(hepMC3WriterConfig), logLevel));

  // Run
  return sequencer.run();
}
