// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/SurfaceMaterialMapper.hpp"
#include "Acts/Plugins/Json/MaterialMapJsonConverter.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Json/JsonMaterialWriter.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackReader.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackWriter.hpp"
#include "ActsExamples/Io/Root/RootMaterialWriter.hpp"
#include "ActsExamples/MaterialMapping/MaterialMapping.hpp"
#include "ActsExamples/MaterialMapping/MaterialMappingOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Propagation/PropagationOptions.hpp"
#include "ActsExamples/Simulation/CommonSimulation.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <memory>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

int runMaterialMapping(int argc, char* argv[],
                       ActsExamples::IBaseDetector& detector) {
  // Setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  ActsExamples::Options::addSequencerOptions(desc);
  ActsExamples::Options::addGeometryOptions(desc);
  ActsExamples::Options::addMaterialOptions(desc);
  ActsExamples::Options::addMaterialMappingOptions(desc);
  ActsExamples::Options::addPropagationOptions(desc);
  ActsExamples::Options::addInputOptions(desc);
  ActsExamples::Options::addOutputOptions(desc,
                                          ActsExamples::OutputFormat::Root |
                                              ActsExamples::OutputFormat::Json |
                                              ActsExamples::OutputFormat::Cbor);

  // Add specific options for this geometry
  detector.addOptions(desc);
  auto vm = ActsExamples::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  ActsExamples::Sequencer sequencer(
      ActsExamples::Options::readSequencerConfig(vm));

  // Get the log level
  auto logLevel = ActsExamples::Options::readLogLevel(vm);

  // The geometry, material and decoration
  auto geometry = ActsExamples::Geometry::build(vm, detector);
  auto tGeometry = geometry.first;
  auto contextDecorators = geometry.second;
  for (auto cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }

  ActsExamples::WhiteBoard contextBoard(
      Acts::getDefaultLogger("contextBoard", logLevel));

  // The geometry context
  ActsExamples::AlgorithmContext context(0, 0, contextBoard);

  /// Decorate the context
  for (auto& cdr : contextDecorators) {
    if (cdr->decorate(context) != ActsExamples::ProcessCode::SUCCESS)
      throw std::runtime_error("Failed to decorate event context");
  }

  /// Default contexts
  Acts::GeometryContext geoContext = context.geoContext;
  Acts::MagneticFieldContext mfContext = context.magFieldContext;

  // Straight line stepper
  using SlStepper = Acts::StraightLineStepper;
  using Propagator = Acts::Propagator<SlStepper, Acts::Navigator>;

  auto mapSurface = vm["mat-mapping-surfaces"].template as<bool>();
  auto mapVolume = vm["mat-mapping-volumes"].template as<bool>();
  auto volumeStep = vm["mat-mapping-volume-stepsize"].template as<float>();
  if (!mapSurface && !mapVolume) {
    return EXIT_FAILURE;
  }
  // ---------------------------------------------------------------------------------
  // Input directory & input file handling
  std::string intputDir = vm["input-dir"].template as<std::string>();
  auto intputFiles = vm["input-files"].template as<std::vector<std::string>>();

  if (vm["input-root"].template as<bool>()) {
    // Read the material step information from a ROOT TTree
    ActsExamples::RootMaterialTrackReader::Config matTrackReaderRootConfig;
    matTrackReaderRootConfig.collection =
        ActsExamples::Simulation::kMaterialTracks;
    matTrackReaderRootConfig.fileList = intputFiles;
    auto matTrackReaderRoot =
        std::make_shared<ActsExamples::RootMaterialTrackReader>(
            matTrackReaderRootConfig, logLevel);
    sequencer.addReader(matTrackReaderRoot);
  }

  /// The material mapping algorithm
  ActsExamples::MaterialMapping::Config mmAlgConfig{geoContext, mfContext};
  mmAlgConfig.collection = ActsExamples::Simulation::kMaterialTracks;
  if (mapSurface) {
    // Get a Navigator
    Acts::Navigator navigator({tGeometry, true, true, true});
    // Make stepper and propagator
    SlStepper stepper;
    Propagator propagator(std::move(stepper), std::move(navigator));
    /// The material surface mapper
    Acts::SurfaceMaterialMapper::Config smmConfig;
    auto smm = std::make_shared<Acts::SurfaceMaterialMapper>(
        smmConfig, std::move(propagator),
        Acts::getDefaultLogger("SurfaceMaterialMapper", logLevel));
    mmAlgConfig.materialSurfaceMapper = smm;
  }
  if (mapVolume) {
    // Get a Navigator
    Acts::Navigator navigator({tGeometry});
    // Make stepper and propagator
    SlStepper stepper;
    Propagator propagator(std::move(stepper), std::move(navigator));
    /// The material volume mapper
    Acts::VolumeMaterialMapper::Config vmmConfig;
    vmmConfig.mappingStep = volumeStep;
    auto vmm = std::make_shared<Acts::VolumeMaterialMapper>(
        vmmConfig, std::move(propagator),
        Acts::getDefaultLogger("VolumeMaterialMapper", logLevel));
    mmAlgConfig.materialVolumeMapper = vmm;
  }
  mmAlgConfig.trackingGeometry = tGeometry;

  // Get the file name from the options
  std::string materialFileName = vm["mat-output-file"].as<std::string>();

  if (not materialFileName.empty() and vm["output-root"].template as<bool>()) {
    // The writer of the indexed material
    ActsExamples::RootMaterialWriter::Config rmwConfig;
    rmwConfig.filePath = materialFileName + ".root";
    // Fullfill the IMaterialWriter interface

    auto rmw =
        std::make_shared<ActsExamples::RootMaterialWriter>(rmwConfig, logLevel);
    mmAlgConfig.materialWriters.push_back(rmw);

    if (mapSurface) {
      // Write the propagation steps as ROOT TTree
      ActsExamples::RootMaterialTrackWriter::Config matTrackWriterRootConfig;
      matTrackWriterRootConfig.filePath = materialFileName + "_tracks.root";
      matTrackWriterRootConfig.collection =
          mmAlgConfig.mappingMaterialCollection;
      matTrackWriterRootConfig.storeSurface = true;
      matTrackWriterRootConfig.storeVolume = true;
      auto matTrackWriterRoot =
          std::make_shared<ActsExamples::RootMaterialTrackWriter>(
              matTrackWriterRootConfig, logLevel);
      sequencer.addWriter(matTrackWriterRoot);
    }
  }

  if (!materialFileName.empty() and (vm["output-json"].template as<bool>() or
                                     vm["output-cbor"].template as<bool>())) {
    /// The name of the output file
    std::string fileName = vm["mat-output-file"].template as<std::string>();
    // the material writer
    Acts::MaterialMapJsonConverter::Config jmConverterCfg;
    jmConverterCfg.processSensitives =
        vm["mat-output-sensitives"].template as<bool>();
    jmConverterCfg.processApproaches =
        vm["mat-output-approaches"].template as<bool>();
    jmConverterCfg.processRepresenting =
        vm["mat-output-representing"].template as<bool>();
    jmConverterCfg.processBoundaries =
        vm["mat-output-boundaries"].template as<bool>();
    jmConverterCfg.processVolumes =
        vm["mat-output-volumes"].template as<bool>();
    jmConverterCfg.context = geoContext;
    // The writer
    ActsExamples::JsonMaterialWriter::Config jmWriterCfg;
    jmWriterCfg.converterCfg = std::move(jmConverterCfg);
    jmWriterCfg.fileName = materialFileName;

    ActsExamples::JsonFormat format = ActsExamples::JsonFormat::NoOutput;
    if (vm["output-json"].template as<bool>()) {
      format = format | ActsExamples::JsonFormat::Json;
    }
    if (vm["output-cbor"].template as<bool>()) {
      format = format | ActsExamples::JsonFormat::Cbor;
    }
    jmWriterCfg.writeFormat = format;

    auto jmw = std::make_shared<ActsExamples::JsonMaterialWriter>(
        std::move(jmWriterCfg), Acts::Logging::INFO);

    mmAlgConfig.materialWriters.push_back(jmw);
  }

  // Create the material mapping
  auto mmAlg = std::make_shared<ActsExamples::MaterialMapping>(mmAlgConfig);

  // Append the Algorithm
  sequencer.addAlgorithm(mmAlg);

  // Initiate the run
  sequencer.run();
  // Return success code
  return 0;
}
