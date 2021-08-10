// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Json/JsonTGeoDetectorConfig.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <cstdlib>
#include <fstream>
#include <vector>

#include <boost/program_options.hpp>

void ActsExamples::Options::addTGeoGeometryOptions(Description& desc) {
  using boost::program_options::value;

  // The options are encoded in a json file, where an empty version containing
  // the required structure can be dumped.
  //
  //   # Unit scalor from ROOT to Acts.
  //   "geo-tgeo-unit-scalor": 1,
  //   # Root world volume to start search from.
  //   "geo-tgeo-worldvolume": "IDET",
  //   # Beam pipe parameters {r, z, t} in [mm]. Beam pipe is automatically
  //     created if the parameters are present.
  //   "geo-tgeo-beampipe-parameters": [29.0,3000.0,0.8]
  //
  // Each detector volume configuration is one logical block which can
  // be repeated as many times as there are usable detector volumes.
  //
  // required per-volume options:
  //
  //   # Detector volume name
  //   "geo-tgeo-volume": "InnerPixels",
  //   # vTolerance interval in r [mm] for automated surface binninng.
  //   "geo-tgeo-sfbin-r-tolerance": [5,5],
  //   # Tolerance interval in phi [rad] for automated surface binning.
  //   "geo-tgeo-sfbin-phi-tolerance": [0.025,0.025],
  //   # Tolerance interval in z [mm] for automated surface binning.
  //   "geo-tgeo-sfbin-z-tolerance": [5,5],
  //   "geo-tgeo-nlayers" false,  # boolean switch whether there are negative layers
  //   "geo-tgeo-clayers" true,  # boolean switch whether there are central layers
  //   "geo-tgeo-players" false,  # boolean switch whether there are positive layers
  //
  // within each volume there can be negative/central/positive layers depending
  // on the which `geo-tgeo-{n,c,p}layers` flags are set to true. If any of
  // them are set, they must be followed by the corresponding layer options. If
  // the `*layers` option is false, an empty json object must be present at the
  // corresponding place in the "Layers" list to keep the ncp association correct.
  //
  // examples: negative and central layers, but not positive layers
  //
  //   "geo-tgeo-nlayers": true,  # Are layers on the negative side?
  //   "geo-tgeo-clayers": true,  # Are layers in the central region?
  //   "geo-tgeo-players": false, # Are layers on the positive side?
  //   "Layers": [
  //       # Config for layers on the negative side
  //       {
  //           # Name identifier of the volume for searching negative layers.
  //           "geo-tgeo-volume-name": "Pixel::Pixel",
  //
  //           # Name identifier for negative sensitive objects
  //           "geo-tgeo-module-name": ["Pixel::siLog"],
  //
  //           # Axes definition for negative sensitive objects.
  //           "geo-tgeo-module-axes": "YZX",
  //
  //           # Radial range(s) for negative layers to restrict the module parsing.
  //           "geo-tgeo-layer-r-range": [0,135],
  //
  //           # R-tolerances (if > 0.) that triggers splitting
  //             of collected surfaces into different negative layers.
  //           "geo-tgeo-layer-z-range": [-3000,-250],
  //
  //           # Longitudinal range(s) for negative layers to restrict the module
  //             parsing.
  //           "geo-tgeo-layer-r-split": -1.0,
  //
  //           # Z-tolerances (if > 0.) that triggers splitting of collected 
  //             surfaces into different negative layers.
  //           "geo-tgeo-layer-z-split": 10.0
  //       },
  //       # Config for layers in the center
  //       {
  //           "geo-tgeo-volume-name": "Pixel::Pixel",
  //           "geo-tgeo-module-name": ["Pixel::siLog"],
  //           "geo-tgeo-module-axes": "YZX",
  //           "geo-tgeo-layer-r-range": [0,135],
  //           "geo-tgeo-layer-z-range": [-250,250],
  //           "geo-tgeo-layer-r-split": 5.0,
  //           "geo-tgeo-layer-z-split": -1.0
  //       },
  //       # Config for layers on the positive side
  //       {
  //       }
  //   ]
  //
  //   "CylinderDisk" # Splitter of type TGeoCylinderDiscSplitter
  //
  //
  //  In case cylinder / disc splitting is on:
  //
  //   "geo-tgeo-cyl-nz-segs"    # number of z segments for cylinder splitting
  //   "geo-tgeo-cyl-nphi-segs"  # number of phi segments for cylinder splitting
  //   "geo-tgeo-disc-nr-segs"   # number of r segments for disc splitting
  //   "geo-tgeo-disc-nphi-segs" # number of phi segments for disc splitting

  auto opt = desc.add_options();
  // required global options
  opt("geo-tgeo-filename", value<std::string>()->default_value(""),
      "Root file name.");
  opt("geo-tgeo-jsonconfig", value<std::string>()->default_value(""),
      "Json config file name.");
}

std::vector<double> ActsExamples::Options::readBeampipeBuilderParam(const std::string& path) {

  nlohmann::json djson;
  if (path.empty()) {
    return {};
  }
  std::ifstream infile(path, std::ifstream::in | std::ifstream::binary);
  // rely on exception for error handling
  infile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  infile >> djson;

  return djson["geo-tgeo-beampipe-parameters"].get<std::vector<double>>();
}

void ActsExamples::Options::from_json(const nlohmann::json& j,
                             Acts::TGeoCylinderDiscSplitter::Config& msc) {

  /// Number of segments in phi for a disc
  msc.cylinderPhiSegments = j["geo-tgeo-cyl-nphi-segs"];
  /// Number of segments in r for a disk
  msc.cylinderLongitudinalSegments = j["geo-tgeo-cyl-nz-segs"];
  /// Number of segments in phi for a disc
  msc.discPhiSegments = j["geo-tgeo-disc-nphi-segs"];
  /// Number of segments in r for a disk
  msc.discRadialSegments = j["geo-tgeo-disc-nr-segs"];
}

void ActsExamples::Options::from_json(const nlohmann::json& j,
                             Acts::TGeoITkModuleSplitter::Config& msc) {

  msc.barrelMap = j["geo-tgeo-barrel-map"].get<std::map<std::string, unsigned int>>();
  msc.discMap = j["geo-tgeo-disc-map"].get<std::map<std::string, std::vector<std::pair<double, double>>>>();
}

void ActsExamples::Options::from_json(const nlohmann::json& j,
                             Acts::TGeoLayerBuilder::LayerConfig& psc) {

  psc.volumeName  = j["geo-tgeo-volume-name"];
  psc.sensorNames = j["geo-tgeo-module-name"].get<std::vector<std::string>>();
  psc.localAxes   = j["geo-tgeo-module-axes"];
  auto r_range = j["geo-tgeo-layer-r-range"].get<std::pair<double, double>>();
  auto z_range = j["geo-tgeo-layer-z-range"].get<std::pair<double, double>>();
  psc.parseRanges = {{Acts::binR, r_range}, {Acts::binZ, z_range}};
  double r_split = j["geo-tgeo-layer-r-split"];
  double z_split = j["geo-tgeo-layer-z-split"];
  if (0 < r_split) {
      psc.splitConfigs.emplace_back(Acts::binR, r_split);
  }
  if (0 < z_split) {
      psc.splitConfigs.emplace_back(Acts::binZ, z_split);
  }
}

std::vector<Acts::TGeoLayerBuilder::Config>
ActsExamples::Options::readTGeoLayerBuilderConfigs(const std::string& path) {
  std::vector<Acts::TGeoLayerBuilder::Config> detLayerConfigs = {};

  nlohmann::json djson;
  if (path.empty()) {
    return detLayerConfigs;
  }
  std::ifstream infile(path, std::ifstream::in | std::ifstream::binary);
  // rely on exception for error handling
  infile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  infile >> djson;

  double unitScalor = djson["geo-tgeo-unit-scalor"];
  std::string worldVolume = djson["geo-tgeo-worldvolume"];

  for (const auto& volume : djson["LayerConfigs"]["Volumes"]) {
    Acts::TGeoLayerBuilder::Config layerBuilderConfig;
    // subdetector selection
    std::string subDetector = volume["geo-tgeo-volume"];
    layerBuilderConfig.configurationName = subDetector;
    layerBuilderConfig.unit = unitScalor;

    // configure surface autobinning
    auto binToleranceR =
      volume["geo-tgeo-sfbin-r-tolerance"].get<std::vector<double>>();
    auto binToleranceZ =
      volume["geo-tgeo-sfbin-z-tolerance"].get<std::vector<double>>();
    auto binTolerancePhi =
      volume["geo-tgeo-sfbin-phi-tolerance"].get<std::vector<double>>();
          std::vector<std::pair<double, double>> binTolerances(
        static_cast<size_t>(Acts::binValues), {0., 0.});
    binTolerances[Acts::binR] = {binToleranceR[0],
                                 binToleranceR[1]};
    binTolerances[Acts::binZ] = {binToleranceZ[0],
                                 binToleranceZ[1]};
    binTolerances[Acts::binPhi] = {binTolerancePhi[0],
                                   binTolerancePhi[1]};
    layerBuilderConfig.autoSurfaceBinning = true;
    layerBuilderConfig.surfaceBinMatcher =
        Acts::SurfaceBinningMatcher(binTolerances);

    auto isLayers = std::vector<bool>{volume["geo-tgeo-nlayers"], 
                                      volume["geo-tgeo-clayers"],
                                      volume["geo-tgeo-players"]};
    size_t ncp = 0;
    for (const auto& layer : volume["Layers"]) {
      if (not isLayers[ncp]) {
        ncp++;
        continue;
      }
      Acts::TGeoLayerBuilder::LayerConfig lConfig;
      from_json(layer, lConfig);
      layerBuilderConfig.layerConfigurations[ncp++].push_back(lConfig);
    }

    for(const auto& splitter : volume["Splitters"]) {
      for (const auto& cdSplitter : splitter["CylinderDisk"]) {
        Acts::TGeoCylinderDiscSplitter::Config cdConfig;
        from_json(cdSplitter, cdConfig);
        layerBuilderConfig.detectorElementSplitter = std::make_shared<Acts::TGeoCylinderDiscSplitter>(cdConfig);
      }

      for (const auto& itkSplitter : splitter["Itk"]) {
        Acts::TGeoITkModuleSplitter::Config itkConfig;
        from_json(itkSplitter, itkConfig);
        layerBuilderConfig.detectorElementSplitter = std::make_shared<Acts::TGeoITkModuleSplitter>(itkConfig);
      }
    }

    detLayerConfigs.push_back(layerBuilderConfig);
  }

  return detLayerConfigs;
}
