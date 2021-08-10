// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/TGeo/TGeoCylinderDiscSplitter.hpp"
#include "Acts/Plugins/TGeo/TGeoITkModuleSplitter.hpp"
#include "Acts/Plugins/TGeo/TGeoLayerBuilder.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"


#include <nlohmann/json.hpp>
#include <vector>

namespace ActsExamples {
namespace Options {

/// Add TGeo detector options prefixed with geo-tgeo.
void addTGeoGeometryOptions(Description& desc);

/// Read the BeamPipe configuration from the user configuration.
std::vector<double> readBeampipeBuilderParam(const std::string& path);

/// Read config for cylinder/disc module splitter
void from_json(const nlohmann::json& j,
                             Acts::TGeoCylinderDiscSplitter::Config& msc);

/// Read config for dedicated itk module splitter
void from_json(const nlohmann::json& j,
                             Acts::TGeoITkModuleSplitter::Config& msc);

/// Read nested layer config for TGeoLayerbuilder
void from_json(const nlohmann::json& j,
                             Acts::TGeoLayerBuilder::LayerConfig& psc);

/// Read the TGeo layer builder configurations from the user configuration.
std::vector<Acts::TGeoLayerBuilder::Config> readTGeoLayerBuilderConfigs(const std::string& path);

}  // namespace Options
}  // namespace ActsExamples
