// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "ActsExamples/Io/Json/JsonTGeoDetectorConfig.hpp"

#include <memory>
#include <vector>

#include <boost/program_options.hpp>

namespace ActsExamples {
namespace TGeo {

/// @brief global method to build the generic tracking geometry
// from a TGeo object.
///
/// It does *currently* not translate the material, this has
/// to be done with a material mapping stage
///
/// @tparam variable_map_t is the variable map
///
/// @param vm is the variable map from the options
std::shared_ptr<const Acts::TrackingGeometry> buildTGeoDetector(
    const boost::program_options::variables_map& vm,
    const Acts::GeometryContext& context,
    std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>&
        detElementStore,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator);

}  // namespace TGeo
}  // namespace ActsExamples
