// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <iostream>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace ActsExamples {

namespace Options {

/// @brief Material mapping options, specially added
///
/// @tparam aopt_t Type of the options object (API bound to boost)
///
/// @param [in] opt_t The options object where the specific digitization
/// options are attached to
template <typename aopt_t>
void addMaterialMappingOptions(aopt_t& opt) {
  opt.add_options()(
      "mat-mapping-emptybins", po::value<bool>()->default_value(true),
      "Empty bin correction (recommended). Corrects for vaccuum/emtpy "
      "assigments.")("mat-mapping-surfaces",
                     po::value<bool>()->default_value(true),
                     "Map the material onto the selected surfaces")(
      "mat-mapping-volumes", po::value<bool>()->default_value(false),
      "Map the material into the selected material volumes")(
      "mat-mapping-volume-stepsize",
      po::value<float>()->default_value(std::numeric_limits<float>::infinity()),
      "Step size of the sampling of volume material for the mapping "
      "(should be smaller than the size of the bins in depth)");
}

}  // namespace Options
}  // namespace ActsExamples
