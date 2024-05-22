// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

class GeoVPhysVol;

#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"

#include <string>

namespace Acts::GeoModelReader {

/// @brief Read the GeoModel from the database
///
/// @param dbPath path to the database
///
/// @return world/top volume of the GeoModel tree in memory
GeoModelTree readFromDb(const std::string& dbPath);

}  // namespace Acts::GeoModelReader
