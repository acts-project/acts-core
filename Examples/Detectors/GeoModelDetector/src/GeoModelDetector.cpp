// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GeoModelDetector/GeoModelDetector.hpp"

#include "Acts/Plugins/GeoModel/GeoModelReader.hpp"
#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"

namespace ActsExamples {

GeoModelDetector::GeoModelDetector(const Config& cfg)
    : Detector(Acts::getDefaultLogger("GeoModelDetector", cfg.logLevel)),
      m_cfg(cfg) {
  m_geoModel = Acts::GeoModelReader::readFromDb(m_cfg.path);
}

}  // namespace ActsExamples
