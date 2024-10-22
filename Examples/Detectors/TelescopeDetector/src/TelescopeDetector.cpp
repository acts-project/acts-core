// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "ActsExamples/TelescopeDetector/BuildTelescopeDetector.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetectorElement.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>

namespace ActsExamples::Telescope {

TelescopeDetector::TelescopeDetector(const Config& cfg)
    : DetectorCommons::Detector(
          Acts::getDefaultLogger("TelescopeDetector", cfg.logLevel)),
      m_cfg(cfg) {}

void TelescopeDetector::buildTrackingGeometry() {
  Acts::GeometryContext gctx;

  if (m_cfg.surfaceType > 1) {
    throw std::invalid_argument(
        "The surface type could either be 0 for plane surface or 1 for disc "
        "surface.");
  }
  if (m_cfg.binValue > 2) {
    throw std::invalid_argument("The axis value could only be 0, 1, or 2.");
  }
  // Check if the bounds values are valid
  if (m_cfg.surfaceType == 1 && m_cfg.bounds[0] >= m_cfg.bounds[1]) {
    throw std::invalid_argument(
        "The minR should be smaller than the maxR for disc surface bounds.");
  }

  if (m_cfg.positions.size() != m_cfg.stereos.size()) {
    throw std::invalid_argument(
        "The number of provided positions must match the number of "
        "provided stereo angles.");
  }

  // Sort the provided distances
  std::vector<double> positions = m_cfg.positions;
  std::vector<double> stereos = m_cfg.stereos;
  std::sort(positions.begin(), positions.end());

  /// Return the telescope detector
  m_trackingGeometry = ActsExamples::Telescope::buildDetector(
      gctx, m_detectorStore, positions, stereos, m_cfg.offsets, m_cfg.bounds,
      m_cfg.thickness, static_cast<TelescopeSurfaceType>(m_cfg.surfaceType),
      static_cast<Acts::BinningValue>(m_cfg.binValue));
}

}  // namespace ActsExamples::Telescope
