// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/DigitizationCoordinatesConverter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"

#include <stdexcept>
#include <utility>

namespace ActsExamples {

DigitizationCoordinatesConverter::DigitizationCoordinatesConverter(
    DigitizationAlgorithm::Config config)
    : m_cfg(std::move(config)) {
  if (m_cfg.surfaceByIdentifier.empty()) {
    throw std::invalid_argument("Missing Surface-GeometryID association map");
  }
}

std::tuple<double, double> DigitizationCoordinatesConverter::globalToLocal(
    std::uint64_t moduleId, double x, double y, double z) const {
  const Acts::GeometryIdentifier moduleGeoId(moduleId);
  auto surfaceItr = m_cfg.surfaceByIdentifier.find(moduleGeoId);
  if (surfaceItr == m_cfg.surfaceByIdentifier.end()) {
    throw std::runtime_error("Surface not found for moduleGeoId");
  }

  // Transform the position into the local surface frame
  const Acts::Surface* surfacePtr = surfaceItr->second;
  const auto& invTransform =
      surfacePtr->transform(Acts::GeometryContext()).inverse();

  const Acts::Vector3 pos(x, y, z);
  Acts::Vector2 pos2Local = (invTransform * pos).segment<2>(0);

  return {pos2Local.x(), pos2Local.y()};
}

std::tuple<double, double, double>
DigitizationCoordinatesConverter::localToGlobal(std::uint64_t moduleId,
                                                double x, double y) const {
  const Acts::GeometryIdentifier moduleGeoId{moduleId};
  auto surfaceItr = m_cfg.surfaceByIdentifier.find(moduleGeoId);
  if (surfaceItr == m_cfg.surfaceByIdentifier.end()) {
    throw std::runtime_error("Surface not found for moduleGeoId");
  }

  // Transform the position into the global frame
  const Acts::Surface* surfacePtr = surfaceItr->second;
  const auto& transform = surfacePtr->transform(Acts::GeometryContext());

  const Acts::Vector3 pos(x, y, 0);
  Acts::Vector3 pos2Global = (transform * pos);

  return {pos2Global.x(), pos2Global.y(), pos2Global.z()};
}

}  // namespace ActsExamples
