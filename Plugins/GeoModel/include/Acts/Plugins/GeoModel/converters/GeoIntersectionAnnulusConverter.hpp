// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Plugins/GeoModel/detail/GenericGeoShapeConverter.hpp"
#include "Acts/Utilities/Result.hpp"

#include <memory>
#include <tuple>

#include <GeoModelKernel/GeoShapeIntersection.h>

namespace Acts {

namespace detail {
struct GeoIntersectionAnnulusConverter {
  /// @brief Convert a GeoBox to a detector element and surface
  ///
  /// @param geoFPV The full physical volume to convert (contains shape)
  /// @param geoIntersection The GeoIntersection to convert
  /// @param absTransform from the GeoPhysVol
  /// @param bool sensitive
  ///
  /// @return The detector element and surface
  std::tuple<std::shared_ptr<GeoModelDetectorElement>, std::shared_ptr<Surface>>
  operator()(const GeoFullPhysVol& geoFPV,
            const GeoShapeIntersection& geoIntersection,
            const Transform3& absTransform, bool sensitive) const;
};
}

/// @brief A dedicated converter for GeoInterseciton that describe annulus bounds
///
/// This is very much tailored to the AnnulusBounds class
using GeoIntersectionAnnulusConverter = detail::GenericGeoShapeConverter<GeoShapeIntersection, detail::GeoIntersectionAnnulusConverter>;

}  // namespace Acts
