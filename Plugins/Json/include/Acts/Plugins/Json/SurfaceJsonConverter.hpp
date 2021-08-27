// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"
#include "Acts/Plugins/Json/SurfaceBoundsJsonConverter.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <string>
#include <vector>

// Custom Json encoder/decoders. Naming is mandated by nlohman::json and thus
// can not match our naming guidelines.
namespace Acts {

using SurfaceAndMaterialWithContext =
    std::tuple<std::shared_ptr<const Acts::Surface>,
               std::shared_ptr<const Acts::ISurfaceMaterial>,
               Acts::GeometryContext>;

static std::vector<std::string> surfaceTypes = {
    "ConeSurface",  "CylinderSurface", "DiscSurface",       "PerigeeSurface",
    "PlaneSurface", "StrawSurface",    "CurvilinearSurface"};

/// Conversion of a pair of surface and material used for the material mapping
void to_json(nlohmann::json& j, const SurfaceAndMaterialWithContext& surface);

/// Non-contextual conversion of a surface
///
/// @note it will take the default context
void to_json(nlohmann::json& j, const Surface& surface);

/// Contextual conversion of a surface
///
/// @param j the json to be filled
/// @param surface the surface to be converted
/// @param gctx the geometry context for this
void toJson(nlohmann::json& j, const Surface& surface,
            const Acts::GeometryContext& gctx);

/// Non-contextual conversion of a surface
///
/// @note it will take the default context
void to_json(nlohmann::json& j, std::shared_ptr<const Surface> surface);

/// Contextual conversion of a surface
///
/// @param j the json to be filled
/// @param surface the surface to be converted
/// @param gctx the geometry context for this
void toJson(nlohmann::json& j, std::shared_ptr<const Surface> surface,
            const Acts::GeometryContext& gctx);

/// Converstion to Surface from jsonn
///
/// @param j the read-in json object
///
/// @return a shared_ptr to a surface object for type polymorphism
std::shared_ptr<Surface> surfaceFromJson(const nlohmann::json& j);

/// Converstion to Surface from json in correct type
///
/// The type is given as a template argument in order to be able
/// to construct the correct fitting types for surfaces.
///
/// @param j the read-in json object
///
/// @return a shared_ptr to a typed surface object for type polymorphism
template <typename surface_t, typename bounds_t>
std::shared_ptr<surface_t> surfaceFromJsonT(const nlohmann::json& j) {
  Transform3 sTransform;
  nlohmann::json jtrf = j["transform"];
  from_json(jtrf, sTransform);
  nlohmann::json jbounds = j["bounds"];
  auto sBounds = surfaceBoundsFromJson<bounds_t>(jbounds);
  return Surface::makeShared<surface_t>(sTransform, std::move(sBounds));
}

}  // namespace Acts
