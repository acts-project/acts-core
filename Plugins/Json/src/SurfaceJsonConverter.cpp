// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/SurfaceJsonConverter.hpp"

#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Plugins/Json/MaterialJsonConverter.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

void Acts::to_json(nlohmann::json& j,
                   const Acts::SurfaceAndMaterialWithContext& surface) {
  toJson(j, std::get<0>(surface), std::get<2>(surface));
  to_json(j, std::get<1>(surface).get());
}

void Acts::to_json(nlohmann::json& j, const Acts::Surface& surface) {
  Acts::GeometryContext gctx;
  toJson(j, surface, gctx);
}

void Acts::to_json(nlohmann::json& j,
                   std::shared_ptr<const Acts::Surface> surface) {
  Acts::GeometryContext gctx;
  toJson(j, *(surface.get()), gctx);
}

void Acts::toJson(nlohmann::json& j,
                  std::shared_ptr<const Acts::Surface> surface,
                  const Acts::GeometryContext& gctx) {
  toJson(j, *(surface.get()), gctx);
}

void Acts::toJson(nlohmann::json& j, const Acts::Surface& surface,
                  const Acts::GeometryContext& gctx) {
  const auto& sBounds = surface.bounds();
  const auto sTransform = surface.transform(gctx);
  j["bounds"] = nlohmann::json(sBounds);
  nlohmann::json trfj;
  to_json(trfj, sTransform);
  j["transform"] = trfj;
  j["type"] = surfaceTypes[surface.type()];
  j["geo_id"] = surface.geometryId().value();
  if (surface.surfaceMaterial() != nullptr) {
    j["material"] = nlohmann::json(surface.surfaceMaterial());
  }
}
std::shared_ptr<Acts::Surface> Acts::surfaceFromJson(const nlohmann::json& j) {
  std::string sType = j["type"];
  std::string bType = j["bounds"]["type"];

  std::shared_ptr<Acts::Surface> mutableSf = nullptr;

  /// Unroll the types
  if (sType == "PlaneSurface") {
    if (bType == "EllipseBounds") {
      mutableSf = surfaceFromJsonT<Acts::PlaneSurface, Acts::EllipseBounds>(j);
    } else if (bType == "RectangleBounds") {
      mutableSf =
          surfaceFromJsonT<Acts::PlaneSurface, Acts::RectangleBounds>(j);
    } else if (bType == "TrapezoidBounds") {
      mutableSf =
          surfaceFromJsonT<Acts::PlaneSurface, Acts::TrapezoidBounds>(j);
    }
  } else if (sType == "DiscSurface") {
    if (bType == "AnnulusBounds") {
      mutableSf = surfaceFromJsonT<Acts::DiscSurface, Acts::AnnulusBounds>(j);
    } else if (bType == "RadialBounds") {
      mutableSf = surfaceFromJsonT<Acts::DiscSurface, Acts::RadialBounds>(j);
    } else if (bType == "DiscTrapezoidBounds") {
      mutableSf =
          surfaceFromJsonT<Acts::DiscSurface, Acts::DiscTrapezoidBounds>(j);
    }
  } else if (sType == "CylinderSurface") {
    mutableSf =
        surfaceFromJsonT<Acts::CylinderSurface, Acts::CylinderBounds>(j);
  } else if (sType == "ConeSurface") {
    mutableSf = surfaceFromJsonT<Acts::ConeSurface, Acts::ConeBounds>(j);
  } else if (sType == "StrawSurface") {
    mutableSf = surfaceFromJsonT<Acts::StrawSurface, Acts::LineBounds>(j);
  } else if (sType == "PerigeeSurface") {
    Transform3 pTransform;
    nlohmann::json trfj = j["transform"];
    ;
    from_json(trfj, pTransform);
    mutableSf = Surface::makeShared<PerigeeSurface>(pTransform);
  }

  if (mutableSf != nullptr) {
    GeometryIdentifier geoID(j["geo_id"]);
    mutableSf->assignGeometryId(geoID);
    // Add material
    if (j.find("material") != j.end() and not j["material"].empty()) {
      const Acts::ISurfaceMaterial* surfaceMaterial;
      from_json(j, surfaceMaterial);
      std::shared_ptr<const ISurfaceMaterial> sharedSurfaceMaterial(
          surfaceMaterial);
      mutableSf->assignSurfaceMaterial(sharedSurfaceMaterial);
    }
  }
  return mutableSf;
}
