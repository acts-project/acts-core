// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/ProtoLayer.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

namespace Acts::detail {

void ProtoLayerBase::measureImpl(const GeometryContext& gctx,
                                 const std::vector<const Surface*>& surfaces,
                                 Extent& extent, const Transform3& transform) {
  for (const auto& sf : surfaces) {
    // To prevent problematic isInsidePolygon check for straw surfaces with only
    // one lseg
    int lseg = (sf->type() != Surface::Straw) ? 1 : 2;
    auto sfPolyhedron = sf->polyhedronRepresentation(gctx, lseg);
    const DetectorElementBase* element = sf->associatedDetectorElement();
    const auto* regSurface = dynamic_cast<const RegularSurface*>(sf);
    if (element != nullptr && regSurface != nullptr) {
      // Take the thickness in account if necessary
      double thickness = element->thickness();
      // We need a translation along and opposite half thickness
      Vector3 sfNormal = regSurface->normal(gctx, sf->center(gctx));
      for (const auto& dT : {-0.5 * thickness, 0.5 * thickness}) {
        Transform3 dtransform = transform * Translation3{dT * sfNormal};
        extent.extend(sfPolyhedron.extent(dtransform));
      }
      continue;
    }
    extent.extend(sfPolyhedron.extent(transform));
  }
}

double ProtoLayerBase::min(AxisDirection aDir, bool addenv) const {
  if (addenv) {
    return extent.min(aDir) - envelope[aDir][0u];
  }
  return extent.min(aDir);
}

double ProtoLayerBase::max(AxisDirection aDir, bool addenv) const {
  if (addenv) {
    return extent.max(aDir) + envelope[aDir][1u];
  }
  return extent.max(aDir);
}

double ProtoLayerBase::medium(AxisDirection aDir, bool addenv) const {
  return 0.5 * (min(aDir, addenv) + max(aDir, addenv));
}

double ProtoLayerBase::range(AxisDirection aDir, bool addenv) const {
  return std::abs(max(aDir, addenv) - min(aDir, addenv));
}

std::ostream& ProtoLayerBase::toStream(std::ostream& sl) const {
  sl << "ProtoLayer with dimensions (min/max)" << std::endl;
  sl << extent.toString();
  return sl;
}

}  // namespace Acts::detail
