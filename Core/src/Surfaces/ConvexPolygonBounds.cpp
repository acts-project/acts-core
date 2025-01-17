// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/ConvexPolygonBounds.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"

#include <optional>
#include <ostream>

namespace Acts {

std::ostream& ConvexPolygonBoundsBase::toStream(std::ostream& sl) const {
  std::vector<Vector2> vtxs = vertices();
  sl << "Acts::ConvexPolygonBounds<" << vtxs.size() << ">: vertices: [x, y]\n";
  for (std::size_t i = 0; i < vtxs.size(); i++) {
    const auto& vtx = vtxs[i];
    if (i > 0) {
      sl << ",";
      sl << "\n";
    }
    sl << "[" << vtx.x() << ", " << vtx.y() << "]";
  }
  return sl;
}

std::vector<double> ConvexPolygonBoundsBase::values() const {
  std::vector<double> values;
  for (const auto& vtx : vertices()) {
    values.push_back(vtx.x());
    values.push_back(vtx.y());
  }
  return values;
}

ConvexPolygonBounds<PolygonDynamic>::ConvexPolygonBounds(
    const std::vector<Vector2>& vertices)
    : m_vertices(vertices.begin(), vertices.end()),
      m_boundingBox(makeBoundingBox(vertices)) {}

bool ConvexPolygonBounds<PolygonDynamic>::inside(
    const Vector2& lposition) const {
  return detail::VerticesHelper::isInsidePolygon(lposition, m_vertices);
}

Vector2 ConvexPolygonBounds<PolygonDynamic>::closestPoint(
    const Vector2& lposition,
    const std::optional<SquareMatrix2>& metric) const {
  return detail::VerticesHelper::computeClosestPointOnPolygon(
      lposition, std::span<const Vector2>(m_vertices.data(), m_vertices.size()),
      metric.value_or(SquareMatrix2::Identity()));
}

std::vector<Vector2> ConvexPolygonBounds<PolygonDynamic>::vertices(
    unsigned int /*lseg*/) const {
  return {m_vertices.begin(), m_vertices.end()};
}

const RectangleBounds& ConvexPolygonBounds<PolygonDynamic>::boundingBox()
    const {
  return m_boundingBox;
}

void ConvexPolygonBounds<PolygonDynamic>::checkConsistency() const
    noexcept(false) {
  convex_impl(m_vertices);
}

}  // namespace Acts
