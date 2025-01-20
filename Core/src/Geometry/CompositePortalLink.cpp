// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CompositePortalLink.hpp"

#include "Acts/Geometry/GridPortalLink.hpp"
#include "Acts/Geometry/PortalError.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Utilities/Axis.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <stdexcept>

#include <boost/algorithm/string/join.hpp>

namespace Acts {

namespace {
std::shared_ptr<RegularSurface> mergedSurface(const Surface& a,
                                              const Surface& b,
                                              AxisDirection direction) {
  if (a.type() != b.type()) {
    throw std::invalid_argument{"Cannot merge surfaces of different types"};
  }

  if (const auto* cylinderA = dynamic_cast<const CylinderSurface*>(&a);
      cylinderA != nullptr) {
    const auto& cylinderB = dynamic_cast<const CylinderSurface&>(b);

    auto [merged, reversed] = cylinderA->mergedWith(cylinderB, direction, true);
    return merged;
  } else if (const auto* discA = dynamic_cast<const DiscSurface*>(&a);
             discA != nullptr) {
    const auto& discB = dynamic_cast<const DiscSurface&>(b);
    auto [merged, reversed] = discA->mergedWith(discB, direction, true);
    return merged;
  } else if (const auto* planeA = dynamic_cast<const PlaneSurface*>(&a);
             planeA != nullptr) {
    throw std::logic_error{"Plane surfaces not implemented yet"};
  } else {
    throw std::invalid_argument{"Unsupported surface type"};
  }
}

std::shared_ptr<RegularSurface> mergePortalLinks(
    const std::vector<std::unique_ptr<PortalLinkBase>>& links,
    AxisDirection direction) {
  assert(std::ranges::all_of(links,
                             [](const auto& link) { return link != nullptr; }));
  assert(!links.empty());

  std::shared_ptr<RegularSurface> result = links.front()->surfacePtr();
  for (auto it = std::next(links.begin()); it != links.end(); ++it) {
    assert(result != nullptr);
    result = mergedSurface(*result, it->get()->surface(), direction);
  }

  return result;
}
}  // namespace

CompositePortalLink::CompositePortalLink(std::unique_ptr<PortalLinkBase> a,
                                         std::unique_ptr<PortalLinkBase> b,
                                         AxisDirection direction, bool flatten)
    : PortalLinkBase(mergedSurface(a->surface(), b->surface(), direction)),
      m_direction{direction} {
  if (!flatten) {
    m_children.push_back(std::move(a));
    m_children.push_back(std::move(b));

  } else {
    auto handle = [&](std::unique_ptr<PortalLinkBase> link) {
      if (auto* composite = dynamic_cast<CompositePortalLink*>(link.get());
          composite != nullptr) {
        m_children.insert(
            m_children.end(),
            std::make_move_iterator(composite->m_children.begin()),
            std::make_move_iterator(composite->m_children.end()));
      } else {
        m_children.push_back(std::move(link));
      }
    };

    handle(std::move(a));
    handle(std::move(b));
  }
}

CompositePortalLink::CompositePortalLink(
    std::vector<std::unique_ptr<PortalLinkBase>> links, AxisDirection direction,
    bool flatten)
    : PortalLinkBase(mergePortalLinks(links, direction)),
      m_direction(direction) {
  if (!flatten) {
    for (auto& child : links) {
      m_children.push_back(std::move(child));
    }

  } else {
    auto handle = [&](std::unique_ptr<PortalLinkBase> link) {
      if (auto* composite = dynamic_cast<CompositePortalLink*>(link.get());
          composite != nullptr) {
        m_children.insert(
            m_children.end(),
            std::make_move_iterator(composite->m_children.begin()),
            std::make_move_iterator(composite->m_children.end()));
      } else {
        m_children.push_back(std::move(link));
      }
    };

    for (auto& child : links) {
      handle(std::move(child));
    }
  }
}

Result<const TrackingVolume*> CompositePortalLink::resolveVolume(
    const GeometryContext& gctx, const Vector2& position,
    double tolerance) const {
  // In this overload, we have to go back to global, because the children have
  // their own local coordinate systems
  Vector3 global = m_surface->localToGlobal(gctx, position);
  auto res = resolveVolume(gctx, global, tolerance);
  if (!res.ok()) {
    return Result<const TrackingVolume*>::failure(res.error());
  }
  return Result<const TrackingVolume*>::success(*res);
}

Result<const TrackingVolume*> CompositePortalLink::resolveVolume(
    const GeometryContext& gctx, const Vector3& position,
    double tolerance) const {
  assert(m_surface->isOnSurface(gctx, position, BoundaryTolerance::None(),
                                tolerance));

  for (const auto& child : m_children) {
    if (child->surface().isOnSurface(gctx, position, BoundaryTolerance::None(),
                                     tolerance)) {
      return child->resolveVolume(gctx, position);
    }
  }

  return Result<const TrackingVolume*>::failure(
      PortalError::PositionNotOnAnyChildPortalLink);
}

std::size_t CompositePortalLink::depth() const {
  std::size_t result = 1;

  for (const auto& child : m_children) {
    if (const auto* composite =
            dynamic_cast<const CompositePortalLink*>(child.get())) {
      result = std::max(result, composite->depth() + 1);
    }
  }

  return result;
}

std::size_t CompositePortalLink::size() const {
  return m_children.size();
}

void CompositePortalLink::toStream(std::ostream& os) const {
  os << "CompositePortalLink";
}

std::unique_ptr<GridPortalLink> CompositePortalLink::makeGrid(
    const GeometryContext& gctx, const Logger& logger) const {
  ACTS_VERBOSE("Attempting to make a grid from a composite portal link");
  std::vector<TrivialPortalLink*> trivialLinks;
  std::ranges::transform(m_children, std::back_inserter(trivialLinks),
                         [](const auto& child) {
                           return dynamic_cast<TrivialPortalLink*>(child.get());
                         });

  if (std::ranges::any_of(trivialLinks,
                          [](auto link) { return link == nullptr; })) {
    ACTS_VERBOSE(
        "Failed to make a grid from a composite portal link -> returning "
        "nullptr");
    return nullptr;
  }

  // we already know all children have surfaces that are compatible, we produced
  // a merged surface for the overall dimensions.

  auto printEdges = [](const auto& edges) -> std::string {
    std::vector<std::string> strings;
    std::ranges::transform(edges, std::back_inserter(strings),
                           [](const auto& v) { return std::to_string(v); });
    return boost::algorithm::join(strings, ", ");
  };

  if (surface().type() == Surface::SurfaceType::Cylinder) {
    ACTS_VERBOSE("Combining composite into cylinder grid");

    if (m_direction != AxisDirection::AxisZ) {
      ACTS_ERROR("Cylinder grid only supports binning in Z direction");
      throw std::runtime_error{"Unsupported binning direction"};
    }

    std::vector<double> edges;
    edges.reserve(m_children.size() + 1);

    const Transform3& groupTransform = m_surface->transform(gctx);
    Transform3 itransform = groupTransform.inverse();

    std::ranges::sort(
        trivialLinks, [&itransform, &gctx](const auto& a, const auto& b) {
          return (itransform * a->surface().transform(gctx)).translation()[eZ] <
                 (itransform * b->surface().transform(gctx)).translation()[eZ];
        });

    for (const auto& [i, child] : enumerate(trivialLinks)) {
      const auto& bounds =
          dynamic_cast<const CylinderBounds&>(child->surface().bounds());
      Transform3 ltransform = itransform * child->surface().transform(gctx);
      double hlZ = bounds.get(CylinderBounds::eHalfLengthZ);
      double minZ = ltransform.translation()[eZ] - hlZ;
      double maxZ = ltransform.translation()[eZ] + hlZ;
      if (i == 0) {
        edges.push_back(minZ);
      }
      edges.push_back(maxZ);
    }

    ACTS_VERBOSE("~> Determined bin edges to be " << printEdges(edges));

    Axis axis{AxisBound, edges};

    auto grid = GridPortalLink::make(m_surface, m_direction, std::move(axis));
    for (const auto& [i, child] : enumerate(trivialLinks)) {
      grid->atLocalBins({i + 1}) = &child->volume();
    }

    return grid;

  } else if (surface().type() == Surface::SurfaceType::Disc) {
    ACTS_VERBOSE("Combining composite into disc grid");

    if (m_direction != AxisDirection::AxisR) {
      ACTS_ERROR("Disc grid only supports binning in R direction");
      throw std::runtime_error{"Unsupported binning direction"};
    }

    std::vector<double> edges;
    edges.reserve(m_children.size() + 1);

    std::ranges::sort(trivialLinks, [](const auto& a, const auto& b) {
      const auto& boundsA =
          dynamic_cast<const RadialBounds&>(a->surface().bounds());
      const auto& boundsB =
          dynamic_cast<const RadialBounds&>(b->surface().bounds());
      return boundsA.get(RadialBounds::eMinR) <
             boundsB.get(RadialBounds::eMinR);
    });

    for (const auto& [i, child] : enumerate(trivialLinks)) {
      const auto& bounds =
          dynamic_cast<const RadialBounds&>(child->surface().bounds());

      if (i == 0) {
        edges.push_back(bounds.get(RadialBounds::eMinR));
      }
      edges.push_back(bounds.get(RadialBounds::eMaxR));
    }

    ACTS_VERBOSE("~> Determined bin edges to be " << printEdges(edges));

    Axis axis{AxisBound, edges};

    auto grid = GridPortalLink::make(m_surface, m_direction, std::move(axis));
    for (const auto& [i, child] : enumerate(trivialLinks)) {
      grid->atLocalBins({i + 1}) = &child->volume();
    }

    return grid;
  } else if (surface().type() == Surface::SurfaceType::Plane) {
    throw std::runtime_error{"Plane surfaces not implemented yet"};
  } else {
    throw std::invalid_argument{"Unsupported surface type"};
  }
}

}  // namespace Acts
