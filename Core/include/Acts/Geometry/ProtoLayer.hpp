// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

#include <iostream>
#include <memory>
#include <utility>
#include <vector>

namespace Acts {

namespace detail {

/// @class ProtoLayerBase
///
/// Base class containing common functionality for ProtoLayer implementations
/// @note This will go away once we remove the Gen1 geometry which assumes this only takes const pointers.
struct ProtoLayerBase {
 public:
  /// The extent of the ProtoLayer
  Extent extent;

  /// The envelope parameters
  ExtentEnvelope envelope = ExtentEnvelope::Zero();

  /// The local transform
  Transform3 transform = Transform3::Identity();

  /// Get the parameters : min
  /// @param aDir The accessed axis direction
  /// @param addenv The steering if enevlope is added or not
  double min(AxisDirection aDir, bool addenv = true) const;

  // Get the  parameters : max
  /// @param aDir The accessed axis direction
  /// @param addenv The steering if enevlope is added or not
  double max(AxisDirection aDir, bool addenv = true) const;

  // Get the  parameters : medium
  /// @param aDir The accessed axis direction
  /// @param addenv The steering if enevlope is added or not
  double medium(AxisDirection aDir, bool addenv = true) const;

  // Get the  parameters : range
  /// @param aDir The accessed axis direction
  /// @param addenv The steering if enevlope is added or not
  double range(AxisDirection aDir, bool addenv = true) const;

  /// Output to ostream
  /// @param sl the input ostream
  std::ostream& toStream(std::ostream& sl) const;

 protected:
  /// Helper method which performs the actual min/max calculation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param surfaces The surfaces to build this protolayer out of
  /// @param extent The extent to modify
  /// @param transform The transform to use
  static void measureImpl(const GeometryContext& gctx,
                          const std::vector<const Surface*>& surfaces,
                          Extent& extent, const Transform3& transform);
};

/// @struct ProtoLayerT
///
/// Encapsulates min/max boundaries that will be turned into a layer.
/// The struct allows this information to be obtained in a consistent
/// way, or be caller provided.
template <bool IsConst>
struct ProtoLayerT : public ProtoLayerBase {
  using SurfacePtr = std::conditional_t<IsConst, const Surface*, Surface*>;
  using SurfaceType = std::conditional_t<IsConst, const Surface, Surface>;

  /// Constructor
  ///
  /// Loops over a provided vector of surface and calculates the various
  /// min/max values in one go. Also takes into account the thickness
  /// of an associated DetectorElement, if it exists.
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param surfaces The vector of surfaces to consider
  /// @param transformIn The local transform to evaluate the sizing in
  ProtoLayerT(const GeometryContext& gctx,
              const std::vector<SurfacePtr>& surfaces,
              const Transform3& transformIn = Transform3::Identity())
      : m_surfaces(surfaces) {
    transform = transformIn;
    std::vector<const Surface*> constSurfaces;
    if constexpr (!IsConst) {
      constSurfaces.reserve(surfaces.size());
      for (auto* sf : surfaces) {
        constSurfaces.push_back(sf);
      }
      measureImpl(gctx, constSurfaces, extent, transform);
    } else {
      measureImpl(gctx, surfaces, extent, transform);
    }
  }

  /// Constructor
  ///
  /// Loops over a provided vector of surface and calculates the various
  /// min/max values in one go. Also takes into account the thickness
  /// of an associated DetectorElement, if it exists.
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param surfaces The vector of surfaces to consider
  /// @param transformIn The local transform to evaluate the sizing in
  ProtoLayerT(const GeometryContext& gctx,
              const std::vector<std::shared_ptr<SurfaceType>>& surfaces,
              const Transform3& transformIn = Transform3::Identity()) {
    transform = transformIn;
    m_surfaces.reserve(surfaces.size());
    for (const auto& sf : surfaces) {
      m_surfaces.push_back(sf.get());
    }
    std::vector<const Surface*> constSurfaces;
    if constexpr (!IsConst) {
      constSurfaces.reserve(surfaces.size());
      for (auto* sf : m_surfaces) {
        constSurfaces.push_back(sf);
      }
      measureImpl(gctx, constSurfaces, extent, transform);
    } else {
      measureImpl(gctx, m_surfaces, extent, transform);
    }
  }

  /// Constructor that accepts non-const shared pointers even when IsConst is
  /// true
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param surfaces The vector of surfaces to consider
  /// @param transformIn The local transform to evaluate the sizing in
  ProtoLayerT(const GeometryContext& gctx,
              const std::vector<std::shared_ptr<Surface>>& surfaces,
              const Transform3& transformIn = Transform3::Identity())
    requires(IsConst)
  {
    transform = transformIn;
    m_surfaces.reserve(surfaces.size());
    for (const auto& sf : surfaces) {
      m_surfaces.push_back(sf.get());
    }
    measureImpl(gctx, m_surfaces, extent, transform);
  }

  ProtoLayerT() = default;

  /// Output stream operator
  /// @param sl the input ostream
  /// @param pl the ProtoLayer to be printed
  /// @return the output ostream
  friend std::ostream& operator<<(std::ostream& sl, const ProtoLayerT& pl) {
    return pl.toStream(sl);
  }

  /// Give access to the surfaces used/assigned to the ProtoLayer
  const std::vector<SurfacePtr>& surfaces() const { return m_surfaces; }

  /// Add a surface, this will also increase the extent
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param surface The surface which is added to the ProtoLayer
  void add(const GeometryContext& gctx, SurfaceType& surface) {
    m_surfaces.push_back(&surface);
    std::vector<const Surface*> constSurfaces;
    if constexpr (!IsConst) {
      constSurfaces.reserve(m_surfaces.size());
      for (auto* sf : m_surfaces) {
        constSurfaces.push_back(sf);
      }
      measureImpl(gctx, constSurfaces, extent, transform);
    } else {
      measureImpl(gctx, m_surfaces, extent, transform);
    }
  }

 private:
  /// Store the list of surfaces used for this proto layer
  std::vector<SurfacePtr> m_surfaces = {};
};

}  // namespace detail

// Forward-declaration friendly class for backward compatibility
struct ProtoLayer : public detail::ProtoLayerT<true> {
  using detail::ProtoLayerT<true>::ProtoLayerT;
};

struct MutableProtoLayer : public detail::ProtoLayerT<false> {
  using detail::ProtoLayerT<false>::ProtoLayerT;
};

}  // namespace Acts
