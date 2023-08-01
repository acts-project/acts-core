// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

namespace Acts {

class RegularSurface : public Surface {
 public:
  // Reuse all constructors from the base class
  using Surface::Surface;

  /// Calculate the normal vector of the surface
  /// This overload requires an on-surface local position
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition is the local position where the normal vector is
  /// constructed
  ///
  /// @return normal vector by value
  virtual Vector3 normal(const GeometryContext& gctx,
                         const Vector2& lposition) const = 0;

  /// Calculate the normal vector of the surface
  /// This overload accepts a global position
  ///
  /// @param position is the global position where the normal vector is
  /// constructed
  /// @note The @p position is required to be on-surface.
  /// 			Use @ref Surface::coerceToSurface(const GeometryContext&, const Vector3&, const Vector3&)
  /// @param gctx The current geometry context object, e.g. alignment
  /// @return normal vector by value
  virtual Vector3 normal(const GeometryContext& gctx,
                         const Vector3& position) const;

  /// Calculate the normal vector of the surface
  /// This overload is fully generic, fulfills the @ref Surface interface and
  /// accepts a global position and a direction. For @c RegularSurface this is
  /// equivalent to the @ref normal(const GeometryContext&, const Vector3&)
  /// overload, ignoring the @p direction
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param pos is the global position where the normal vector is constructed
  /// @param direction is the direction of the normal vector (ignored for @c RegularSurface)
  Vector3 normal(const GeometryContext& gctx, const Vector3& pos,
                 const Vector3& direction) const override;
};
}  // namespace Acts
