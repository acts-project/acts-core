// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceConcept.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Result.hpp"

#include <memory>
#include <string>

namespace Acts {

class DetectorElementBase;
class PlanarBounds;
class SurfaceBounds;

/// @class PlaneSurface
///
/// Class for a planaer in the TrackingGeometry.
///
/// The PlaneSurface extends the Surface class with the possibility to
/// convert local to global positions (vice versa).
///
/// @image html PlaneSurface.png
///
class PlaneSurface : public RegularSurface {
  friend class Surface;

 protected:
  /// Copy Constructor
  ///
  /// @param other is the source surface for the copy
  PlaneSurface(const PlaneSurface& other);

  /// Copy constructor - with shift
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param other is the source cone surface
  /// @param transform is the additional transform applied after copying
  PlaneSurface(const GeometryContext& gctx, const PlaneSurface& other,
               const Transform3& transform);

  /// Constructor from DetectorElementBase : Element proxy
  ///
  /// @param pbounds are the provided planar bounds
  /// @param detelement is the linked detector element to this surface
  PlaneSurface(std::shared_ptr<const PlanarBounds> pbounds,
               const DetectorElementBase& detelement);

  /// Constructor for Planes with (optional) shared bounds object
  ///
  /// @param transform transform in 3D that positions this surface
  /// @param pbounds bounds object to describe the actual surface area
  PlaneSurface(const Transform3& transform,
               std::shared_ptr<const PlanarBounds> pbounds = nullptr);

 public:
  /// Assignment operator
  ///
  /// @param other The source PlaneSurface for assignment
  PlaneSurface& operator=(const PlaneSurface& other);

  // Use overloads from `RegularSurface`
  using RegularSurface::globalToLocal;
  using RegularSurface::localToGlobal;
  using RegularSurface::normal;

  /// Get the normal vector of this surface at a given local position
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition is the local position is ignored
  ///
  /// return a Vector3 by value
  Vector3 normal(const GeometryContext& gctx,
                 const Vector2& lposition) const final;

  /// Get the normal vector of this surface at a given global position
  /// @note The @p position is required to be on-surface.
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the global positiono (for @ref PlaneSurface this is ignored)
  /// @return The normal vector
  Vector3 normal(const GeometryContext& gctx,
                 const Vector3& position) const final;

  /// Get the normal vector, independent of the location
  /// @param gctx The current geometry context object, e.g. alignment
  /// @return The normal vector
  Vector3 normal(const GeometryContext& gctx) const;

  /// The axis position is the position calculated
  /// for a certain axis type
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param aDir is the axis direction of reference position request
  ///
  /// @return position that can be used for this axis
  Vector3 referencePosition(const GeometryContext& gctx,
                            AxisDirection aDir) const final;

  /// Return the surface type
  SurfaceType type() const override;

  /// Return method for bounds object of this surfrace
  const SurfaceBounds& bounds() const override;

  /// Local to global transformation
  ///
  /// @note For planar surfaces the momentum direction is ignored in the local to global
  /// transformation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition local 2D position in specialized surface frame
  ///
  /// @return the global position by value
  Vector3 localToGlobal(const GeometryContext& gctx,
                        const Vector2& lposition) const override;

  /// Global to local transformation
  ///
  /// @note For planar surfaces the momentum direction is ignored in the global to local
  /// transformation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param tolerance optional tolerance within which a point is considered
  /// valid on surface
  ///
  /// @return a Result<Vector2> which can be !ok() if the operation fails
  Result<Vector2> globalToLocal(
      const GeometryContext& gctx, const Vector3& position,
      double tolerance = s_onSurfaceTolerance) const override;

  /// Method that calculates the correction due to incident angle
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position global 3D position (ignored for @ref PlaneSurface)
  /// @param direction global 3D momentum direction (ignored for @ref PlaneSurface)
  /// @return a double representing the scaling factor
  double pathCorrection(const GeometryContext& gctx, const Vector3& position,
                        const Vector3& direction) const final;

  /// @brief Straight line intersection
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The start position of the intersection attempt
  /// @param direction The direction of the intersection attempt,
  /// (@note expected to be normalized)
  /// @param boundaryTolerance The boundary check directive
  /// @param tolerance the tolerance used for the intersection
  ///
  /// <b>mathematical motivation:</b>
  ///
  /// the equation of the plane is given by: <br>
  /// @f$ \vec n \cdot \vec x = \vec n \cdot \vec p,@f$ <br>
  /// where @f$ \vec n = (n_{x}, n_{y}, n_{z})@f$ denotes the normal vector of
  /// the plane,  @f$ \vec p = (p_{x}, p_{y}, p_{z})@f$ one specific point
  /// on the plane and @f$ \vec x = (x,y,z) @f$ all possible points
  /// on the plane.<br>
  ///
  /// Given a line with:<br>
  /// @f$ \vec l(u) = \vec l_{1} + u \cdot \vec v @f$, <br>
  /// the solution for @f$ u @f$ can be written:
  /// @f$ u = \frac{\vec n (\vec p - \vec l_{1})}{\vec n \vec v}@f$ <br>
  /// If the denominator is 0 then the line lies:
  /// - either in the plane
  /// - perpendicular to the normal of the plane
  ///
  /// @return the @c SurfaceMultiIntersection object
  SurfaceMultiIntersection intersect(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction,
      const BoundaryTolerance& boundaryTolerance =
          BoundaryTolerance::Infinite(),
      double tolerance = s_onSurfaceTolerance) const final;

  /// Return a Polyhedron for the surfaces
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param quarterSegments is the number of segments used to describe curved
  /// segments in a quarter of the phi range. If it is 1, then only the extrema
  /// points in phi are inserted next to the segment corners.
  ///
  /// @note for planar surfaces without curved segments @c quarterSegments is ignored
  ///
  /// @return A list of vertices and a face/facett description of it
  Polyhedron polyhedronRepresentation(
      const GeometryContext& gctx, unsigned int quarterSegments) const override;

  /// Return properly formatted class name for screen output
  std::string name() const override;

  /// Calculate the derivative of bound track parameters local position w.r.t.
  /// position in local 3D Cartesian coordinates
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The position of the parameters in global
  ///
  /// @return Derivative of bound local position w.r.t. position in local 3D
  /// cartesian coordinates
  ActsMatrix<2, 3> localCartesianToBoundLocalDerivative(
      const GeometryContext& gctx, const Vector3& position) const final;

  /// Merge two plane surfaces into a single one.
  /// @note The surfaces need to be *compatible*, i.e. have bounds
  ///       that align along merging direction, and have the same bound size
  ///       along the non-merging direction
  /// @param other The other plane surface to merge with
  /// @param direction The direction: either @c AxisX or @c AxisY
  /// @param logger The logger to use
  /// @return The merged plane surface and a boolean indicating if surfaces are reversed
  /// @note The returned boolean is `false` if `this` is *left* or
  ///       *counter-clockwise* of @p other, and `true` if not.
  std::pair<std::shared_ptr<PlaneSurface>, bool> mergedWith(
      const PlaneSurface& other, AxisDirection direction,
      const Logger& logger = getDummyLogger()) const;

 protected:
  /// the bounds of this surface
  std::shared_ptr<const PlanarBounds> m_bounds;

 private:
};

static_assert(RegularSurfaceConcept<PlaneSurface>,
              "PlaneSurface does not fulfill RegularSurfaceConcept");

}  // namespace Acts
