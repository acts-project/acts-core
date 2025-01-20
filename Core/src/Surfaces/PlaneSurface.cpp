// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PlaneSurface.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/SurfaceError.hpp"
#include "Acts/Surfaces/SurfaceMergingException.hpp"
#include "Acts/Surfaces/detail/FacesHelper.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Acts {

PlaneSurface::PlaneSurface(const PlaneSurface& other)
    : GeometryObject(), RegularSurface(other), m_bounds(other.m_bounds) {}

PlaneSurface::PlaneSurface(const GeometryContext& gctx,
                           const PlaneSurface& other,
                           const Transform3& transform)
    : GeometryObject(),
      RegularSurface(gctx, other, transform),
      m_bounds(other.m_bounds) {}

PlaneSurface::PlaneSurface(std::shared_ptr<const PlanarBounds> pbounds,
                           const DetectorElementBase& detelement)
    : RegularSurface(detelement), m_bounds(std::move(pbounds)) {
  // surfaces representing a detector element must have bounds
  throw_assert(m_bounds, "PlaneBounds must not be nullptr");
}

PlaneSurface::PlaneSurface(const Transform3& transform,
                           std::shared_ptr<const PlanarBounds> pbounds)
    : RegularSurface(transform), m_bounds(std::move(pbounds)) {}

PlaneSurface& PlaneSurface::operator=(const PlaneSurface& other) {
  if (this != &other) {
    Surface::operator=(other);
    m_bounds = other.m_bounds;
  }
  return *this;
}

Surface::SurfaceType PlaneSurface::type() const {
  return Surface::Plane;
}

Vector3 PlaneSurface::localToGlobal(const GeometryContext& gctx,
                                    const Vector2& lposition) const {
  return transform(gctx) * Vector3(lposition[0], lposition[1], 0.);
}

Result<Vector2> PlaneSurface::globalToLocal(const GeometryContext& gctx,
                                            const Vector3& position,
                                            double tolerance) const {
  Vector3 loc3Dframe = transform(gctx).inverse() * position;
  if (std::abs(loc3Dframe.z()) > std::abs(tolerance)) {
    return Result<Vector2>::failure(SurfaceError::GlobalPositionNotOnSurface);
  }
  return Result<Vector2>::success({loc3Dframe.x(), loc3Dframe.y()});
}

std::string PlaneSurface::name() const {
  return "Acts::PlaneSurface";
}

const SurfaceBounds& PlaneSurface::bounds() const {
  if (m_bounds) {
    return (*m_bounds.get());
  }
  return s_noBounds;
}

Polyhedron PlaneSurface::polyhedronRepresentation(
    const GeometryContext& gctx, unsigned int quarterSegments) const {
  // Prepare vertices and faces
  std::vector<Vector3> vertices;
  bool exactPolyhedron = true;

  // If you have bounds you can create a polyhedron representation
  if (m_bounds) {
    auto vertices2D = m_bounds->vertices(quarterSegments);
    vertices.reserve(vertices2D.size() + 1);
    for (const auto& v2D : vertices2D) {
      vertices.push_back(transform(gctx) * Vector3(v2D.x(), v2D.y(), 0.));
    }
    bool isEllipse = bounds().type() == SurfaceBounds::eEllipse;
    bool innerExists = false, coversFull = false;
    if (isEllipse) {
      exactPolyhedron = false;
      auto vStore = bounds().values();
      innerExists = vStore[EllipseBounds::eInnerRx] > s_epsilon &&
                    vStore[EllipseBounds::eInnerRy] > s_epsilon;
      coversFull = std::abs(vStore[EllipseBounds::eHalfPhiSector] -
                            std::numbers::pi) < s_epsilon;
    }
    // All of those can be described as convex
    // @todo same as for Discs: coversFull is not the right criterium
    // for triangulation
    if (!isEllipse || !innerExists || !coversFull) {
      auto [faces, triangularMesh] =
          detail::FacesHelper::convexFaceMesh(vertices);
      return Polyhedron(vertices, faces, triangularMesh, exactPolyhedron);
    } else {
      // Two concentric rings, we use the pure concentric method momentarily,
      // but that creates too  many unneccesarry faces, when only two
      // are needed to describe the mesh, @todo investigate merging flag
      auto [faces, triangularMesh] =
          detail::FacesHelper::cylindricalFaceMesh(vertices);
      return Polyhedron(vertices, faces, triangularMesh, exactPolyhedron);
    }
  }
  throw std::domain_error(
      "Polyhedron representation of boundless surface not possible.");
}

Vector3 PlaneSurface::normal(const GeometryContext& gctx,
                             const Vector2& /*lpos*/) const {
  return normal(gctx);
}

Vector3 PlaneSurface::normal(const GeometryContext& gctx,
                             const Vector3& /*pos*/) const {
  return normal(gctx);
}

Vector3 PlaneSurface::normal(const GeometryContext& gctx) const {
  return transform(gctx).linear().col(2);
}

Vector3 PlaneSurface::referencePosition(const GeometryContext& gctx,
                                        AxisDirection /*aDir*/) const {
  return center(gctx);
}

double PlaneSurface::pathCorrection(const GeometryContext& gctx,
                                    const Vector3& /*position*/,
                                    const Vector3& direction) const {
  // We can ignore the global position here
  return 1. / std::abs(normal(gctx).dot(direction));
}

SurfaceMultiIntersection PlaneSurface::intersect(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const BoundaryTolerance& boundaryTolerance,
    double tolerance) const {
  // Get the contextual transform
  const auto& gctxTransform = transform(gctx);
  // Use the intersection helper for planar surfaces
  auto intersection =
      PlanarHelper::intersect(gctxTransform, position, direction, tolerance);
  auto status = intersection.status();
  // Evaluate boundary check if requested (and reachable)
  if (intersection.status() != IntersectionStatus::unreachable) {
    // Built-in local to global for speed reasons
    const auto& tMatrix = gctxTransform.matrix();
    // Create the reference vector in local
    const Vector3 vecLocal(intersection.position() - tMatrix.block<3, 1>(0, 3));
    if (!insideBounds(tMatrix.block<3, 2>(0, 0).transpose() * vecLocal,
                      boundaryTolerance)) {
      status = IntersectionStatus::unreachable;
    }
  }
  return {{Intersection3D(intersection.position(), intersection.pathLength(),
                          status),
           Intersection3D::invalid()},
          this};
}

ActsMatrix<2, 3> PlaneSurface::localCartesianToBoundLocalDerivative(
    const GeometryContext& /*gctx*/, const Vector3& /*position*/) const {
  const ActsMatrix<2, 3> loc3DToLocBound = ActsMatrix<2, 3>::Identity();
  return loc3DToLocBound;
}

std::pair<std::shared_ptr<PlaneSurface>, bool> PlaneSurface::mergedWith(
    const PlaneSurface& other, AxisDirection direction,
    const Logger& logger) const {
  ACTS_VERBOSE("Merging plane surfaces in " << axisDirectionName(direction)
                                            << " direction");

  if (m_associatedDetElement != nullptr ||
      other.m_associatedDetElement != nullptr) {
    throw SurfaceMergingException(getSharedPtr(), other.getSharedPtr(),
                                  "PlaneSurface::merge: surfaces are "
                                  "associated with a detector element");
  }

  assert(m_transform != nullptr && other.m_transform != nullptr);

  Transform3 otherLocal = m_transform->inverse() * *other.m_transform;

  // TODO: Is it a good tolerance?
  constexpr auto tolerance = s_onSurfaceTolerance;

  // Surface cannot have any relative rotation
  if ((otherLocal.rotation().matrix() - RotationMatrix3::Identity()).norm() >
      tolerance) {
    ACTS_ERROR("PlaneSurface::merge: surfaces have relative rotation");
    throw SurfaceMergingException(
        getSharedPtr(), other.getSharedPtr(),
        "PlaneSurface::merge: surfaces have relative rotation");
  }

  const RectangleBounds* thisBounds =
      dynamic_cast<const RectangleBounds*>(&bounds());
  const RectangleBounds* otherBounds =
      dynamic_cast<const RectangleBounds*>(&other.bounds());

  if (thisBounds == nullptr || otherBounds == nullptr) {
    throw SurfaceMergingException(
        getSharedPtr(), other.getSharedPtr(),
        "PlaneSurface::merge: only Rectangle Bounds are supported");
  }

  if (direction != AxisDirection::AxisX && direction != AxisDirection::AxisY) {
    throw SurfaceMergingException(getSharedPtr(), other.getSharedPtr(),
                                  "PlaneSurface::merge: invalid direction " +
                                      axisDirectionName(direction));
  }

  bool mergeX = direction == AxisDirection::AxisX;

  double thisHalfMerge =
      mergeX ? thisBounds->halfLengthX() : thisBounds->halfLengthY();
  double otherHalfMerge =
      mergeX ? otherBounds->halfLengthX() : otherBounds->halfLengthY();

  double thisHalfNonMerge =
      mergeX ? thisBounds->halfLengthY() : thisBounds->halfLengthX();
  double otherHalfNonMerge =
      mergeX ? otherBounds->halfLengthY() : otherBounds->halfLengthX();

  if (std::abs(thisHalfNonMerge - otherHalfNonMerge) > tolerance) {
    ACTS_ERROR(
        "PlaneSurface::merge: surfaces have different non-merging lengths");
    throw SurfaceMergingException(
        getSharedPtr(), other.getSharedPtr(),
        "PlaneSurface::merge: surfaces have different non-merging lengths");
  }
  Vector3 otherTranslation = otherLocal.translation();

  // No translation in non-merging direction/z is allowed
  double nonMergeShift = mergeX ? otherTranslation.y() : otherTranslation.x();

  if (std::abs(nonMergeShift) > tolerance ||
      std::abs(otherTranslation.z()) > tolerance) {
    ACTS_ERROR(
        "PlaneSurface::merge: surfaces have relative translation in y/z");
    throw SurfaceMergingException(
        getSharedPtr(), other.getSharedPtr(),
        "PlaneSurface::merge: surfaces have relative translation in y/z");
  }

  double mergeShift = mergeX ? otherTranslation.x() : otherTranslation.y();

  double thisMinMerge = -thisHalfMerge;
  double thisMaxMerge = thisHalfMerge;

  double otherMinMerge = mergeShift - otherHalfMerge;
  double otherMaxMerge = mergeShift + otherHalfMerge;

  // Surfaces have to "touch" along merging direction
  if (std::abs(thisMaxMerge - otherMinMerge) > tolerance &&
      std::abs(thisMinMerge - otherMaxMerge) > tolerance) {
    ACTS_ERROR(
        "PlaneSurface::merge: surfaces have incompatible merge bound location");
    throw SurfaceMergingException(
        getSharedPtr(), other.getSharedPtr(),
        "PlaneSurface::merge: surfaces have incompatible merge bound location");
  }

  double newMaxMerge = std::max(thisMaxMerge, otherMaxMerge);
  double newMinMerge = std::min(thisMinMerge, otherMinMerge);

  double newHalfMerge = (newMaxMerge - newMinMerge) / 2;

  double newMidMerge = (newMaxMerge + newMinMerge) / 2;

  auto newBounds =
      mergeX
          ? std::make_shared<RectangleBounds>(newHalfMerge, thisHalfNonMerge)
          : std::make_shared<RectangleBounds>(thisHalfNonMerge, newHalfMerge);

  Vector3 unitDir = mergeX ? Vector3::UnitX() : Vector3::UnitY();
  Transform3 newTransform = *m_transform * Translation3{unitDir * newMidMerge};
  return {Surface::makeShared<PlaneSurface>(newTransform, newBounds),
          mergeShift < 0};
}

}  // namespace Acts
