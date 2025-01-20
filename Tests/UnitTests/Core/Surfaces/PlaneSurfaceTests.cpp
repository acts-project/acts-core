// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/SurfaceMergingException.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

#include <cmath>
#include <memory>
#include <numbers>
#include <string>

using namespace Acts::UnitLiterals;

namespace Acts::Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

BOOST_AUTO_TEST_SUITE(PlaneSurfaces)
/// Unit test for creating compliant/non-compliant PlaneSurface object
BOOST_AUTO_TEST_CASE(PlaneSurfaceConstruction) {
  /// Test default construction
  // default construction is deleted

  // bounds object, rectangle type
  auto rBounds = std::make_shared<const RectangleBounds>(3., 4.);
  /// Constructor with transform and bounds
  Translation3 translation{0., 1., 2.};
  auto pTransform = Transform3(translation);

  /// Constructor with transform
  BOOST_CHECK_EQUAL(
      Surface::makeShared<PlaneSurface>(pTransform, rBounds)->type(),
      Surface::Plane);

  /// Copy constructor
  auto planeSurfaceObject =
      Surface::makeShared<PlaneSurface>(pTransform, rBounds);
  auto copiedPlaneSurface =
      Surface::makeShared<PlaneSurface>(*planeSurfaceObject);
  BOOST_CHECK_EQUAL(copiedPlaneSurface->type(), Surface::Plane);
  BOOST_CHECK(*copiedPlaneSurface == *planeSurfaceObject);

  /// Copied and transformed
  auto copiedTransformedPlaneSurface = Surface::makeShared<PlaneSurface>(
      tgContext, *planeSurfaceObject, pTransform);
  BOOST_CHECK_EQUAL(copiedTransformedPlaneSurface->type(), Surface::Plane);

  /// Construct with nullptr bounds
  DetectorElementStub detElem;
  BOOST_CHECK_THROW(
      auto nullBounds = Surface::makeShared<PlaneSurface>(nullptr, detElem),
      AssertionFailureException);
}

/// Unit test for testing PlaneSurface properties
BOOST_AUTO_TEST_CASE(PlaneSurfaceProperties) {
  // bounds object, rectangle type
  auto rBounds = std::make_shared<const RectangleBounds>(3., 4.);

  /// Test clone method
  Translation3 translation{0., 1., 2.};
  auto pTransform = Transform3(translation);
  auto planeSurfaceObject =
      Surface::makeShared<PlaneSurface>(pTransform, rBounds);
  // Is it in the right place?
  Translation3 translation2{0., 2., 4.};
  auto pTransform2 = Transform3(translation2);
  auto planeSurfaceObject2 =
      Surface::makeShared<PlaneSurface>(pTransform2, rBounds);

  /// Test type (redundant)
  BOOST_CHECK_EQUAL(planeSurfaceObject->type(), Surface::Plane);

  /// Test referencePosition
  Vector3 referencePosition{0., 1., 2.};
  BOOST_CHECK_EQUAL(
      planeSurfaceObject->referencePosition(tgContext, AxisDirection::AxisX),
      referencePosition);

  /// Test referenceFrame
  Vector3 arbitraryGlobalPosition{2., 2., 2.};
  Vector3 momentum{1.e6, 1.e6, 1.e6};
  RotationMatrix3 expectedFrame;
  expectedFrame << 1., 0., 0., 0., 1., 0., 0., 0., 1.;

  CHECK_CLOSE_OR_SMALL(planeSurfaceObject->referenceFrame(
                           tgContext, arbitraryGlobalPosition, momentum),
                       expectedFrame, 1e-6, 1e-9);

  /// Test normal, given 3D position
  Vector3 normal3D(0., 0., 1.);
  BOOST_CHECK_EQUAL(planeSurfaceObject->normal(tgContext), normal3D);

  /// Test bounds
  BOOST_CHECK_EQUAL(planeSurfaceObject->bounds().type(),
                    SurfaceBounds::eRectangle);

  /// Test localToGlobal
  Vector2 localPosition{1.5, 1.7};
  Vector3 globalPosition =
      planeSurfaceObject->localToGlobal(tgContext, localPosition, momentum);
  // expected position is the translated one
  Vector3 expectedPosition{1.5 + translation.x(), 1.7 + translation.y(),
                           translation.z()};

  CHECK_CLOSE_REL(globalPosition, expectedPosition, 1e-2);

  /// Testing globalToLocal
  localPosition =
      planeSurfaceObject->globalToLocal(tgContext, globalPosition, momentum)
          .value();
  Vector2 expectedLocalPosition{1.5, 1.7};

  CHECK_CLOSE_REL(localPosition, expectedLocalPosition, 1e-2);

  Vector3 globalPositionOff =
      globalPosition +
      planeSurfaceObject->normal(tgContext, localPosition) * 0.1;

  BOOST_CHECK(
      planeSurfaceObject->globalToLocal(tgContext, globalPositionOff, momentum)
          .error());
  BOOST_CHECK(planeSurfaceObject
                  ->globalToLocal(tgContext, globalPositionOff, momentum, 0.05)
                  .error());
  BOOST_CHECK(planeSurfaceObject
                  ->globalToLocal(tgContext, globalPositionOff, momentum, 0.2)
                  .ok());

  /// Test isOnSurface
  Vector3 offSurface{0, 1, -2.};
  BOOST_CHECK(planeSurfaceObject->isOnSurface(
      tgContext, globalPosition, momentum, BoundaryTolerance::None()));
  BOOST_CHECK(planeSurfaceObject->isOnSurface(tgContext, globalPosition,
                                              BoundaryTolerance::None()));
  BOOST_CHECK(!planeSurfaceObject->isOnSurface(tgContext, offSurface, momentum,
                                               BoundaryTolerance::None()));
  BOOST_CHECK(!planeSurfaceObject->isOnSurface(tgContext, offSurface,
                                               BoundaryTolerance::None()));

  /// Test intersection
  Vector3 direction{0., 0., 1.};
  auto sfIntersection = planeSurfaceObject
                            ->intersect(tgContext, offSurface, direction,
                                        BoundaryTolerance::None())
                            .closest();
  Intersection3D expectedIntersect{Vector3{0, 1, 2}, 4.,
                                   IntersectionStatus::reachable};
  BOOST_CHECK(sfIntersection.isValid());
  BOOST_CHECK_EQUAL(sfIntersection.position(), expectedIntersect.position());
  BOOST_CHECK_EQUAL(sfIntersection.pathLength(),
                    expectedIntersect.pathLength());
  BOOST_CHECK_EQUAL(sfIntersection.object(), planeSurfaceObject.get());

  /// Test pathCorrection
  CHECK_CLOSE_REL(planeSurfaceObject->pathCorrection(tgContext, offSurface,
                                                     momentum.normalized()),
                  std::numbers::sqrt3, 0.01);

  /// Test name
  BOOST_CHECK_EQUAL(planeSurfaceObject->name(),
                    std::string("Acts::PlaneSurface"));

  /// Test dump
  boost::test_tools::output_test_stream dumpOutput;
  dumpOutput << planeSurfaceObject->toStream(tgContext);
  BOOST_CHECK(dumpOutput.is_equal(
      "Acts::PlaneSurface\n"
      "     Center position  (x, y, z) = (0.0000, 1.0000, 2.0000)\n"
      "     Rotation:             colX = (1.000000, 0.000000, 0.000000)\n"
      "                           colY = (0.000000, 1.000000, 0.000000)\n"
      "                           colZ = (0.000000, 0.000000, 1.000000)\n"
      "     Bounds  : Acts::RectangleBounds:  (hlX, hlY) = (3.0000000, "
      "4.0000000)\n"
      "(lower left, upper right):\n"
      "-3.0000000 -4.0000000\n"
      "3.0000000 4.0000000"));
}

BOOST_AUTO_TEST_CASE(PlaneSurfaceEqualityOperators) {
  // rectangle bounds
  auto rBounds = std::make_shared<const RectangleBounds>(3., 4.);
  Translation3 translation{0., 1., 2.};
  auto pTransform = Transform3(translation);
  auto planeSurfaceObject =
      Surface::makeShared<PlaneSurface>(pTransform, rBounds);
  auto planeSurfaceObject2 =
      Surface::makeShared<PlaneSurface>(pTransform, rBounds);

  /// Test equality operator
  BOOST_CHECK(*planeSurfaceObject == *planeSurfaceObject2);

  BOOST_TEST_CHECKPOINT(
      "Create and then assign a PlaneSurface object to the existing one");

  /// Test assignment
  auto assignedPlaneSurface =
      Surface::makeShared<PlaneSurface>(Transform3::Identity(), nullptr);
  *assignedPlaneSurface = *planeSurfaceObject;

  /// Test equality of assigned to original
  BOOST_CHECK(*assignedPlaneSurface == *planeSurfaceObject);
}

/// Unit test for testing PlaneSurface extent via Polyhedron representation
BOOST_AUTO_TEST_CASE(PlaneSurfaceExtent) {
  // First test - non-rotated
  static const Transform3 planeZX =
      AngleAxis3(-std::numbers::pi / 2., Vector3::UnitX()) *
      AngleAxis3(-std::numbers::pi / 2., Vector3::UnitZ()) *
      Transform3::Identity();

  double rHx = 2.;
  double rHy = 4.;
  double yPs = 3.;
  auto rBounds = std::make_shared<RectangleBounds>(rHx, rHy);

  auto plane = Surface::makeShared<PlaneSurface>(
      Transform3(Translation3(Vector3(0., yPs, 0.)) * planeZX), rBounds);

  auto planeExtent = plane->polyhedronRepresentation(tgContext, 1).extent();

  CHECK_CLOSE_ABS(planeExtent.min(AxisDirection::AxisZ), -rHx,
                  s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtent.max(AxisDirection::AxisZ), rHx,
                  s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtent.min(AxisDirection::AxisX), -rHy,
                  s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtent.max(AxisDirection::AxisX), rHy,
                  s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtent.min(AxisDirection::AxisY), yPs,
                  s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtent.max(AxisDirection::AxisY), yPs,
                  s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtent.min(AxisDirection::AxisR), yPs,
                  s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtent.max(AxisDirection::AxisR), std::hypot(yPs, rHy),
                  s_onSurfaceTolerance);

  // Now rotate
  double alpha = 0.123;
  auto planeRot = Surface::makeShared<PlaneSurface>(
      Transform3(Translation3(Vector3(0., yPs, 0.)) *
                 AngleAxis3(alpha, Vector3(0., 0., 1.)) * planeZX),
      rBounds);

  auto planeExtentRot =
      planeRot->polyhedronRepresentation(tgContext, 1).extent();
  CHECK_CLOSE_ABS(planeExtentRot.min(AxisDirection::AxisZ), -rHx,
                  s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtentRot.max(AxisDirection::AxisZ), rHx,
                  s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtentRot.min(AxisDirection::AxisX),
                  -rHy * std::cos(alpha), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtentRot.max(AxisDirection::AxisX),
                  rHy * std::cos(alpha), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtentRot.min(AxisDirection::AxisY),
                  yPs - rHy * std::sin(alpha), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtentRot.max(AxisDirection::AxisY),
                  yPs + rHy * std::sin(alpha), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtentRot.min(AxisDirection::AxisR),
                  yPs * std::cos(alpha), s_onSurfaceTolerance);
}

BOOST_AUTO_TEST_CASE(RotatedTrapezoid) {
  const double shortHalfX = 100.;
  const double longHalfX = 200.;
  const double halfY = 300.;
  const double rotAngle = 45._degree;

  Vector2 edgePoint{longHalfX - 10., halfY};

  std::shared_ptr<TrapezoidBounds> bounds =
      std::make_shared<TrapezoidBounds>(shortHalfX, longHalfX, halfY);

  BOOST_CHECK(bounds->inside(edgePoint, BoundaryTolerance::None()));
  BOOST_CHECK(!bounds->inside(Eigen::Rotation2D(-rotAngle) * edgePoint,
                              BoundaryTolerance::None()));

  std::shared_ptr<TrapezoidBounds> rotatedBounds =
      std::make_shared<TrapezoidBounds>(shortHalfX, longHalfX, halfY, rotAngle);

  BOOST_CHECK(!rotatedBounds->inside(edgePoint, BoundaryTolerance::None()));
  BOOST_CHECK(rotatedBounds->inside(Eigen::Rotation2D(-rotAngle) * edgePoint,
                                    BoundaryTolerance::None()));
}

/// Unit test for testing PlaneSurface alignment derivatives
BOOST_AUTO_TEST_CASE(PlaneSurfaceAlignment) {
  // bounds object, rectangle type
  auto rBounds = std::make_shared<const RectangleBounds>(3., 4.);
  // Test clone method
  Translation3 translation{0., 1., 2.};
  const double rotationAngle = std::numbers::pi / 2.;
  AngleAxis3 rotation(rotationAngle, Vector3::UnitY());
  RotationMatrix3 rotationMat = rotation.toRotationMatrix();

  auto pTransform = Transform3{translation * rotationMat};
  auto planeSurfaceObject =
      Surface::makeShared<PlaneSurface>(pTransform, rBounds);

  // The local frame z axis
  const Vector3 localZAxis = rotationMat.col(2);
  // Check the local z axis is aligned to global x axis
  CHECK_CLOSE_ABS(localZAxis, Vector3(1., 0., 0.), 1e-15);

  // Define the track (local) position and direction
  Vector2 localPosition{1, 2};
  Vector3 momentum{1, 0, 0};
  Vector3 direction = momentum.normalized();
  // Get the global position
  Vector3 globalPosition =
      planeSurfaceObject->localToGlobal(tgContext, localPosition, momentum);

  // (a) Test the derivative of path length w.r.t. alignment parameters
  const AlignmentToPathMatrix& alignToPath =
      planeSurfaceObject->alignmentToPathDerivative(tgContext, globalPosition,
                                                    direction);
  // The expected results
  AlignmentToPathMatrix expAlignToPath = AlignmentToPathMatrix::Zero();
  expAlignToPath << 1, 0, 0, 2, -1, 0;

  // Check if the calculated derivative is as expected
  CHECK_CLOSE_ABS(alignToPath, expAlignToPath, 1e-10);

  // (b) Test the derivative of bound track parameters local position w.r.t.
  // position in local 3D Cartesian coordinates
  const auto& loc3DToLocBound =
      planeSurfaceObject->localCartesianToBoundLocalDerivative(tgContext,
                                                               globalPosition);
  // For plane surface, this should be identity matrix
  CHECK_CLOSE_ABS(loc3DToLocBound, (ActsMatrix<2, 3>::Identity()), 1e-10);

  // (c) Test the derivative of bound parameters (only test loc0, loc1 here)
  // w.r.t. alignment parameters
  FreeVector derivatives = FreeVector::Zero();
  derivatives.head<3>() = direction;
  const AlignmentToBoundMatrix& alignToBound =
      planeSurfaceObject->alignmentToBoundDerivative(tgContext, globalPosition,
                                                     direction, derivatives);
  const AlignmentToPathMatrix alignToloc0 =
      alignToBound.block<1, 6>(eBoundLoc0, eAlignmentCenter0);
  const AlignmentToPathMatrix alignToloc1 =
      alignToBound.block<1, 6>(eBoundLoc1, eAlignmentCenter0);
  // The expected results
  AlignmentToPathMatrix expAlignToloc0;
  expAlignToloc0 << 0, 0, 1, 0, 0, 2;
  AlignmentToPathMatrix expAlignToloc1;
  expAlignToloc1 << 0, -1, 0, 0, 0, -1;
  // Check if the calculated derivatives are as expected
  CHECK_CLOSE_ABS(alignToloc0, expAlignToloc0, 1e-10);
  CHECK_CLOSE_ABS(alignToloc1, expAlignToloc1, 1e-10);
}

BOOST_AUTO_TEST_SUITE(PlaneSurfaceMerging)

auto logger = Acts::getDefaultLogger("UnitTests", Acts::Logging::VERBOSE);

// Create a test context
GeometryContext gctx = GeometryContext();

auto rBounds = std::make_shared<const RectangleBounds>(1., 2.);

BOOST_AUTO_TEST_CASE(SurfaceOverlap) {
  // Correct orientation, overlapping along merging direction
  Translation3 offsetX{4., 0., 0.};
  Translation3 offsetY{0., 2., 0.};

  Transform3 base(Translation3::Identity());
  Transform3 otherX = base * offsetX;
  Transform3 otherY = base * offsetY;

  auto plane = Surface::makeShared<PlaneSurface>(base, rBounds);
  auto planeX = Surface::makeShared<PlaneSurface>(otherX, rBounds);
  auto planeY = Surface::makeShared<PlaneSurface>(otherY, rBounds);

  BOOST_CHECK_THROW(plane->mergedWith(*planeX, Acts::AxisDirection::AxisX),
                    SurfaceMergingException);
  BOOST_CHECK_THROW(plane->mergedWith(*planeY, Acts::AxisDirection::AxisY),
                    SurfaceMergingException);

  BOOST_CHECK_THROW(planeX->mergedWith(*plane, Acts::AxisDirection::AxisX),
                    SurfaceMergingException);
  BOOST_CHECK_THROW(planeY->mergedWith(*plane, Acts::AxisDirection::AxisY),
                    SurfaceMergingException);
}

BOOST_AUTO_TEST_CASE(SurfaceMisalignmentShift) {
  // Correct orientation, not aligned along orthogonal to merging direction
  Translation3 offsetX{2., 1., 0.};
  Translation3 offsetY{-1., 4., 0.};
  Translation3 offsetZ{0., 4., 1.};

  Transform3 base(Translation3::Identity());
  Transform3 otherX = base * offsetX;
  Transform3 otherY = base * offsetY;
  Transform3 otherZ = base * offsetZ;

  auto plane = Surface::makeShared<PlaneSurface>(base, rBounds);
  auto planeX = Surface::makeShared<PlaneSurface>(otherX, rBounds);
  auto planeY = Surface::makeShared<PlaneSurface>(otherY, rBounds);
  auto planeZ = Surface::makeShared<PlaneSurface>(otherZ, rBounds);

  BOOST_CHECK_THROW(plane->mergedWith(*planeX, Acts::AxisDirection::AxisX),
                    SurfaceMergingException);
  BOOST_CHECK_THROW(plane->mergedWith(*planeY, Acts::AxisDirection::AxisY),
                    SurfaceMergingException);
  BOOST_CHECK_THROW(plane->mergedWith(*planeZ, Acts::AxisDirection::AxisX),
                    SurfaceMergingException);

  BOOST_CHECK_THROW(planeX->mergedWith(*plane, Acts::AxisDirection::AxisX),
                    SurfaceMergingException);
  BOOST_CHECK_THROW(planeY->mergedWith(*plane, Acts::AxisDirection::AxisY),
                    SurfaceMergingException);
  BOOST_CHECK_THROW(planeZ->mergedWith(*plane, Acts::AxisDirection::AxisX),
                    SurfaceMergingException);
}

BOOST_AUTO_TEST_CASE(SurfaceMisalignedAngle) {
  // Correct positioning, rotated in different directions
  Translation3 offsetX{2., 0., 0.};
  Translation3 offsetY{0., 4., 0.};

  double angle = std::numbers::pi / 12;
  Transform3 base(Translation3::Identity());
  Transform3 otherX = base * offsetX * AngleAxis3(angle, Vector3::UnitZ());
  Transform3 otherY = base * offsetY * AngleAxis3(angle, Vector3::UnitY());
  Transform3 otherZ = base * offsetY * AngleAxis3(angle, Vector3::UnitZ());

  auto plane = Surface::makeShared<PlaneSurface>(base, rBounds);
  auto planeX = Surface::makeShared<PlaneSurface>(otherX, rBounds);
  auto planeY = Surface::makeShared<PlaneSurface>(otherY, rBounds);
  auto planeZ = Surface::makeShared<PlaneSurface>(otherZ, rBounds);

  BOOST_CHECK_THROW(plane->mergedWith(*planeX, Acts::AxisDirection::AxisX),
                    SurfaceMergingException);
  BOOST_CHECK_THROW(plane->mergedWith(*planeY, Acts::AxisDirection::AxisY),
                    SurfaceMergingException);
  BOOST_CHECK_THROW(plane->mergedWith(*planeZ, Acts::AxisDirection::AxisY),
                    SurfaceMergingException);

  BOOST_CHECK_THROW(planeX->mergedWith(*plane, Acts::AxisDirection::AxisX),
                    SurfaceMergingException);
  BOOST_CHECK_THROW(planeY->mergedWith(*plane, Acts::AxisDirection::AxisY),
                    SurfaceMergingException);
  BOOST_CHECK_THROW(planeZ->mergedWith(*plane, Acts::AxisDirection::AxisY),
                    SurfaceMergingException);
}

BOOST_AUTO_TEST_CASE(SurfaceDifferentBounds) {
  // Correct orientation and alignment, different bounds lengths along
  // orthogonal to merging direction
  Translation3 offset{2., 0., 0.};

  Transform3 base(Translation3::Identity());
  Transform3 other = base * offset;

  auto plane = Surface::makeShared<PlaneSurface>(base, rBounds);

  auto rBoundsOther = std::make_shared<const RectangleBounds>(2., 4.);
  auto planeOther = Surface::makeShared<PlaneSurface>(other, rBoundsOther);

  BOOST_CHECK_THROW(plane->mergedWith(*planeOther, Acts::AxisDirection::AxisX),
                    SurfaceMergingException);
}

BOOST_AUTO_TEST_CASE(XYDirection) {
  double angle = std::numbers::pi / 12;
  Translation3 offsetX{2., 0., 0.};
  Translation3 offsetY{0., 4., 0.};

  Transform3 base =
      AngleAxis3(angle, Vector3::UnitX()) * Translation3::Identity();
  Transform3 otherX = base * offsetX;
  Transform3 otherY = base * offsetY;

  auto plane = Surface::makeShared<PlaneSurface>(base, rBounds);
  auto planeX = Surface::makeShared<PlaneSurface>(otherX, rBounds);
  auto planeY = Surface::makeShared<PlaneSurface>(otherY, rBounds);

  BOOST_CHECK_THROW(plane->mergedWith(*planeX, Acts::AxisDirection::AxisZ),
                    SurfaceMergingException);

  auto expectedBoundsX = std::make_shared<const RectangleBounds>(2, 2);
  auto [planeXMerged, reversedX] =
      plane->mergedWith(*planeX, Acts::AxisDirection::AxisX, *logger);
  BOOST_REQUIRE_NE(planeXMerged, nullptr);
  BOOST_CHECK(!reversedX);
  BOOST_CHECK_EQUAL(planeXMerged->bounds(), *expectedBoundsX);
  BOOST_CHECK_EQUAL(planeXMerged->center(gctx), base * Vector3::UnitX() * 1);

  auto expectedBoundsY = std::make_shared<const RectangleBounds>(1, 4);
  auto [planeYMerged, reversedY] =
      plane->mergedWith(*planeY, Acts::AxisDirection::AxisY, *logger);
  BOOST_REQUIRE_NE(planeYMerged, nullptr);
  BOOST_CHECK(!reversedY);
  BOOST_CHECK_EQUAL(planeYMerged->bounds(), *expectedBoundsY);
  BOOST_CHECK_EQUAL(planeYMerged->center(gctx), base * Vector3::UnitY() * 2);

  auto [planeXMerged2, reversedX2] =
      planeX->mergedWith(*plane, Acts::AxisDirection::AxisX, *logger);
  BOOST_REQUIRE_NE(planeXMerged2, nullptr);
  BOOST_CHECK(planeXMerged->bounds() == planeXMerged2->bounds());
  BOOST_CHECK(reversedX2);
  BOOST_CHECK_EQUAL(planeXMerged2->bounds(), *expectedBoundsX);
  BOOST_CHECK_EQUAL(planeXMerged2->center(gctx), base * Vector3::UnitX() * 1);

  auto [planeYMerged2, reversedY2] =
      planeY->mergedWith(*plane, Acts::AxisDirection::AxisY, *logger);
  BOOST_REQUIRE_NE(planeYMerged2, nullptr);
  BOOST_CHECK(planeYMerged->bounds() == planeYMerged2->bounds());
  BOOST_CHECK(reversedY2);
  BOOST_CHECK_EQUAL(planeYMerged2->bounds(), *expectedBoundsY);
  BOOST_CHECK_EQUAL(planeYMerged2->center(gctx), base * Vector3::UnitY() * 2);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
