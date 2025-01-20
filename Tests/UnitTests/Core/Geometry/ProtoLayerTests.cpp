// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/RangeXD.hpp"

#include <cmath>
#include <memory>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

namespace Acts::Test::Layers {

GeometryContext tgContext = GeometryContext();

BOOST_AUTO_TEST_SUITE(Geometry)

// Test both const and non-const versions
template <bool IsConst>
void testProtoLayer() {
  using enum AxisDirection;
  using TestProtoLayer = detail::ProtoLayerT<IsConst>;
  using SurfacePtr = typename TestProtoLayer::SurfacePtr;

  // Create a proto layer with 4 surfaces on the x/y grid
  auto recBounds = std::make_shared<RectangleBounds>(3., 6.);

  // Planar definitions to help construct the boundary surfaces
  static const Transform3 planeYZ =
      AngleAxis3(std::numbers::pi / 2., Vector3::UnitY()) *
      AngleAxis3(std::numbers::pi / 2., Vector3::UnitZ()) *
      Transform3::Identity();
  static const Transform3 planeZX =
      AngleAxis3(-std::numbers::pi / 2., Vector3::UnitX()) *
      AngleAxis3(-std::numbers::pi / 2., Vector3::UnitZ()) *
      Transform3::Identity();

  std::vector<
      std::shared_ptr<std::conditional_t<IsConst, const Surface, Surface>>>
      surfaceStore;
  surfaceStore.reserve(100);

  auto createProtoLayer = [&](const Transform3& trf,
                              bool shared = false) -> TestProtoLayer {
    auto atNegX = Surface::makeShared<PlaneSurface>(
        Transform3(trf * Translation3(Vector3(-3., 0., 0.)) * planeYZ),
        recBounds);

    auto atPosX = Surface::makeShared<PlaneSurface>(
        Transform3(trf * Translation3(Vector3(3., 0., 0.)) * planeYZ),
        recBounds);

    auto atNegY = Surface::makeShared<PlaneSurface>(
        Transform3(trf * Translation3(Vector3(0., -3, 0.)) * planeZX),
        recBounds);

    auto atPosY = Surface::makeShared<PlaneSurface>(
        Transform3(trf * Translation3(Vector3(0., 3., 0.)) * planeZX),
        recBounds);

    std::vector<
        std::shared_ptr<std::conditional_t<IsConst, const Surface, Surface>>>
        sharedSurfaces = {atNegX, atNegY, atPosX, atPosY};
    surfaceStore.insert(surfaceStore.begin(), sharedSurfaces.begin(),
                        sharedSurfaces.end());
    if (!shared) {
      std::vector<SurfacePtr> surfaces = {atNegX.get(), atNegY.get(),
                                          atPosX.get(), atPosY.get()};

      return TestProtoLayer(tgContext, surfaces);
    }
    return TestProtoLayer(tgContext, sharedSurfaces);
  };

  // Test 0 - check constructor with surfaces and shared surfaces
  auto pLayerSf = createProtoLayer(Transform3::Identity());
  auto pLayerSfShared = createProtoLayer(Transform3::Identity());

  BOOST_CHECK(pLayerSf.extent.range() == pLayerSfShared.extent.range());
  BOOST_CHECK(pLayerSf.envelope == pLayerSfShared.envelope);

  // CHECK That you have 4 surfaces
  BOOST_CHECK_EQUAL(pLayerSf.surfaces().size(), 4);
  // Add one surface from a detector element (to test thickness)
  auto rB = std::make_shared<RectangleBounds>(30., 60.);

  // Create the detector element
  auto addSurface =
      Surface::makeShared<PlaneSurface>(Transform3::Identity(), rB);

  pLayerSf.add(tgContext, *addSurface.get());
  // CHECK That if you now have 5 surfaces
  BOOST_CHECK_EQUAL(pLayerSf.surfaces().size(), 5);

  // That should invalidate the ranges
  BOOST_CHECK(!(pLayerSf.extent.range() == pLayerSfShared.extent.range()));

  // Test 1 - identity transform
  auto protoLayer = createProtoLayer(Transform3::Identity());

  CHECK_CLOSE_ABS(protoLayer.range(AxisX), 12., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.medium(AxisX), 0., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.min(AxisX), -6., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.max(AxisX), 6., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.range(AxisY), 6., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.medium(AxisY), 0., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.min(AxisY), -3., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.max(AxisY), 3., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.range(AxisZ), 12., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.medium(AxisZ), 0., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.min(AxisZ), -6., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.max(AxisZ), 6., 1e-8);
  CHECK_CLOSE_ABS(protoLayer.max(AxisR), std::hypot(3, 6), 1e-8);
  CHECK_CLOSE_ABS(protoLayer.min(AxisR), 3., 1e-8);

  // Test 2 - rotate around Z-Axis, should leave R, Z untouched,
  // only preserves medium values
  auto protoLayerRot = createProtoLayer(AngleAxis3(-0.345, Vector3::UnitZ()) *
                                        Transform3::Identity());

  BOOST_CHECK_NE(protoLayer.min(AxisX), -6.);
  CHECK_CLOSE_ABS(protoLayerRot.medium(AxisX), 0., 1e-8);
  CHECK_CLOSE_ABS(protoLayerRot.medium(AxisY), 0., 1e-8);
  CHECK_CLOSE_ABS(protoLayerRot.range(AxisZ), 12., 1e-8);
  CHECK_CLOSE_ABS(protoLayerRot.medium(AxisZ), 0., 1e-8);
  CHECK_CLOSE_ABS(protoLayerRot.min(AxisZ), -6., 1e-8);
  CHECK_CLOSE_ABS(protoLayerRot.max(AxisZ), 6., 1e-8);
  CHECK_CLOSE_ABS(protoLayerRot.min(AxisR), 3., 1e-8);
  CHECK_CLOSE_ABS(protoLayerRot.max(AxisR), std::hypot(3, 6), 1e-8);

  std::stringstream sstream;
  protoLayerRot.toStream(sstream);
  std::string oString = R"(ProtoLayer with dimensions (min/max)
Extent in space :
  - value :     AxisX | range = [-6.66104, 6.66104]
  - value :     AxisY | range = [-4.85241, 4.85241]
  - value :     AxisZ | range = [-6, 6]
  - value :     AxisR | range = [3, 6.7082]
  - value :   AxisPhi | range = [-3.02295, 2.33295]
  - value :  AxisRPhi | range = [-20.2785, 15.6499]
  - value : AxisTheta | range = [0.61548, 2.52611]
  - value :   AxisEta | range = [-1.14622, 1.14622]
  - value :   AxisMag | range = [7.34847, 7.34847]
)";
  BOOST_CHECK_EQUAL(sstream.str(), oString);
}

BOOST_AUTO_TEST_CASE(ProtoLayerTests) {
  testProtoLayer<true>();  // Test const version
}

BOOST_AUTO_TEST_CASE(ProtoLayerTestsNonConst) {
  testProtoLayer<false>();  // Test non-const version
}

BOOST_AUTO_TEST_CASE(OrientedLayer) {
  using enum AxisDirection;
  using namespace Acts::UnitLiterals;

  Transform3 base = Transform3::Identity();

  auto recBounds = std::make_shared<RectangleBounds>(3_mm, 6_mm);

  std::vector<std::unique_ptr<DetectorElementBase>> detectorElements;

  auto makeFan = [&](double yrot, double thickness = 0) {
    detectorElements.clear();

    std::size_t nSensors = 8;
    double deltaPhi = 2 * std::numbers::pi / nSensors;
    double r = 20_mm;
    std::vector<std::shared_ptr<const Surface>> surfaces;
    for (std::size_t i = 0; i < nSensors; i++) {
      // Create a fan of sensors

      Transform3 trf = base * AngleAxis3{yrot, Vector3::UnitY()} *
                       AngleAxis3{deltaPhi * i, Vector3::UnitZ()} *
                       Translation3(Vector3::UnitX() * r);

      auto& element = detectorElements.emplace_back(
          std::make_unique<DetectorElementStub>(trf, recBounds, thickness));

      surfaces.push_back(element->surface().getSharedPtr());
    }
    return surfaces;
  };

  std::vector<std::shared_ptr<const Surface>> surfaces = makeFan(0_degree);

  ProtoLayer protoLayer(tgContext, surfaces);

  BOOST_CHECK_EQUAL(protoLayer.surfaces().size(), 8);
  BOOST_CHECK_CLOSE(protoLayer.min(AxisX), -23_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.max(AxisX), 23_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.min(AxisY), -23_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.max(AxisY), 23_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.min(AxisZ), 0_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.max(AxisZ), 0_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.min(AxisR), 17_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.max(AxisR), 23.769728648_mm, 1e-8);

  surfaces = makeFan(45_degree);

  // Do NOT provide rotation matrix: sizing will be affected
  protoLayer = {tgContext, surfaces};

  BOOST_CHECK_EQUAL(protoLayer.surfaces().size(), 8);
  BOOST_CHECK_CLOSE(protoLayer.min(AxisX), -16.26345596_mm, 1e-4);
  BOOST_CHECK_CLOSE(protoLayer.max(AxisX), 16.26345596_mm, 1e-4);
  BOOST_CHECK_CLOSE(protoLayer.min(AxisY), -23_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.max(AxisY), 23_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.min(AxisZ), -16.26345596_mm, 1e-4);
  BOOST_CHECK_CLOSE(protoLayer.max(AxisZ), 16.26345596_mm, 1e-4);

  protoLayer = {tgContext, surfaces,
                Transform3{AngleAxis3{45_degree, Vector3::UnitY()}}.inverse()};

  BOOST_CHECK_EQUAL(protoLayer.surfaces().size(), 8);
  BOOST_CHECK_CLOSE(protoLayer.range(AxisX), 46_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.min(AxisX), -23_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.max(AxisX), 23_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.range(AxisY), 46_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.min(AxisY), -23_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.max(AxisY), 23_mm, 1e-8);
  CHECK_SMALL(protoLayer.range(AxisZ), 1e-14);
  CHECK_SMALL(protoLayer.min(AxisZ), 1e-14);
  CHECK_SMALL(protoLayer.max(AxisZ), 1e-14);

  surfaces = makeFan(0_degree, 10_mm);

  protoLayer = {tgContext, surfaces};

  BOOST_CHECK_EQUAL(protoLayer.surfaces().size(), 8);
  BOOST_CHECK_CLOSE(protoLayer.range(AxisX), 46_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.min(AxisX), -23_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.max(AxisX), 23_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.range(AxisY), 46_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.min(AxisY), -23_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.max(AxisY), 23_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.range(AxisZ), 10_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.min(AxisZ), -5_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.max(AxisZ), 5_mm, 1e-8);

  surfaces = makeFan(45_degree, 10_mm);

  protoLayer = {tgContext, surfaces,
                Transform3{AngleAxis3{45_degree, Vector3::UnitY()}}.inverse()};

  BOOST_CHECK_EQUAL(protoLayer.surfaces().size(), 8);
  BOOST_CHECK_CLOSE(protoLayer.range(AxisX), 46_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.min(AxisX), -23_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.max(AxisX), 23_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.range(AxisY), 46_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.min(AxisY), -23_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.max(AxisY), 23_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.range(AxisZ), 10_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.min(AxisZ), -5_mm, 1e-8);
  BOOST_CHECK_CLOSE(protoLayer.max(AxisZ), 5_mm, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test::Layers
