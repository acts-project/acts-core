// This file is part of the Acts project.
//
// Copyright (C) 2020-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Visualization/ObjVisualization3D.hpp"

#include <iostream>
#include <string>
#include <vector>

#include "TrackingGeometryView3DBase.hpp"
#include "Visualization3DTester.hpp"

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(Visualization)

/// This tests if the corresponding obj output is well formatted
BOOST_AUTO_TEST_CASE(TrackingGeometryView3DObj) {
  ObjVisualization3D obj;
  // Standard test
  bool triangulate = false;
  auto objTest = TrackingGeometryView3DTest::run(obj, triangulate, "");
  auto objErrors = testObjString(objTest, triangulate);
  std::cout << "Surfaces Obj Test    : " << objTest.size()
            << " characters written with " << objErrors.size() << " errors."
            << std::endl;
  BOOST_CHECK(objErrors.empty());
  for (const auto& objerr : objErrors) {
    std::cout << objerr << std::endl;
  }
  // Triangular mesh test
  triangulate = true;
  auto objTest3M = TrackingGeometryView3DTest::run(obj, triangulate, "_3M");
  auto objErrors3M = testObjString(objTest3M, triangulate);
  std::cout << "Surfaces Obj Test 3M : " << objTest3M.size()
            << " characters written with " << objErrors3M.size() << " errors."
            << std::endl;
  BOOST_CHECK(objErrors3M.empty());
  for (const auto& objerr : objErrors3M) {
    std::cout << objerr << std::endl;
  }
}

/*
/// This tests if the corresponding ply output is well formatted
BOOST_AUTO_TEST_CASE(TrackingGeometryView3DPly) {
  PlyVisualization3D ply;
  // Standard test
  bool triangulate = false;
  auto plyTest = TrackingGeometryView3DTest::run(ply, triangulate, "");
  auto plyErrors = testPlyString(plyTest, triangulate);
  std::cout << "Surfaces Ply Test    : " << plyTest.size()
            << " characters written with " << plyErrors.size() << " errors."
            << std::endl;
  BOOST_CHECK_EQUAL(plyErrors.size(), 0);
  for (const auto& plyerr : plyErrors) {
    std::cout << plyerr << std::endl;
  }
  // Triangular mesh test
  triangulate = true;
  auto plyTest3M = TrackingGeometryView3DTest::run(ply, triangulate, "_3M");
  auto plyErrors3M = testPlyString(plyTest3M, triangulate);
  std::cout << "Surfaces Ply Test 3M : " << plyTest3M.size()
            << " characters written with " << plyErrors3M.size() << " errors."
            << std::endl;
  BOOST_CHECK_EQUAL(plyErrors3M.size(), 0);
  for (const auto& plyerr : plyErrors3M) {
    std::cout << plyerr << std::endl;
  }
}
*/

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
