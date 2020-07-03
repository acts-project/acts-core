// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/Volume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

namespace tt = boost::test_tools;

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_CASE(VolumeTest) {
	using namespace Acts::UnitLiterals;
	double eps = std::numeric_limits<double>::epsilon();
	
	// Build a translation
	Vector3D translation{1_mm, 2_mm, 3_mm};
	
	// Build a translation
	ActsMatrixD<3, 3> rotation = RotationMatrix3D::Identity();
	double rotationAngle = 60_degree;
	Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
	Vector3D yPos(0., 1., 0.);
	Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
	rotation.col(0) = xPos;
	rotation.col(1) = yPos;
	rotation.col(2) = zPos;
	
	// Build a transform
	Transform3D transform(Transform3D::Identity() * rotation);
	transform.translation() = translation;
	// Build the bounds
	CuboidVolumeBounds bounds(4_mm, 5_mm, 6_mm);
	
	// Build and test the volume
	Volume volume(std::make_shared<const Transform3D>(transform), std::make_shared<const CuboidVolumeBounds>(bounds));
	BOOST_TEST(volume.transform().matrix() == transform.matrix());
	CHECK_CLOSE_ABS(volume.itransform().matrix(), transform.inverse().matrix(), eps);
	BOOST_TEST(volume.center() == translation);
	auto vBounds = static_cast<const decltype(bounds)*>(&volume.volumeBounds());
	BOOST_TEST(*vBounds == bounds);
	
	// Build and test a shifted volume
	Transform3D shift(Transform3D::Identity());
	Vector3D shiftTranslation{-4_mm, -5_mm, -6_mm};
	shift.translation() = shiftTranslation;
	Volume volumeShift(volume, &shift);
	BOOST_TEST(volumeShift.center() == volume.center() + shiftTranslation);
	BOOST_TEST(volumeShift.transform().rotation() == volume.transform().rotation());
	
	// Test a bounding box from the volume
	// Vector3D envelope{7_mm, 8_mm, 9_mm};
	// BOOST_TEST(volume.boundingBox(envelope) == bounds.boundingBox(&transform, envelope, &volume));
	
	// Inside/Outside check
	BOOST_TEST(volume.inside(translation));
	BOOST_TEST(!volume.inside({10_mm, 2_mm, 3_mm})); 
	BOOST_TEST(volume.inside({10_mm, 2_mm, 3_mm}, 2_mm));
	
	// Binning test
	GeometryContext gctx;
	BOOST_TEST(volume.binningPosition(gctx, binX) == volume.center());
}
}  // namespace Test
}  // namespace Acts
