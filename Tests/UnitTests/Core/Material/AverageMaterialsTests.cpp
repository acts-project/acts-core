// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/detail/AverageMaterials.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"

namespace {
using Acts::detail::averageMaterials;

const Acts::MaterialProperties vacuum = Acts::MaterialProperties();
const Acts::MaterialProperties percent = Acts::Test::makePercentSlab();
const Acts::MaterialProperties unit = Acts::Test::makeUnitSlab();

}  // namespace

BOOST_AUTO_TEST_SUITE(AverageMaterials)

// average two identical slabs

BOOST_AUTO_TEST_CASE(VacuumSlabs) {
  auto slab = averageMaterials(vacuum, vacuum);
  BOOST_CHECK(not slab.material());
  BOOST_CHECK_EQUAL(slab.thickness(), 0);
  BOOST_CHECK_EQUAL(slab.thicknessInX0(), 0);
  BOOST_CHECK_EQUAL(slab.thicknessInL0(), 0);
}

BOOST_AUTO_TEST_CASE(PercentSlabs) {
  auto slab = averageMaterials(percent, percent);
  // combining two identical slabs must give the same average material
  BOOST_CHECK(slab.material());
  BOOST_CHECK_EQUAL(slab.material(), percent);
  // thickness-like properties must double
  BOOST_CHECK_EQUAL(slab.thickness(), 2 * percent.thickness());
  BOOST_CHECK_EQUAL(slab.thicknessInX0(), 2 * percent.thicknessInX0());
  BOOST_CHECK_EQUAL(slab.thicknessInL0(), 2 * percent.thicknessInL0());
}

BOOST_AUTO_TEST_CASE(UnitSlabs) {
  auto slab = averageMaterials(unit, unit);
  // combining two identical slabs must give the same average material
  BOOST_CHECK(slab.material());
  BOOST_CHECK_EQUAL(slab.material(), unit);
  // thickness-like properties must double
  BOOST_CHECK_EQUAL(slab.thickness(), 2 * unit.thickness());
  BOOST_CHECK_EQUAL(slab.thicknessInX0(), 2 * unit.thicknessInX0());
  BOOST_CHECK_EQUAL(slab.thicknessInL0(), 2 * unit.thicknessInL0());
}

// average a non-vacuum and an infinitely thin vacuum slab

BOOST_AUTO_TEST_CASE(PercentVacuumSlabs) {
  {
    auto slab = averageMaterials(percent, vacuum);
    BOOST_CHECK(slab.material());
    BOOST_CHECK_EQUAL(slab.material(), percent);
    BOOST_CHECK_EQUAL(slab.thickness(), percent.thickness());
    BOOST_CHECK_EQUAL(slab.thicknessInX0(), percent.thicknessInX0());
    BOOST_CHECK_EQUAL(slab.thicknessInL0(), percent.thicknessInL0());
  }
  // reverse input order
  {
    auto slab = averageMaterials(vacuum, percent);
    BOOST_CHECK(slab.material());
    BOOST_CHECK_EQUAL(slab.material(), percent);
    BOOST_CHECK_EQUAL(slab.thickness(), percent.thickness());
    BOOST_CHECK_EQUAL(slab.thicknessInX0(), percent.thicknessInX0());
    BOOST_CHECK_EQUAL(slab.thicknessInL0(), percent.thicknessInL0());
  }
}

BOOST_AUTO_TEST_CASE(UnitVacuumSlabs) {
  {
    auto slab = averageMaterials(unit, vacuum);
    BOOST_CHECK(slab.material());
    BOOST_CHECK_EQUAL(slab.material(), unit);
    BOOST_CHECK_EQUAL(slab.thickness(), unit.thickness());
    BOOST_CHECK_EQUAL(slab.thicknessInX0(), unit.thicknessInX0());
    BOOST_CHECK_EQUAL(slab.thicknessInL0(), unit.thicknessInL0());
  }
  // reverse input order
  {
    auto slab = averageMaterials(vacuum, unit);
    BOOST_CHECK(slab.material());
    BOOST_CHECK_EQUAL(slab.material(), unit);
    BOOST_CHECK_EQUAL(slab.thickness(), unit.thickness());
    BOOST_CHECK_EQUAL(slab.thicknessInX0(), unit.thicknessInX0());
    BOOST_CHECK_EQUAL(slab.thicknessInL0(), unit.thicknessInL0());
  }
}

BOOST_AUTO_TEST_SUITE_END()
