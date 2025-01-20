// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

#include <cstddef>
#include <vector>

namespace Acts::Test {
BOOST_AUTO_TEST_SUITE(Axes)

BOOST_AUTO_TEST_CASE(equidistant_axis) {
  Axis a(0.0, 10.0, 10u);

  // general binning properties
  BOOST_CHECK_EQUAL(a.getNBins(), 10u);
  BOOST_CHECK_EQUAL(a.getMax(), 10.);
  BOOST_CHECK_EQUAL(a.getMin(), 0.);
  BOOST_CHECK_EQUAL(a.getBinWidth(), 1.);

  // bin index calculation
  BOOST_CHECK_EQUAL(a.getBin(-0.3), 0u);
  BOOST_CHECK_EQUAL(a.getBin(-0.), 1u);
  BOOST_CHECK_EQUAL(a.getBin(0.), 1u);
  BOOST_CHECK_EQUAL(a.getBin(0.7), 1u);
  BOOST_CHECK_EQUAL(a.getBin(1), 2u);
  BOOST_CHECK_EQUAL(a.getBin(1.2), 2u);
  BOOST_CHECK_EQUAL(a.getBin(2.), 3u);
  BOOST_CHECK_EQUAL(a.getBin(2.7), 3u);
  BOOST_CHECK_EQUAL(a.getBin(3.), 4u);
  BOOST_CHECK_EQUAL(a.getBin(3.6), 4u);
  BOOST_CHECK_EQUAL(a.getBin(4.), 5u);
  BOOST_CHECK_EQUAL(a.getBin(4.98), 5u);
  BOOST_CHECK_EQUAL(a.getBin(5.), 6u);
  BOOST_CHECK_EQUAL(a.getBin(5.12), 6u);
  BOOST_CHECK_EQUAL(a.getBin(6.), 7u);
  BOOST_CHECK_EQUAL(a.getBin(6.00001), 7u);
  BOOST_CHECK_EQUAL(a.getBin(7.), 8u);
  BOOST_CHECK_EQUAL(a.getBin(7.5), 8u);
  BOOST_CHECK_EQUAL(a.getBin(8.), 9u);
  BOOST_CHECK_EQUAL(a.getBin(8.1), 9u);
  BOOST_CHECK_EQUAL(a.getBin(9.), 10u);
  BOOST_CHECK_EQUAL(a.getBin(9.999), 10u);
  BOOST_CHECK_EQUAL(a.getBin(10.), 11u);
  BOOST_CHECK_EQUAL(a.getBin(100.3), 11u);

  // lower bin boundaries
  BOOST_CHECK_EQUAL(a.getBinLowerBound(1), 0.);
  BOOST_CHECK_EQUAL(a.getBinLowerBound(2), 1.);
  BOOST_CHECK_EQUAL(a.getBinLowerBound(3), 2.);
  BOOST_CHECK_EQUAL(a.getBinLowerBound(4), 3.);
  BOOST_CHECK_EQUAL(a.getBinLowerBound(5), 4.);
  BOOST_CHECK_EQUAL(a.getBinLowerBound(6), 5.);
  BOOST_CHECK_EQUAL(a.getBinLowerBound(7), 6.);
  BOOST_CHECK_EQUAL(a.getBinLowerBound(8), 7.);
  BOOST_CHECK_EQUAL(a.getBinLowerBound(9), 8.);
  BOOST_CHECK_EQUAL(a.getBinLowerBound(10), 9.);

  // upper bin boundaries
  BOOST_CHECK_EQUAL(a.getBinUpperBound(1), 1.);
  BOOST_CHECK_EQUAL(a.getBinUpperBound(2), 2.);
  BOOST_CHECK_EQUAL(a.getBinUpperBound(3), 3.);
  BOOST_CHECK_EQUAL(a.getBinUpperBound(4), 4.);
  BOOST_CHECK_EQUAL(a.getBinUpperBound(5), 5.);
  BOOST_CHECK_EQUAL(a.getBinUpperBound(6), 6.);
  BOOST_CHECK_EQUAL(a.getBinUpperBound(7), 7.);
  BOOST_CHECK_EQUAL(a.getBinUpperBound(8), 8.);
  BOOST_CHECK_EQUAL(a.getBinUpperBound(9), 9.);
  BOOST_CHECK_EQUAL(a.getBinUpperBound(10), 10.);

  // bin centers
  BOOST_CHECK_EQUAL(a.getBinCenter(1), 0.5);
  BOOST_CHECK_EQUAL(a.getBinCenter(2), 1.5);
  BOOST_CHECK_EQUAL(a.getBinCenter(3), 2.5);
  BOOST_CHECK_EQUAL(a.getBinCenter(4), 3.5);
  BOOST_CHECK_EQUAL(a.getBinCenter(5), 4.5);
  BOOST_CHECK_EQUAL(a.getBinCenter(6), 5.5);
  BOOST_CHECK_EQUAL(a.getBinCenter(7), 6.5);
  BOOST_CHECK_EQUAL(a.getBinCenter(8), 7.5);
  BOOST_CHECK_EQUAL(a.getBinCenter(9), 8.5);
  BOOST_CHECK_EQUAL(a.getBinCenter(10), 9.5);

  // inside check
  BOOST_CHECK(!a.isInside(-0.2));
  BOOST_CHECK(a.isInside(0.));
  BOOST_CHECK(a.isInside(3.));
  BOOST_CHECK(!a.isInside(10.));
  BOOST_CHECK(!a.isInside(12.));
}

BOOST_AUTO_TEST_CASE(variable_axis) {
  Axis a({0, 0.5, 3, 4.5, 6});

  // general binning properties
  BOOST_CHECK_EQUAL(a.getNBins(), 4u);
  BOOST_CHECK_EQUAL(a.getMax(), 6.);
  BOOST_CHECK_EQUAL(a.getMin(), 0.);

  // bin index calculation
  BOOST_CHECK_EQUAL(a.getBin(-0.3), 0u);
  BOOST_CHECK_EQUAL(a.getBin(-0.), 1u);
  BOOST_CHECK_EQUAL(a.getBin(0.), 1u);
  BOOST_CHECK_EQUAL(a.getBin(0.3), 1u);
  BOOST_CHECK_EQUAL(a.getBin(0.5), 2u);
  BOOST_CHECK_EQUAL(a.getBin(1.2), 2u);
  BOOST_CHECK_EQUAL(a.getBin(2.7), 2u);
  BOOST_CHECK_EQUAL(a.getBin(3.), 3u);
  BOOST_CHECK_EQUAL(a.getBin(4.49999), 3u);
  BOOST_CHECK_EQUAL(a.getBin(4.5), 4u);
  BOOST_CHECK_EQUAL(a.getBin(5.12), 4u);
  BOOST_CHECK_EQUAL(a.getBin(6.), 5u);
  BOOST_CHECK_EQUAL(a.getBin(6.00001), 5u);
  BOOST_CHECK_EQUAL(a.getBin(7.5), 5u);

  // lower bin boundaries
  BOOST_CHECK_EQUAL(a.getBinLowerBound(1), 0.);
  BOOST_CHECK_EQUAL(a.getBinLowerBound(2), 0.5);
  BOOST_CHECK_EQUAL(a.getBinLowerBound(3), 3.);
  BOOST_CHECK_EQUAL(a.getBinLowerBound(4), 4.5);

  // upper bin boundaries
  BOOST_CHECK_EQUAL(a.getBinUpperBound(1), 0.5);
  BOOST_CHECK_EQUAL(a.getBinUpperBound(2), 3.);
  BOOST_CHECK_EQUAL(a.getBinUpperBound(3), 4.5);
  BOOST_CHECK_EQUAL(a.getBinUpperBound(4), 6.);

  // bin centers
  BOOST_CHECK_EQUAL(a.getBinCenter(1), 0.25);
  BOOST_CHECK_EQUAL(a.getBinCenter(2), 1.75);
  BOOST_CHECK_EQUAL(a.getBinCenter(3), 3.75);
  BOOST_CHECK_EQUAL(a.getBinCenter(4), 5.25);

  // inside check
  BOOST_CHECK(!a.isInside(-0.2));
  BOOST_CHECK(a.isInside(0.));
  BOOST_CHECK(a.isInside(3.));
  BOOST_CHECK(!a.isInside(6.));
  BOOST_CHECK(!a.isInside(12.));
}

BOOST_AUTO_TEST_CASE(open_axis) {
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> a(0, 10, 10);

  // normal inside
  BOOST_CHECK_EQUAL(a.getBin(0.5), 1u);
  BOOST_CHECK_EQUAL(a.getBin(9.5), 10u);

  // out of bounds, but is open
  // -> should clamp
  BOOST_CHECK_EQUAL(a.getBin(-0.5), 1u);
  BOOST_CHECK_EQUAL(a.getBin(10.5), 10u);

  Axis<AxisType::Variable, AxisBoundaryType::Bound> b(
      {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10});

  // normal inside
  BOOST_CHECK_EQUAL(b.getBin(0.5), 1u);
  BOOST_CHECK_EQUAL(b.getBin(9.5), 10u);

  // out of bounds, but is open
  // -> should clamp
  BOOST_CHECK_EQUAL(b.getBin(-0.5), 1u);
  BOOST_CHECK_EQUAL(b.getBin(10.5), 10u);
}

BOOST_AUTO_TEST_CASE(closed_axis) {
  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> a(0, 10, 10);

  // normal inside
  BOOST_CHECK_EQUAL(a.getBin(0.5), 1u);
  BOOST_CHECK_EQUAL(a.getBin(9.5), 10u);

  // out of bounds, but is closed
  // -> should wrap to opposite side bin
  BOOST_CHECK_EQUAL(a.getBin(-0.5), 10u);
  BOOST_CHECK_EQUAL(a.getBin(10.5), 1u);

  Axis<AxisType::Variable, AxisBoundaryType::Closed> b(
      {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10});

  // normal inside
  BOOST_CHECK_EQUAL(b.getBin(0.5), 1u);
  BOOST_CHECK_EQUAL(b.getBin(9.5), 10u);

  // out of bounds, but is closed
  // -> should wrap to opposite side bin
  BOOST_CHECK_EQUAL(b.getBin(-0.5), 10u);
  BOOST_CHECK_EQUAL(b.getBin(10.5), 1u);
}

BOOST_AUTO_TEST_CASE(neighborhood) {
  using bins_t = std::vector<std::size_t>;
  Axis<AxisType::Equidistant, AxisBoundaryType::Open> a1(0.0, 1.0, 10u);

  BOOST_CHECK(a1.neighborHoodIndices(0, 1).collect() == bins_t({0, 1}));
  BOOST_CHECK(a1.neighborHoodIndices(1, 1).collect() == bins_t({0, 1, 2}));
  BOOST_CHECK(a1.neighborHoodIndices(11, 1).collect() == bins_t({10, 11}));
  BOOST_CHECK(a1.neighborHoodIndices(10, 1).collect() == bins_t({9, 10, 11}));
  BOOST_CHECK(a1.neighborHoodIndices(5, 1).collect() == bins_t({4, 5, 6}));
  BOOST_CHECK(a1.neighborHoodIndices(5, {-1, 0}).collect() == bins_t({4, 5}));
  BOOST_CHECK(a1.neighborHoodIndices(5, {0, 1}).collect() == bins_t({5, 6}));

  BOOST_CHECK(a1.neighborHoodIndices(0, 2).collect() == bins_t({0, 1, 2}));
  BOOST_CHECK(a1.neighborHoodIndices(1, 2).collect() == bins_t({0, 1, 2, 3}));
  BOOST_CHECK(a1.neighborHoodIndices(11, 2).collect() == bins_t({9, 10, 11}));
  BOOST_CHECK(a1.neighborHoodIndices(10, 2).collect() ==
              bins_t({8, 9, 10, 11}));
  BOOST_CHECK(a1.neighborHoodIndices(5, 2).collect() ==
              bins_t({3, 4, 5, 6, 7}));

  Axis<AxisType::Variable, AxisBoundaryType::Open> a2(
      {0.0, 2.0, 4.0, 9.0, 10.0});
  BOOST_CHECK(a2.neighborHoodIndices(0, 1).collect() == bins_t({0, 1}));
  BOOST_CHECK(a2.neighborHoodIndices(1, 1).collect() == bins_t({0, 1, 2}));
  BOOST_CHECK(a2.neighborHoodIndices(5, 1).collect() == bins_t({4, 5}));
  BOOST_CHECK(a2.neighborHoodIndices(4, 1).collect() == bins_t({3, 4, 5}));
  BOOST_CHECK(a2.neighborHoodIndices(4, {-1, 0}).collect() == bins_t({3, 4}));
  BOOST_CHECK(a2.neighborHoodIndices(2, 1).collect() == bins_t({1, 2, 3}));
  BOOST_CHECK(a2.neighborHoodIndices(2, {0, 1}).collect() == bins_t({2, 3}));

  BOOST_CHECK(a2.neighborHoodIndices(0, 2).collect() == bins_t({0, 1, 2}));
  BOOST_CHECK(a2.neighborHoodIndices(1, 2).collect() == bins_t({0, 1, 2, 3}));
  BOOST_CHECK(a2.neighborHoodIndices(5, 2).collect() == bins_t({3, 4, 5}));
  BOOST_CHECK(a2.neighborHoodIndices(4, 2).collect() == bins_t({2, 3, 4, 5}));
  BOOST_CHECK(a2.neighborHoodIndices(3, 2).collect() ==
              bins_t({1, 2, 3, 4, 5}));

  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> a3(0.0, 1.0, 10u);

  BOOST_CHECK(a3.neighborHoodIndices(0, 1).collect() == bins_t({}));
  BOOST_CHECK(a3.neighborHoodIndices(1, 1).collect() == bins_t({1, 2}));
  BOOST_CHECK(a3.neighborHoodIndices(11, 1).collect() == bins_t({}));
  BOOST_CHECK(a3.neighborHoodIndices(10, 1).collect() == bins_t({9, 10}));
  BOOST_CHECK(a3.neighborHoodIndices(10, {0, 1}).collect() == bins_t({10}));
  BOOST_CHECK(a3.neighborHoodIndices(5, 1).collect() == bins_t({4, 5, 6}));
  BOOST_CHECK(a3.neighborHoodIndices(5, {-1, 0}).collect() == bins_t({4, 5}));
  BOOST_CHECK(a3.neighborHoodIndices(5, {0, 1}).collect() == bins_t({5, 6}));

  BOOST_CHECK(a3.neighborHoodIndices(0, 2).collect() == bins_t({}));
  BOOST_CHECK(a3.neighborHoodIndices(1, 2).collect() == bins_t({1, 2, 3}));
  BOOST_CHECK(a3.neighborHoodIndices(11, 2).collect() == bins_t({}));
  BOOST_CHECK(a3.neighborHoodIndices(10, 2).collect() == bins_t({8, 9, 10}));
  BOOST_CHECK(a3.neighborHoodIndices(5, 2).collect() ==
              bins_t({3, 4, 5, 6, 7}));

  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> a4(0.0, 1.0, 10u);

  BOOST_CHECK(a4.neighborHoodIndices(0, 1).collect() == bins_t({}));
  BOOST_CHECK(a4.neighborHoodIndices(1, 1).collect() == bins_t({10, 1, 2}));
  BOOST_CHECK(a4.neighborHoodIndices(11, 1).collect() == bins_t({}));
  BOOST_CHECK(a4.neighborHoodIndices(10, 1).collect() == bins_t({9, 10, 1}));
  BOOST_CHECK(a4.neighborHoodIndices(10, {0, 1}).collect() == bins_t({10, 1}));
  BOOST_CHECK(a4.neighborHoodIndices(5, 1).collect() == bins_t({4, 5, 6}));
  BOOST_CHECK(a4.neighborHoodIndices(5, {-1, 0}).collect() == bins_t({4, 5}));
  BOOST_CHECK(a4.neighborHoodIndices(5, {0, 1}).collect() == bins_t({5, 6}));

  BOOST_CHECK(a4.neighborHoodIndices(0, 2).collect() == bins_t({}));
  BOOST_CHECK(a4.neighborHoodIndices(1, 2).collect() ==
              bins_t({9, 10, 1, 2, 3}));
  BOOST_CHECK(a4.neighborHoodIndices(11, 2).collect() == bins_t({}));
  BOOST_CHECK(a4.neighborHoodIndices(10, 2).collect() ==
              bins_t({8, 9, 10, 1, 2}));
  BOOST_CHECK(a4.neighborHoodIndices(5, 2).collect() ==
              bins_t({3, 4, 5, 6, 7}));
  BOOST_CHECK(a4.neighborHoodIndices(3, 2).collect() ==
              bins_t({1, 2, 3, 4, 5}));

  Axis<AxisType::Variable, AxisBoundaryType::Bound> a5(
      {0.0, 2.0, 4.0, 9.0, 9.5, 10.0});
  BOOST_CHECK(a5.neighborHoodIndices(0, 1).collect() == bins_t({}));
  BOOST_CHECK(a5.neighborHoodIndices(1, 1).collect() == bins_t({1, 2}));
  BOOST_CHECK(a5.neighborHoodIndices(6, 1).collect() == bins_t({}));
  BOOST_CHECK(a5.neighborHoodIndices(5, 1).collect() == bins_t({4, 5}));
  BOOST_CHECK(a5.neighborHoodIndices(5, {0, 1}).collect() == bins_t({5}));
  BOOST_CHECK(a5.neighborHoodIndices(2, 1).collect() == bins_t({1, 2, 3}));
  BOOST_CHECK(a5.neighborHoodIndices(2, {-1, 0}).collect() == bins_t({1, 2}));
  BOOST_CHECK(a5.neighborHoodIndices(2, {0, 1}).collect() == bins_t({2, 3}));

  BOOST_CHECK(a5.neighborHoodIndices(0, 2).collect() == bins_t({}));
  BOOST_CHECK(a5.neighborHoodIndices(1, 2).collect() == bins_t({1, 2, 3}));
  BOOST_CHECK(a5.neighborHoodIndices(6, 2).collect() == bins_t({}));
  BOOST_CHECK(a5.neighborHoodIndices(5, 2).collect() == bins_t({3, 4, 5}));
  BOOST_CHECK(a5.neighborHoodIndices(3, 2).collect() ==
              bins_t({1, 2, 3, 4, 5}));

  Axis<AxisType::Variable, AxisBoundaryType::Closed> a6(
      {0.0, 2.0, 4.0, 9.0, 9.5, 10.0});
  BOOST_CHECK(a6.neighborHoodIndices(0, 1).collect() == bins_t({}));
  BOOST_CHECK(a6.neighborHoodIndices(1, 1).collect() == bins_t({5, 1, 2}));
  BOOST_CHECK(a6.neighborHoodIndices(6, 1).collect() == bins_t({}));
  BOOST_CHECK(a6.neighborHoodIndices(5, 1).collect() == bins_t({4, 5, 1}));
  BOOST_CHECK(a6.neighborHoodIndices(5, {0, 1}).collect() == bins_t({5, 1}));
  BOOST_CHECK(a6.neighborHoodIndices(2, 1).collect() == bins_t({1, 2, 3}));
  BOOST_CHECK(a6.neighborHoodIndices(2, {-1, 0}).collect() == bins_t({1, 2}));
  BOOST_CHECK(a6.neighborHoodIndices(2, {0, 1}).collect() == bins_t({2, 3}));

  BOOST_CHECK(a6.neighborHoodIndices(0, 2).collect() == bins_t({}));
  BOOST_CHECK(a6.neighborHoodIndices(1, 2).collect() ==
              bins_t({4, 5, 1, 2, 3}));
  BOOST_CHECK(a6.neighborHoodIndices(6, 2).collect() == bins_t({}));
  BOOST_CHECK(a6.neighborHoodIndices(5, 2).collect() ==
              bins_t({3, 4, 5, 1, 2}));
  BOOST_CHECK(a6.neighborHoodIndices(3, 2).collect() ==
              bins_t({1, 2, 3, 4, 5}));
  BOOST_CHECK(a6.neighborHoodIndices(3, {0, 2}).collect() == bins_t({3, 4, 5}));

  BOOST_CHECK(a6.neighborHoodIndices(1, 3).collect() ==
              bins_t({1, 2, 3, 4, 5}));
  BOOST_CHECK(a6.neighborHoodIndices(5, 3).collect() ==
              bins_t({1, 2, 3, 4, 5}));
}

BOOST_AUTO_TEST_CASE(wrapBin) {
  Axis<AxisType::Equidistant, AxisBoundaryType::Open> a1(0.0, 1.0, 10u);
  BOOST_CHECK_EQUAL(a1.wrapBin(0), 0u);
  BOOST_CHECK_EQUAL(a1.wrapBin(1), 1u);
  BOOST_CHECK_EQUAL(a1.wrapBin(-1), 0u);
  BOOST_CHECK_EQUAL(a1.wrapBin(10), 10u);
  BOOST_CHECK_EQUAL(a1.wrapBin(11), 11u);
  BOOST_CHECK_EQUAL(a1.wrapBin(12), 11u);

  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> a2(0.0, 1.0, 10u);
  BOOST_CHECK_EQUAL(a2.wrapBin(0), 1u);
  BOOST_CHECK_EQUAL(a2.wrapBin(1), 1u);
  BOOST_CHECK_EQUAL(a2.wrapBin(-1), 1u);
  BOOST_CHECK_EQUAL(a2.wrapBin(10), 10u);
  BOOST_CHECK_EQUAL(a2.wrapBin(11), 10u);
  BOOST_CHECK_EQUAL(a2.wrapBin(12), 10u);

  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> a3(0.0, 1.0, 10u);
  BOOST_CHECK_EQUAL(a3.wrapBin(0), 10u);
  BOOST_CHECK_EQUAL(a3.wrapBin(1), 1u);
  BOOST_CHECK_EQUAL(a3.wrapBin(-1), 9u);
  BOOST_CHECK_EQUAL(a3.wrapBin(10), 10u);
  BOOST_CHECK_EQUAL(a3.wrapBin(11), 1u);
  BOOST_CHECK_EQUAL(a3.wrapBin(12), 2u);

  Axis<AxisType::Variable, AxisBoundaryType::Open> a4(
      {0.0, 2.0, 4.0, 9.0, 10.0});
  BOOST_CHECK_EQUAL(a4.wrapBin(0), 0u);
  BOOST_CHECK_EQUAL(a4.wrapBin(1), 1u);
  BOOST_CHECK_EQUAL(a4.wrapBin(-1), 0u);
  BOOST_CHECK_EQUAL(a4.wrapBin(4), 4u);
  BOOST_CHECK_EQUAL(a4.wrapBin(5), 5u);
  BOOST_CHECK_EQUAL(a4.wrapBin(6), 5u);

  Axis<AxisType::Variable, AxisBoundaryType::Bound> a5(
      {0.0, 2.0, 4.0, 9.0, 9.5, 10.0});
  BOOST_CHECK_EQUAL(a5.wrapBin(0), 1u);
  BOOST_CHECK_EQUAL(a5.wrapBin(1), 1u);
  BOOST_CHECK_EQUAL(a5.wrapBin(-1), 1u);
  BOOST_CHECK_EQUAL(a5.wrapBin(4), 4u);
  BOOST_CHECK_EQUAL(a5.wrapBin(5), 5u);
  BOOST_CHECK_EQUAL(a5.wrapBin(6), 5u);

  Axis<AxisType::Variable, AxisBoundaryType::Closed> a6(
      {0.0, 2.0, 4.0, 9.0, 9.5, 10.0});
  BOOST_CHECK_EQUAL(a6.wrapBin(0), 5u);
  BOOST_CHECK_EQUAL(a6.wrapBin(1), 1u);
  BOOST_CHECK_EQUAL(a6.wrapBin(-1), 4u);
  BOOST_CHECK_EQUAL(a6.wrapBin(4), 4u);
  BOOST_CHECK_EQUAL(a6.wrapBin(5), 5u);
  BOOST_CHECK_EQUAL(a6.wrapBin(6), 1u);
  BOOST_CHECK_EQUAL(a6.wrapBin(7), 2u);
}

BOOST_AUTO_TEST_CASE(AxisTypeDeduction) {
  auto eqOpen = Axis{0.0, 10., 10};
  static_assert(
      std::is_same_v<decltype(eqOpen),
                     Axis<AxisType::Equidistant, AxisBoundaryType::Open>>);
  auto eqBound = Axis{AxisBound, 0.0, 10., 10};
  static_assert(
      std::is_same_v<decltype(eqBound),
                     Axis<AxisType::Equidistant, AxisBoundaryType::Bound>>);
  auto eqClosed = Axis{AxisClosed, 0.0, 10., 10};
  static_assert(
      std::is_same_v<decltype(eqClosed),
                     Axis<AxisType::Equidistant, AxisBoundaryType::Closed>>);

  auto varOpen = Axis{{0, 1, 2., 3, 4}};
  static_assert(
      std::is_same_v<decltype(varOpen),
                     Axis<AxisType::Variable, AxisBoundaryType::Open>>);
  auto varBound = Axis{AxisBound, {0, 1, 2., 3, 4}};
  static_assert(
      std::is_same_v<decltype(varBound),
                     Axis<AxisType::Variable, AxisBoundaryType::Bound>>);
  auto varClosed = Axis{AxisClosed, {0, 1, 2., 3, 4}};
  static_assert(
      std::is_same_v<decltype(varClosed),
                     Axis<AxisType::Variable, AxisBoundaryType::Closed>>);
}

BOOST_AUTO_TEST_CASE(AxisVisit) {
  using enum AxisBoundaryType;
  using enum AxisType;

  auto eqOpen = Axis{0.0, 10., 10};
  eqOpen.visit([](const auto& axis) {
    BOOST_CHECK((
        std::is_same_v<std::decay_t<decltype(axis)>, Axis<Equidistant, Open>>));
  });

  auto eqBound = Axis{AxisBound, 0.0, 10., 10};
  eqBound.visit([](const auto& axis) {
    BOOST_CHECK((std::is_same_v<std::decay_t<decltype(axis)>,
                                Axis<Equidistant, Bound>>));
  });

  auto eqClosed = Axis{AxisClosed, 0.0, 10., 10};
  eqClosed.visit([](const auto& axis) {
    BOOST_CHECK((std::is_same_v<std::decay_t<decltype(axis)>,
                                Axis<Equidistant, Closed>>));
  });

  auto varOpen = Axis{{0, 1, 2., 3, 4}};
  varOpen.visit([](const auto& axis) {
    BOOST_CHECK(
        (std::is_same_v<std::decay_t<decltype(axis)>, Axis<Variable, Open>>));
  });

  auto varBound = Axis{AxisBound, {0, 1, 2., 3, 4}};
  varBound.visit([](const auto& axis) {
    BOOST_CHECK(
        (std::is_same_v<std::decay_t<decltype(axis)>, Axis<Variable, Bound>>));
  });

  auto varClosed = Axis{AxisClosed, {0, 1, 2., 3, 4}};
  varClosed.visit([](const auto& axis) {
    BOOST_CHECK(
        (std::is_same_v<std::decay_t<decltype(axis)>, Axis<Variable, Closed>>));
  });

  std::vector<double> edges =
      varClosed.visit([](const auto& axis) { return axis.getBinEdges(); });
  BOOST_CHECK_EQUAL(edges.size(), varClosed.getBinEdges().size());

  // Test return values from visit method with type-dependent values
  int typeValue = eqOpen.visit([](const auto& axis) {
    if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                 Axis<Equidistant, Open>>) {
      return 1;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Equidistant, Bound>>) {
      return 2;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Equidistant, Closed>>) {
      return 3;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Variable, Open>>) {
      return 4;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Variable, Bound>>) {
      return 5;
    } else {
      return 6;  // Variable, Closed
    }
  });
  BOOST_CHECK_EQUAL(typeValue, 1);  // Should be Equidistant, Open

  typeValue = eqBound.visit([](const auto& axis) {
    if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                 Axis<Equidistant, Open>>) {
      return 1;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Equidistant, Bound>>) {
      return 2;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Equidistant, Closed>>) {
      return 3;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Variable, Open>>) {
      return 4;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Variable, Bound>>) {
      return 5;
    } else {
      return 6;  // Variable, Closed
    }
  });
  BOOST_CHECK_EQUAL(typeValue, 2);  // Should be Equidistant, Bound

  typeValue = eqClosed.visit([](const auto& axis) {
    if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                 Axis<Equidistant, Open>>) {
      return 1;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Equidistant, Bound>>) {
      return 2;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Equidistant, Closed>>) {
      return 3;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Variable, Open>>) {
      return 4;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Variable, Bound>>) {
      return 5;
    } else {
      return 6;  // Variable, Closed
    }
  });
  BOOST_CHECK_EQUAL(typeValue, 3);  // Should be Equidistant, Closed

  typeValue = varOpen.visit([](const auto& axis) {
    if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                 Axis<Equidistant, Open>>) {
      return 1;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Equidistant, Bound>>) {
      return 2;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Equidistant, Closed>>) {
      return 3;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Variable, Open>>) {
      return 4;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Variable, Bound>>) {
      return 5;
    } else {
      return 6;  // Variable, Closed
    }
  });
  BOOST_CHECK_EQUAL(typeValue, 4);  // Should be Variable, Open

  typeValue = varBound.visit([](const auto& axis) {
    if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                 Axis<Equidistant, Open>>) {
      return 1;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Equidistant, Bound>>) {
      return 2;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Equidistant, Closed>>) {
      return 3;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Variable, Open>>) {
      return 4;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Variable, Bound>>) {
      return 5;
    } else {
      return 6;  // Variable, Closed
    }
  });
  BOOST_CHECK_EQUAL(typeValue, 5);  // Should be Variable, Bound

  typeValue = varClosed.visit([](const auto& axis) {
    if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                 Axis<Equidistant, Open>>) {
      return 1;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Equidistant, Bound>>) {
      return 2;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Equidistant, Closed>>) {
      return 3;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Variable, Open>>) {
      return 4;
    } else if constexpr (std::is_same_v<std::decay_t<decltype(axis)>,
                                        Axis<Variable, Bound>>) {
      return 5;
    } else {
      return 6;  // Variable, Closed
    }
  });
  BOOST_CHECK_EQUAL(typeValue, 6);  // Should be Variable, Closed

  // Test return value using axis properties
  double minValue =
      eqOpen.visit([](const auto& axis) { return axis.getMin(); });
  BOOST_CHECK_EQUAL(minValue, 0.0);

  double maxValue =
      eqBound.visit([](const auto& axis) { return axis.getMax(); });
  BOOST_CHECK_EQUAL(maxValue, 10.0);

  std::size_t nBins =
      varClosed.visit([](const auto& axis) { return axis.getNBins(); });
  BOOST_CHECK_EQUAL(nBins, 4u);
}

BOOST_AUTO_TEST_CASE(IAxis_Factories) {
  using enum AxisType;
  using enum AxisBoundaryType;

  // Equidistan: Bound, Open, Closed
  auto eb = IAxis::create(Bound, 0.0, 10., 10);
  BOOST_CHECK_EQUAL(eb->getType(), Equidistant);
  BOOST_CHECK_EQUAL(eb->getBoundaryType(), Bound);

  auto eo = IAxis::create(Open, 0.0, 10., 10);
  BOOST_CHECK_EQUAL(eo->getType(), Equidistant);
  BOOST_CHECK_EQUAL(eo->getBoundaryType(), Open);

  auto ec = IAxis::create(Closed, 0.0, 10., 10);
  BOOST_CHECK_EQUAL(ec->getType(), Equidistant);
  BOOST_CHECK_EQUAL(ec->getBoundaryType(), Closed);

  // Variable: Bound, Open, Closed
  auto vb = IAxis::create(Bound, {0, 1, 2., 3, 4});
  BOOST_CHECK_EQUAL(vb->getType(), Variable);
  BOOST_CHECK_EQUAL(vb->getBoundaryType(), Bound);

  auto vo = IAxis::create(Open, {0, 1, 2., 3, 4});
  BOOST_CHECK_EQUAL(vo->getType(), Variable);
  BOOST_CHECK_EQUAL(vo->getBoundaryType(), Open);

  auto vc = IAxis::create(Closed, {0, 1, 2., 3, 4});
  BOOST_CHECK_EQUAL(vc->getType(), Variable);
  BOOST_CHECK_EQUAL(vc->getBoundaryType(), Closed);

  // Invalid constructors
  // min > max
  BOOST_CHECK_THROW(IAxis::create(Bound, 10., 0., 3.), std::invalid_argument);
  // nBins = 0
  BOOST_CHECK_THROW(IAxis::create(Bound, 0., 10., 0.), std::invalid_argument);
  // #edges < 2
  BOOST_CHECK_THROW(IAxis::create(Bound, std::vector<double>{2.}),
                    std::invalid_argument);
  // edges not ordered
  BOOST_CHECK_THROW(IAxis::create(Bound, std::vector<double>{2., 1.5, 1.}),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(Output) {
  std::stringstream ss;

  Axis a{AxisBound, 0.0, 10., 10};
  Axis b{AxisBound, {0.0, 10., 11}};

  ss << a;

  BOOST_CHECK_EQUAL(ss.str(), "Axis<Equidistant, Bound>(0, 10, 10)");

  ss.str("");

  const IAxis& ia = a;

  ss << ia;

  BOOST_CHECK_EQUAL(ss.str(), "Axis<Equidistant, Bound>(0, 10, 10)");

  ss.str("");

  ss << b;

  BOOST_CHECK_EQUAL(ss.str(), "Axis<Variable, Bound>(0, 10, 11)");
}

BOOST_AUTO_TEST_CASE(Equality) {
  Axis a{AxisBound, 0.0, 10., 10};
  Axis b{AxisClosed, 0.0, 10., 10};

  BOOST_CHECK_EQUAL(a, a);
  BOOST_CHECK_NE(a, b);

  const IAxis& ia = a;
  const IAxis& ib = b;

  BOOST_CHECK_EQUAL(ia, ia);
  BOOST_CHECK_NE(ia, ib);
  BOOST_CHECK_NE(ia, b);
  BOOST_CHECK_NE(a, ib);
  BOOST_CHECK_EQUAL(a, ia);
  BOOST_CHECK_EQUAL(b, ib);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
