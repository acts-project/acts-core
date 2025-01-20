// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

BOOST_AUTO_TEST_SUITE(ProtoAxis)

BOOST_AUTO_TEST_CASE(EquidistantProtoAxis) {
  using enum Acts::AxisBoundaryType;
  using enum Acts::AxisDirection;
  using enum Acts::AxisType;

  // Bound, equidistant axis
  Acts::ProtoAxis epab(AxisX, Bound, 0.0, 1.0, 10);

  // Direct access
  BOOST_CHECK_EQUAL(epab.getAxisDirection(), AxisX);
  BOOST_CHECK(!epab.isAutorange());

  // Access via IAxis
  BOOST_CHECK(epab.getAxis().isEquidistant());

  BOOST_CHECK(!epab.getAxis().isVariable());

  auto edges = epab.getAxis().getBinEdges();
  BOOST_CHECK_EQUAL(edges.size(), 11);

  BOOST_CHECK_EQUAL(epab.getAxis().getType(), Equidistant);

  BOOST_CHECK_EQUAL(epab.getAxis().getBoundaryType(), Bound);

  BOOST_CHECK_EQUAL(epab.getAxis().getNBins(), 10);

  CHECK_CLOSE_ABS(epab.getAxis().getMin(), 0.0, 1e-15);

  CHECK_CLOSE_ABS(epab.getAxis().getMax(), 1.0, 1e-15);

  std::string rString =
      "ProtoAxis: 10 bins in AxisX, equidistant within [0, 1]";
  std::string oString = epab.toString();
  BOOST_CHECK_EQUAL(rString, oString);

  // Create a grid from a single proto axis
  auto grid1D = Acts::makeGrid<double>(epab);
  BOOST_CHECK(grid1D != nullptr);
  BOOST_CHECK_EQUAL(grid1D->axes().size(), 1);
  auto axis1D =
      dynamic_cast<const Acts::Axis<Acts::AxisType::Equidistant, Bound>*>(
          grid1D->axes().front());
  BOOST_CHECK(axis1D != nullptr);

  // Open, equidistant axis
  Acts::ProtoAxis epao(AxisY, Open, 0., 2.0, 10.);

  // Create a 2D grid from a two proto axes
  auto grid2D = Acts::makeGrid<double>(epab, epao);
  BOOST_CHECK(grid2D != nullptr);
  auto grid2Daxes = grid2D->axes();
  BOOST_CHECK_EQUAL(grid2Daxes.size(), 2);
  auto axis2D1 =
      dynamic_cast<const Acts::Axis<Acts::AxisType::Equidistant, Bound>*>(
          grid2Daxes[0]);
  BOOST_CHECK(axis2D1 != nullptr);
  auto axis2D2 =
      dynamic_cast<const Acts::Axis<Acts::AxisType::Equidistant, Open>*>(
          grid2Daxes[1]);
  BOOST_CHECK(axis2D2 != nullptr);

  BOOST_CHECK_EQUAL(epao.getAxis().getBoundaryType(), Open);

  // Invalid constructor, min > max
  BOOST_CHECK_THROW(Acts::ProtoAxis(AxisZ, Bound, 1.0, 0.0, 10),
                    std::invalid_argument);

  // Invalid constructor, nbins < 1
  BOOST_CHECK_THROW(Acts::ProtoAxis(AxisZ, Bound, 0.0, 1.0, 0),
                    std::invalid_argument);

  // Invalid constructor, closed with something else than phi or rphi
  BOOST_CHECK_THROW(Acts::ProtoAxis(AxisZ, Closed, 0.0, 1.0, 10),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(AutorangeProtoAxis) {
  using enum Acts::AxisBoundaryType;
  using enum Acts::AxisDirection;
  using enum Acts::AxisType;

  // Bound, equidistant axis, autorange
  Acts::ProtoAxis epa(AxisX, Bound, 10);

  // Direct access
  BOOST_CHECK_EQUAL(epa.getAxisDirection(), AxisX);

  BOOST_CHECK(epa.isAutorange());

  // Access via IAxis
  BOOST_CHECK(epa.getAxis().isEquidistant());

  BOOST_CHECK(!epa.getAxis().isVariable());

  BOOST_CHECK_EQUAL(epa.getAxis().getType(), Equidistant);

  BOOST_CHECK_EQUAL(epa.getAxis().getBoundaryType(), Bound);

  BOOST_CHECK_EQUAL(epa.getAxis().getNBins(), 10);

  std::string rString =
      "ProtoAxis: 10 bins in AxisX, equidistant within automatic range";
  std::string oString = epa.toString();
  BOOST_CHECK_EQUAL(rString, oString);

  // Invalid constructor, closed with something else than phi or rphi
  BOOST_CHECK_THROW(Acts::ProtoAxis(AxisZ, Closed, 10), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(VariableProtoAxis) {
  using enum Acts::AxisBoundaryType;
  using enum Acts::AxisDirection;
  using enum Acts::AxisType;

  // Bound, equidistant axis
  Acts::ProtoAxis vpab(AxisX, Bound, {0.0, 1.0, 10});

  // Direct access
  BOOST_CHECK_EQUAL(vpab.getAxisDirection(), AxisX);

  BOOST_CHECK(!vpab.isAutorange());

  // Access via IAxis
  BOOST_CHECK(!vpab.getAxis().isEquidistant());

  BOOST_CHECK(vpab.getAxis().isVariable());

  auto edges = vpab.getAxis().getBinEdges();
  BOOST_CHECK_EQUAL(edges.size(), 3);

  BOOST_CHECK_EQUAL(vpab.getAxis().getType(), Variable);

  BOOST_CHECK_EQUAL(vpab.getAxis().getBoundaryType(), Bound);

  BOOST_CHECK_EQUAL(vpab.getAxis().getNBins(), 2);

  CHECK_CLOSE_ABS(vpab.getAxis().getMin(), 0.0, 1e-15);

  CHECK_CLOSE_ABS(vpab.getAxis().getMax(), 10.0, 1e-15);

  std::string rString = "ProtoAxis: 2 bins in AxisX, variable within [0, 10]";
  std::string oString = vpab.toString();
  BOOST_CHECK_EQUAL(rString, oString);

  // Invalid constructor, min > max
  BOOST_CHECK_THROW(Acts::ProtoAxis(AxisZ, Bound, std::vector<double>{2.}),
                    std::invalid_argument);

  // Invalid constructor, nbins < 1
  BOOST_CHECK_THROW(Acts::ProtoAxis(AxisZ, Bound, {3., 2., 1}),
                    std::invalid_argument);

  // Invalid constructor, closed with something else than phi or rphi
  std::vector<Acts::AxisDirection> invalidDirections = {
      AxisX, AxisY, AxisZ, AxisR, AxisEta, AxisTheta, AxisMag};
  for (const auto& adir : invalidDirections) {
    BOOST_CHECK_THROW(Acts::ProtoAxis(adir, Closed, {0., 1., 2., 3.}),
                      std::invalid_argument);
  }
}

BOOST_AUTO_TEST_SUITE_END()
