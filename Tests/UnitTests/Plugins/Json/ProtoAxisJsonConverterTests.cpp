// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Plugins/Json/ProtoAxisJsonConverter.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

BOOST_AUTO_TEST_SUITE(ProtoAxisJsonConversion)

BOOST_AUTO_TEST_CASE(EquidistantProtoAxisJsonConversion) {
  using enum Acts::AxisBoundaryType;
  using enum Acts::AxisDirection;
  using enum Acts::AxisType;

  // Bound, equidistant axis
  Acts::ProtoAxis epab(AxisX, Bound, 0.0, 1.0, 10);

  nlohmann::json jProtoAxis = Acts::ProtoAxisJsonConverter::toJson(epab);

  BOOST_CHECK(jProtoAxis.contains("axis"));
  BOOST_CHECK(jProtoAxis.contains("axis_dir"));
  BOOST_CHECK(jProtoAxis.contains("autorange"));

  Acts::ProtoAxis epabRead = Acts::ProtoAxisJsonConverter::fromJson(jProtoAxis);

  BOOST_CHECK_EQUAL(epabRead.getAxisDirection(), epab.getAxisDirection());
  BOOST_CHECK_EQUAL(epabRead.getAxis(), epab.getAxis());
  BOOST_CHECK_EQUAL(epabRead.isAutorange(), epab.isAutorange());
  BOOST_CHECK_EQUAL(epabRead.toString(), epab.toString());
}

BOOST_AUTO_TEST_CASE(AutorangeProtoAxisJsonConversion) {
  using enum Acts::AxisBoundaryType;
  using enum Acts::AxisDirection;
  using enum Acts::AxisType;

  // Bound, equidistant axis, autorange
  Acts::ProtoAxis epa(AxisX, Bound, 10);

  nlohmann::json jProtoAxis = Acts::ProtoAxisJsonConverter::toJson(epa);

  BOOST_CHECK(jProtoAxis.contains("axis"));
  BOOST_CHECK(jProtoAxis.contains("axis_dir"));
  BOOST_CHECK(jProtoAxis.contains("autorange"));

  Acts::ProtoAxis epaRead = Acts::ProtoAxisJsonConverter::fromJson(jProtoAxis);

  BOOST_CHECK_EQUAL(epaRead.getAxisDirection(), epa.getAxisDirection());
  BOOST_CHECK_EQUAL(epaRead.getAxis(), epa.getAxis());
  BOOST_CHECK_EQUAL(epaRead.isAutorange(), epa.isAutorange());
  BOOST_CHECK_EQUAL(epaRead.toString(), epa.toString());
}

BOOST_AUTO_TEST_CASE(VariableProtoAxisJsonConversion) {
  using enum Acts::AxisBoundaryType;
  using enum Acts::AxisDirection;
  using enum Acts::AxisType;

  // Bound, equidistant axis
  Acts::ProtoAxis vpab(AxisX, Bound, {0.0, 1.0, 10});

  nlohmann::json jProtoAxis = Acts::ProtoAxisJsonConverter::toJson(vpab);
  BOOST_CHECK(jProtoAxis.contains("axis"));
  BOOST_CHECK(jProtoAxis.contains("axis_dir"));
  BOOST_CHECK(jProtoAxis.contains("autorange"));

  Acts::ProtoAxis vpabRead = Acts::ProtoAxisJsonConverter::fromJson(jProtoAxis);

  BOOST_CHECK_EQUAL(vpabRead.getAxisDirection(), vpab.getAxisDirection());
  BOOST_CHECK_EQUAL(vpabRead.getAxis(), vpab.getAxis());
  BOOST_CHECK_EQUAL(vpabRead.isAutorange(), vpab.isAutorange());
  BOOST_CHECK_EQUAL(vpabRead.toString(), vpab.toString());
}

BOOST_AUTO_TEST_SUITE_END()
