// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/detail/GenerateParameters.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "ActsExamples/EventData/Measurement.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

using namespace Acts;
using namespace Acts::detail::Test;
using namespace ActsExamples;
namespace bd = boost::unit_test::data;

namespace {
constexpr BoundIndices boundIndices[] = {
    eBoundLoc0, eBoundLoc1, eBoundTime, eBoundPhi, eBoundTheta, eBoundQOverP,
};
constexpr Acts::GeometryIdentifier geoId{1};
// fix seed for reproducible tests
std::default_random_engine rng(123);
}  // namespace

// the underlying subspace implementation is already tested exhaustively in a
// separate unit test. here we only test concrete extreme cases and
// measurement-specific functionality.

BOOST_AUTO_TEST_SUITE(EventDataMeasurement)

BOOST_DATA_TEST_CASE(VariableBoundOne, bd::make(boundIndices), index) {
  MeasurementContainer container;

  auto [params, cov] = generateParametersCovariance<double, 1u>(rng);

  FixedBoundMeasurementProxy<1> meas = container.makeMeasurement<1>(geoId);
  meas.setSubspaceIndices(std::array{index});
  meas.parameters() = params;
  meas.covariance() = cov;

  BOOST_CHECK_EQUAL(meas.size(), 1);
  for (auto i : boundIndices) {
    BOOST_CHECK_EQUAL(meas.contains(i), i == index);
  }
  BOOST_CHECK_EQUAL(meas.parameters(), params);
  BOOST_CHECK_EQUAL(meas.covariance(), cov);
  BOOST_CHECK_EQUAL(meas.geometryId(), geoId);
}

BOOST_DATA_TEST_CASE(VariableBoundOneEmplace, bd::make(boundIndices), index) {
  MeasurementContainer container;

  auto [params, cov] = generateParametersCovariance<double, 1u>(rng);

  FixedBoundMeasurementProxy<1> meas =
      container.emplaceMeasurement<1>(geoId, std::array{index}, params, cov);

  BOOST_CHECK_EQUAL(meas.size(), 1);
  for (auto i : boundIndices) {
    BOOST_CHECK_EQUAL(meas.contains(i), i == index);
  }
  BOOST_CHECK_EQUAL(meas.parameters(), params);
  BOOST_CHECK_EQUAL(meas.covariance(), cov);
  BOOST_CHECK_EQUAL(meas.geometryId(), geoId);
}

BOOST_AUTO_TEST_CASE(VariableBoundAll) {
  MeasurementContainer container;

  auto [params, cov] = generateParametersCovariance<double, eBoundSize>(rng);

  FixedBoundMeasurementProxy<eBoundSize> meas =
      container.makeMeasurement<eBoundSize>(geoId);
  meas.setSubspaceIndices(std::array{eBoundLoc0, eBoundLoc1, eBoundTime,
                                     eBoundPhi, eBoundTheta, eBoundQOverP});
  meas.parameters() = params;
  meas.covariance() = cov;

  BOOST_CHECK_EQUAL(meas.size(), eBoundSize);
  for (auto i : boundIndices) {
    BOOST_CHECK(meas.contains(i));
  }
  BOOST_CHECK_EQUAL(meas.parameters(), params);
  BOOST_CHECK_EQUAL(meas.covariance(), cov);
  BOOST_CHECK_EQUAL(meas.geometryId(), geoId);
}

BOOST_AUTO_TEST_CASE(VariableBoundAllEmplace) {
  MeasurementContainer container;

  auto [params, cov] = generateParametersCovariance<double, eBoundSize>(rng);

  FixedBoundMeasurementProxy<eBoundSize> meas =
      container.emplaceMeasurement<eBoundSize>(
          geoId,
          std::array{eBoundLoc0, eBoundLoc1, eBoundTime, eBoundPhi, eBoundTheta,
                     eBoundQOverP},
          params, cov);

  BOOST_CHECK_EQUAL(meas.size(), eBoundSize);
  for (auto i : boundIndices) {
    BOOST_CHECK(meas.contains(i));
  }
  BOOST_CHECK_EQUAL(meas.parameters(), params);
  BOOST_CHECK_EQUAL(meas.covariance(), cov);
  BOOST_CHECK_EQUAL(meas.geometryId(), geoId);
}

BOOST_AUTO_TEST_CASE(VariableBoundReassign) {
  MeasurementContainer container;

  // generate w/ two parameter
  auto [params1, cov1] = generateParametersCovariance<double, 2u>(rng);

  VariableBoundMeasurementProxy meas = container.makeMeasurement(2, geoId);
  meas.setSubspaceIndices(std::array{eBoundPhi, eBoundTheta});
  meas.parameters() = params1;
  meas.covariance() = cov1;

  BOOST_CHECK_EQUAL(meas.size(), 2);
  BOOST_CHECK(!meas.contains(eBoundLoc0));
  BOOST_CHECK(!meas.contains(eBoundLoc1));
  BOOST_CHECK(!meas.contains(eBoundTime));
  BOOST_CHECK(meas.contains(eBoundPhi));
  BOOST_CHECK(meas.contains(eBoundTheta));
  BOOST_CHECK(!meas.contains(eBoundQOverP));

  // reassign w/ all parameters
  auto [paramsN, covN] = generateParametersCovariance<double, eBoundSize>(rng);

  meas = container.makeMeasurement(eBoundSize, geoId);
  meas.setSubspaceIndices(std::array{eBoundLoc0, eBoundLoc1, eBoundTime,
                                     eBoundPhi, eBoundTheta, eBoundQOverP});
  meas.parameters() = paramsN;
  meas.covariance() = covN;

  BOOST_CHECK_EQUAL(meas.size(), eBoundSize);
  BOOST_CHECK(meas.contains(eBoundLoc0));
  BOOST_CHECK(meas.contains(eBoundLoc1));
  BOOST_CHECK(meas.contains(eBoundTime));
  BOOST_CHECK(meas.contains(eBoundPhi));
  BOOST_CHECK(meas.contains(eBoundTheta));
  BOOST_CHECK(meas.contains(eBoundQOverP));
}

BOOST_AUTO_TEST_SUITE_END()
