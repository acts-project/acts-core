// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Tests/CommonHelpers/BenchmarkTools.hpp"
#include "Acts/Utilities/QuickMath.hpp"

#include <cmath>
#include <random>

namespace bdata = boost::unit_test::data;

namespace Acts::Test {

/// @brief Another fast power function @see `fastPow`
// Taken from
// https://martin.ankerl.com/2007/02/11/optimized-exponential-functions-for-java
/// @param a the base
/// @param b the exponent
constexpr double fastPowAnother(double a, double b) {
  // enable only on IEEE 754
  static_assert(std::numeric_limits<double>::is_iec559);

  union {
    double f;
    std::int64_t i;
  } u = {};

  u.i = static_cast<std::int64_t>(
      9076650 * (a - 1) / (a + 1 + 4 * std::sqrt(a)) * b + 1072632447);
  u.i <<= 32;

  // result seems broken?
  return u.f;
}

// Some randomness & number crunching
unsigned int ntests = 10;
unsigned int nrepts = 10000;

BOOST_DATA_TEST_CASE(
    benchmark_pow_25,
    bdata::random(
        (bdata::engine = std::mt19937(), bdata::seed = 21,
         bdata::distribution = std::uniform_real_distribution<double>(-4, 4))) ^
        bdata::xrange(ntests),
    aExp, index) {
  (void)index;

  double a = std::pow(10, aExp);
  double b = 0.25;

  std::cout << std::endl
            << "Benchmarking a=" << a << ", b=" << b << "..." << std::endl;
  std::cout << "- void: "
            << Acts::Test::microBenchmark([&] { return 0; }, nrepts)
            << std::endl;
  std::cout << "- std::pow: "
            << Acts::Test::microBenchmark([&] { return std::pow(a, b); },
                                          nrepts)
            << std::endl;
  std::cout << "- std::exp: "
            << Acts::Test::microBenchmark(
                   [&] { return std::exp(std::log(a) * b); }, nrepts)
            << std::endl;
  std::cout << "- std::sqrt: "
            << Acts::Test::microBenchmark(
                   [&] { return std::sqrt(std::sqrt(a)); }, nrepts)
            << std::endl;
  std::cout << "- fastPow: "
            << Acts::Test::microBenchmark([&] { return fastPow(a, b); }, nrepts)
            << std::endl;
  std::cout << "- fastPowMorePrecise: "
            << Acts::Test::microBenchmark(
                   [&] { return fastPowMorePrecise(a, b); }, nrepts)
            << std::endl;
  std::cout << "- fastPowAnother: "
            << Acts::Test::microBenchmark([&] { return fastPowAnother(a, b); },
                                          nrepts)
            << std::endl;
}

}  // namespace Acts::Test
