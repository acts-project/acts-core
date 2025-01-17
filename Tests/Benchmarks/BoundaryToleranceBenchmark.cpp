// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/ConvexPolygonBounds.hpp"
#include "Acts/Tests/CommonHelpers/BenchmarkTools.hpp"

#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

using namespace Acts;

int main(int /*argc*/, char** /*argv[]*/) {
  // === PROBLEM DATA ===

  // Trapezoidal area of interest
  ConvexPolygonBounds<PolygonDynamic> poly(
      {{0.4, 0.25}, {0.6, 0.25}, {0.8, 0.75}, {0.2, 0.75}});

  // Covariance matrix which specifies "soft" boundary check tolerance
  SquareMatrix2 cov;
  cov << 0.2, 0.02, 0.15, 0.02;

  // Random tests cover the ((0, 0), (1, 1)) rectangle. 20% of that area is
  // covered by the trapezoid and most of it is covered by a 3*sigma tolerance
  // given the above covariance matrix.
  std::mt19937 rng(42);
  std::uniform_real_distribution<double> axis(0, 1);
  auto random_point = [&]() -> Vector2 {
    return Vector2(axis(rng), axis(rng));
  };

  // This point is inside the area
  const Vector2 center(0.5, 0.5);
  // This point is still inside the area, but close to an edge
  const Vector2 edge_inside(0.401, 0.251);
  // This point is just a bit outside, should be considered "in" by tolerance
  const Vector2 edge_outside(0.399, 0.249);
  // This point should always be considered outside the area
  const Vector2 far_away(-1000., -1000.);

  // === BENCHMARKS ===

  // Number of benchmark runs
  constexpr int NTESTS = 5'000;

  // Some checks are much slower, so we tune down benchmark iterations
  constexpr int NTESTS_SLOW = NTESTS / 10;

  // Conversely, no-op tests are so fast that we need to tune up iterations
  constexpr int NTESTS_NOOP = NTESTS * 10;

  // We use this to switch between iteration counts
  enum class Mode { None, FastOutside, SlowOutside };

  // Benchmark output display
  auto print_bench_header = [](const std::string& check_name) {
    std::cout << check_name << ":" << std::endl;
  };
  auto print_bench_result = [](const std::string& bench_name,
                               const Acts::Test::MicroBenchmarkResult& res) {
    std::cout << "- " << bench_name << ": " << res << std::endl;
  };

  // Benchmark runner
  auto run_bench = [&](auto&& iteration, int num_iters,
                       const std::string& bench_name) {
    auto bench_result = Acts::Test::microBenchmark(iteration, num_iters);
    print_bench_result(bench_name, bench_result);
  };
  auto run_bench_with_inputs = [&](auto&& iterationWithArg, auto&& inputs,
                                   const std::string& bench_name) {
    auto bench_result = Acts::Test::microBenchmark(iterationWithArg, inputs);
    print_bench_result(bench_name, bench_result);
  };
  auto run_all_benches = [&](const BoundaryTolerance& check,
                             const std::string& check_name, const Mode mode) {
    // Announce a set of benchmarks
    print_bench_header(check_name);

    // Pre-determined "interesting" test points
    int num_inside_points = 0;
    int num_outside_points = 0;
    switch (mode) {
      case Mode::None:
        num_inside_points = NTESTS_NOOP;
        num_outside_points = NTESTS_NOOP;
        break;
      case Mode::FastOutside:
        num_inside_points = NTESTS;
        num_outside_points = NTESTS;
        break;
      case Mode::SlowOutside:
        num_inside_points = NTESTS;
        num_outside_points = NTESTS_SLOW;
      default:  // do nothing
        break;
    };
    run_bench([&] { return poly.inside(center, check); }, num_inside_points,
              "Center");
    run_bench([&] { return poly.inside(edge_inside, check); },
              num_inside_points, "Inside edge");
    run_bench([&] { return poly.inside(edge_outside, check); },
              num_outside_points, "Outside edge");
    run_bench([&] { return poly.inside(far_away, check); }, num_outside_points,
              "Far away");

    // Pre-rolled random points
    std::vector<Vector2> points(num_outside_points);
    std::generate(points.begin(), points.end(), random_point);
    run_bench_with_inputs(
        [&](const auto& point) { return poly.inside(point, check); }, points,
        "Random");
  };

  // Benchmark scenarios
  run_all_benches(BoundaryTolerance::Infinite(), "No check", Mode::None);
  run_all_benches(BoundaryTolerance::None(), "No tolerance", Mode::FastOutside);
  run_all_benches(BoundaryTolerance::AbsoluteBound(0.6, 0.45), "Abs. tolerance",
                  Mode::SlowOutside);
  run_all_benches(BoundaryTolerance::Chi2Bound(cov, 3.0), "Cov. tolerance",
                  Mode::SlowOutside);

  return 0;
}
