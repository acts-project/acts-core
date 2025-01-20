// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Tests/CommonHelpers/BenchmarkTools.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <iostream>

#include <boost/program_options.hpp>

namespace Acts::Test {

namespace po = boost::program_options;
using namespace Acts;
using namespace Acts::UnitLiterals;

struct BenchmarkStepper {
  unsigned int toys{};
  double ptInGeV{};
  double BzInT{};
  double maxPathInM{};
  unsigned int lvl{};
  bool withCov{};

  std::optional<int> parseOptions(int argc, char* argv[]) {
    try {
      po::options_description desc("Allowed options");
      // clang-format off
      desc.add_options()
        ("help", "produce help message")
        ("toys",po::value<unsigned int>(&toys)->default_value(20000),"number of tracks to propagate")
        ("pT",po::value<double>(&ptInGeV)->default_value(1),"transverse momentum in GeV")
        ("B",po::value<double>(&BzInT)->default_value(2),"z-component of B-field in T")
        ("path",po::value<double>(&maxPathInM)->default_value(5),"maximum path length in m")
        ("cov",po::value<bool>(&withCov)->default_value(true),"propagation with covariance matrix")
        ("verbose",po::value<unsigned int>(&lvl)->default_value(Acts::Logging::INFO),"logging level");
      // clang-format on
      po::variables_map vm;
      po::store(po::parse_command_line(argc, argv, desc), vm);
      po::notify(vm);

      if (vm.contains("help")) {
        std::cout << desc << std::endl;
        return 0;
      }
    } catch (std::exception& e) {
      std::cerr << "error: " << e.what() << std::endl;
      return 1;
    }

    return std::nullopt;
  }

  std::unique_ptr<MagneticFieldProvider> makeField() const {
    return std::make_unique<ConstantBField>(
        Vector3{0, 0, BzInT * UnitConstants::T});
  }

  template <typename Stepper>
  void run(Stepper stepper, const std::string& name) const {
    using Propagator = Propagator<Stepper>;
    using PropagatorOptions = typename Propagator::template Options<>;
    using Covariance = BoundSquareMatrix;

    // Create a test context
    GeometryContext tgContext = GeometryContext();
    MagneticFieldContext mfContext = MagneticFieldContext();

    ACTS_LOCAL_LOGGER(getDefaultLogger(name, Acts::Logging::Level(lvl)));

    // print information about profiling setup
    ACTS_INFO("propagating " << toys << " tracks with pT = " << ptInGeV
                             << "GeV in a " << BzInT << "T B-field");

    Propagator propagator(std::move(stepper));

    PropagatorOptions options(tgContext, mfContext);
    options.pathLimit = maxPathInM * UnitConstants::m;

    Vector4 pos4(0, 0, 0, 0);
    Vector3 dir(1, 0, 0);
    Covariance cov;
    // clang-format off
    cov << 10_mm, 0, 0, 0, 0, 0,
            0, 10_mm, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0,
            0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1_e / 10_GeV, 0,
            0, 0, 0, 0, 0, 0;
    // clang-format on

    std::optional<Covariance> covOpt = std::nullopt;
    if (withCov) {
      covOpt = cov;
    }
    CurvilinearTrackParameters pars(pos4, dir, +1 / ptInGeV, covOpt,
                                    ParticleHypothesis::pion());

    double totalPathLength = 0;
    std::size_t numSteps = 0;
    std::size_t numStepTrials = 0;
    std::size_t numIters = 0;
    const auto propagationBenchResult = Acts::Test::microBenchmark(
        [&] {
          auto state = propagator.makeState(pars, options);
          propagator.initialize(state);
          auto tmp = propagator.propagate(state);
          auto r = propagator.makeResult(state, tmp, options, true).value();
          if (totalPathLength == 0.) {
            ACTS_DEBUG("reached position "
                       << r.endParameters->position(tgContext).transpose()
                       << " in " << r.steps << " steps");
          }
          totalPathLength += r.pathLength;
          numSteps += r.steps;
          numStepTrials += state.stepping.nStepTrials;
          ++numIters;
          return r;
        },
        1, toys);

    ACTS_INFO("Execution stats: " << propagationBenchResult);
    ACTS_INFO("average path length = " << totalPathLength / numIters / 1_mm
                                       << "mm");
    ACTS_INFO("average number of steps = " << 1.0 * numSteps / numIters);
    ACTS_INFO("step efficiency = " << 1.0 * numSteps / numStepTrials);
  }
};

}  // namespace Acts::Test
