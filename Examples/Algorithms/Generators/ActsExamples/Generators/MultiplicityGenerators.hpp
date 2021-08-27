// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2018-03-13
/// @author Moritz Kiehn <msmk@cern.ch>

#pragma once

#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Generators/EventGenerator.hpp"

#include <random>

namespace ActsExamples {

struct FixedMultiplicityGenerator
    : public EventGenerator::MultiplicityGenerator {
  size_t n = 1;

  FixedMultiplicityGenerator(size_t _n) : n{_n} {}
  FixedMultiplicityGenerator() = default;

  size_t operator()(RandomEngine& /* unused */) const override { return n; }
};

struct PoissonMultiplicityGenerator
    : public EventGenerator::MultiplicityGenerator {
  double mean = 1;
  PoissonMultiplicityGenerator(double _mean) : mean{_mean} {}
  PoissonMultiplicityGenerator() = default;

  size_t operator()(RandomEngine& rng) const override {
    return (0 < mean) ? std::poisson_distribution<size_t>(mean)(rng) : 0;
  }
};

}  // namespace ActsExamples
