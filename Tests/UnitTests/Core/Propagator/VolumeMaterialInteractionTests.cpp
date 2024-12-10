// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/detail/VolumeMaterialInteraction.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"

#include <memory>

using namespace Acts::UnitLiterals;

namespace Acts::Test {

/// @brief Simplified stepper state
struct StepperState {
  ParticleHypothesis particleHypothesis = ParticleHypothesis::pion();
  Vector3 pos, dir;
  double t = 0, p = 0, q = 0;
  bool covTransport = false;
  double absCharge = UnitConstants::e;
};

/// @brief Simplified navigator
struct NaivgatorState {
  TrackingVolume* currentVolume = nullptr;
};

/// @brief Simplified propagator state
struct State {
  struct {
    Direction direction = Direction::Forward();
  } options;

  StepperState stepping;
  NaivgatorState navigation;
};

/// @brief Simplified stepper
struct Stepper {
  Stepper() = default;

  Vector3 position(const StepperState& state) const { return state.pos; }

  double time(const StepperState& state) const { return state.t; }

  Vector3 direction(const StepperState& state) const { return state.dir; }

  double qOverP(const StepperState& state) const { return state.q / state.p; }

  double absoluteMomentum(const StepperState& state) const { return state.p; }

  double charge(const StepperState& state) const { return state.q; };

  ParticleHypothesis particleHypothesis(const StepperState& state) const {
    return state.particleHypothesis;
  };
};

/// @brief Simplified navigator
struct Navigator {
  const TrackingVolume* currentVolume(const NaivgatorState& state) const {
    return state.currentVolume;
  }
};

BOOST_AUTO_TEST_CASE(volume_material_interaction_test) {
  // Create a Tracking Volume
  auto htrans = Transform3(Translation3{-10., -10., 0.});
  auto bound = std::make_shared<CuboidVolumeBounds>(1_m, 1_m, 1_m);
  auto mat = makeSilicon();
  auto volMat = std::make_shared<const HomogeneousVolumeMaterial>(mat);
  auto volume = std::make_shared<TrackingVolume>(htrans, bound, volMat);

  // Create a propagator state
  State state;
  state.stepping.particleHypothesis =
      ParticleHypothesis(static_cast<PdgParticle>(11), 10., AnyCharge{9.});
  state.stepping.pos = Vector3(1., 2., 3.);
  state.stepping.dir = Vector3(4., 5., 6.);
  state.stepping.t = 7.;
  state.stepping.p = 8.;
  state.stepping.q = 9.;
  state.stepping.absCharge = std::abs(state.stepping.q);
  state.stepping.covTransport = true;
  state.options.direction = Direction::Backward();
  state.navigation.currentVolume = volume.get();

  Stepper stepper;
  Navigator navigator;

  // Build the VolumeMaterialInteraction & test assignments
  detail::VolumeMaterialInteraction volMatInt(volume.get(), state, stepper);
  BOOST_CHECK_EQUAL(volMatInt.volume.trackingVolume, volume.get());
  BOOST_CHECK_EQUAL(volMatInt.pos, stepper.position(state.stepping));
  BOOST_CHECK_EQUAL(volMatInt.time, stepper.time(state.stepping));
  BOOST_CHECK_EQUAL(volMatInt.dir, stepper.direction(state.stepping));
  BOOST_CHECK_EQUAL(volMatInt.momentum,
                    stepper.absoluteMomentum(state.stepping));
  BOOST_CHECK_EQUAL(volMatInt.absQ, std::abs(stepper.charge(state.stepping)));
  CHECK_CLOSE_ABS(volMatInt.qOverP, stepper.qOverP(state.stepping), 1e-6);
  BOOST_CHECK_EQUAL(volMatInt.mass, state.stepping.particleHypothesis.mass());
  BOOST_CHECK_EQUAL(volMatInt.absPdg,
                    state.stepping.particleHypothesis.absolutePdg());
  BOOST_CHECK_EQUAL(volMatInt.performCovarianceTransport,
                    state.stepping.covTransport);
  BOOST_CHECK_EQUAL(volMatInt.navDir, state.options.direction);

  // Evaluate the material
  bool result = volMatInt.evaluateMaterialSlab(state, navigator);
  BOOST_CHECK(result);
  BOOST_CHECK_EQUAL(volMatInt.slab.material(), mat);
  BOOST_CHECK_EQUAL(volMatInt.slab.thickness(), 1.);
  BOOST_CHECK_EQUAL(volMatInt.pathCorrection, 0.);

  // Evaluate the material without a tracking volume
  state.navigation.currentVolume = nullptr;
  result = volMatInt.evaluateMaterialSlab(state, navigator);
  BOOST_CHECK(!result);
  BOOST_CHECK_EQUAL(volMatInt.pathCorrection, 0.);
}

}  // namespace Acts::Test
