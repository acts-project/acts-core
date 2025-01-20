// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/GenericCurvilinearTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/Assertions.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Result.hpp"

#include <algorithm>
#include <array>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <type_traits>
#include <utility>

using namespace Acts::UnitLiterals;

namespace Acts::Test {

using Acts::VectorHelpers::makeVector4;
using Covariance = BoundSquareMatrix;
using Jacobian = BoundMatrix;
using Stepper = Acts::AtlasStepper;

// epsilon for floating point comparisons
static constexpr auto eps = 1024 * std::numeric_limits<double>::epsilon();

// propagation settings
static constexpr auto stepSize = 10_mm;
static constexpr Direction navDir = Direction::Backward();
static auto magneticField =
    std::make_shared<ConstantBField>(Vector3(0.1_T, -0.2_T, 2_T));

// initial parameter state
static const Vector4 pos4(1_mm, -1_mm, 2_mm, 2_ns);
static const Vector3 pos = pos4.segment<3>(ePos0);
static const auto time = pos4[eTime];
static const Vector3 unitDir = Vector3(-2, 2, 1).normalized();
static constexpr auto absMom = 1_GeV;
static constexpr auto charge = -1_e;
static const auto particleHypothesis = ParticleHypothesis::pion();
static const Covariance cov = Covariance::Identity();

// context objects
static const GeometryContext geoCtx;
static const MagneticFieldContext magCtx;

BOOST_AUTO_TEST_SUITE(AtlasStepper)

// test state construction from parameters w/o covariance
BOOST_AUTO_TEST_CASE(ConstructState) {
  CurvilinearTrackParameters cp(pos4, unitDir, charge / absMom, std::nullopt,
                                particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options, cp);

  BOOST_CHECK(!state.covTransport);
  BOOST_CHECK_EQUAL(state.covariance, nullptr);
  BOOST_CHECK_EQUAL(state.pVector[0], pos.x());
  BOOST_CHECK_EQUAL(state.pVector[1], pos.y());
  BOOST_CHECK_EQUAL(state.pVector[2], pos.z());
  BOOST_CHECK_EQUAL(state.pVector[3], time);
  CHECK_CLOSE_ABS(state.pVector[4], unitDir.x(), eps);
  CHECK_CLOSE_ABS(state.pVector[5], unitDir.y(), eps);
  CHECK_CLOSE_ABS(state.pVector[6], unitDir.z(), eps);
  BOOST_CHECK_EQUAL(state.pVector[7], charge / absMom);
  BOOST_CHECK_EQUAL(state.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(state.stepSize.value(), stepSize);
  BOOST_CHECK_EQUAL(state.previousStepSize, 0.);
}

// test state construction from parameters w/ covariance
BOOST_AUTO_TEST_CASE(ConstructStateWithCovariance) {
  CurvilinearTrackParameters cp(pos4, unitDir, charge / absMom, cov,
                                particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options, cp);

  BOOST_CHECK(state.covTransport);
  BOOST_CHECK_EQUAL(*state.covariance, cov);
  BOOST_CHECK_EQUAL(state.pVector[0], pos.x());
  BOOST_CHECK_EQUAL(state.pVector[1], pos.y());
  BOOST_CHECK_EQUAL(state.pVector[2], pos.z());
  BOOST_CHECK_EQUAL(state.pVector[3], time);
  CHECK_CLOSE_ABS(state.pVector[4], unitDir.x(), eps);
  CHECK_CLOSE_ABS(state.pVector[5], unitDir.y(), eps);
  CHECK_CLOSE_ABS(state.pVector[6], unitDir.z(), eps);
  BOOST_CHECK_EQUAL(state.pVector[7], charge / absMom);
  BOOST_CHECK_EQUAL(state.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(state.stepSize.value(), stepSize);
  BOOST_CHECK_EQUAL(state.previousStepSize, 0.);
}

// test stepper getters for particle state
BOOST_AUTO_TEST_CASE(Getters) {
  CurvilinearTrackParameters cp(pos4, unitDir, charge / absMom, cov,
                                particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options, cp);

  CHECK_CLOSE_ABS(stepper.position(state), pos, eps);
  CHECK_CLOSE_ABS(stepper.time(state), time, eps);
  CHECK_CLOSE_ABS(stepper.direction(state), unitDir, eps);
  CHECK_CLOSE_ABS(stepper.absoluteMomentum(state), absMom, eps);
  BOOST_CHECK_EQUAL(stepper.charge(state), charge);
}

// test stepper update methods with bound state as input
BOOST_AUTO_TEST_CASE(UpdateFromBound) {
  CurvilinearTrackParameters cp(pos4, unitDir, charge / absMom, cov,
                                particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options, cp);

  auto newPos4 = (pos4 + Vector4(1_mm, 2_mm, 3_mm, 20_ns)).eval();
  auto newPos = newPos4.segment<3>(ePos0);
  auto newTime = newPos4[eTime];
  auto newUnitDir = (unitDir + Vector3(1, -1, -1)).normalized();
  auto newAbsMom = 0.9 * absMom;

  // example surface and bound parameters at the updated position
  auto plane = CurvilinearSurface(newPos, newUnitDir).planeSurface();
  auto params =
      BoundTrackParameters::create(plane, geoCtx, newPos4, newUnitDir,
                                   charge / absMom, cov, particleHypothesis)
          .value();
  FreeVector freeParams;
  freeParams[eFreePos0] = newPos4[ePos0];
  freeParams[eFreePos1] = newPos4[ePos1];
  freeParams[eFreePos2] = newPos4[ePos2];
  freeParams[eFreeTime] = newPos4[eTime];
  freeParams[eFreeDir0] = newUnitDir[eMom0];
  freeParams[eFreeDir1] = newUnitDir[eMom1];
  freeParams[eFreeDir2] = newUnitDir[eMom2];
  freeParams[eFreeQOverP] = charge / newAbsMom;

  // WARNING for some reason there seems to be an additional flag that makes
  //         the update method not do anything when it is set. Why?
  state.state_ready = false;
  BOOST_CHECK(params.covariance().has_value());
  stepper.update(state, freeParams, params.parameters(), *params.covariance(),
                 *plane);
  CHECK_CLOSE_ABS(stepper.position(state), newPos, eps);
  CHECK_CLOSE_ABS(stepper.time(state), newTime, eps);
  CHECK_CLOSE_ABS(stepper.direction(state), newUnitDir, eps);
  CHECK_CLOSE_ABS(stepper.absoluteMomentum(state), newAbsMom, eps);
  BOOST_CHECK_EQUAL(stepper.charge(state), charge);
}

// test stepper update methods with individual components as input
BOOST_AUTO_TEST_CASE(UpdateFromComponents) {
  CurvilinearTrackParameters cp(pos4, unitDir, charge / absMom, cov,
                                particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options, cp);

  auto newPos = (pos + Vector3(1_mm, 2_mm, 3_mm)).eval();
  auto newTime = time + 20_ns;
  auto newUnitDir = (unitDir + Vector3(1, -1, -1)).normalized();
  auto newAbsMom = 0.9 * absMom;

  stepper.update(state, newPos, newUnitDir, charge / newAbsMom, newTime);
  CHECK_CLOSE_ABS(stepper.position(state), newPos, eps);
  CHECK_CLOSE_ABS(stepper.time(state), newTime, eps);
  CHECK_CLOSE_ABS(stepper.direction(state), newUnitDir, eps);
  CHECK_CLOSE_ABS(stepper.absoluteMomentum(state), newAbsMom, eps);
  BOOST_CHECK_EQUAL(stepper.charge(state), charge);
}

// test building a bound state object from the stepper state
BOOST_AUTO_TEST_CASE(BuildBound) {
  CurvilinearTrackParameters cp(pos4, unitDir, charge / absMom, cov,
                                particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options, cp);

  // example surface at the current state position
  auto plane = CurvilinearSurface(pos, unitDir).planeSurface();

  auto&& [pars, jac, pathLength] = stepper.boundState(state, *plane).value();
  // check parameters
  CHECK_CLOSE_ABS(pars.position(geoCtx), pos, eps);
  CHECK_CLOSE_ABS(pars.time(), time, eps);
  CHECK_CLOSE_ABS(pars.momentum(), absMom * unitDir, eps);
  BOOST_CHECK_EQUAL(pars.charge(), charge);
  BOOST_CHECK(pars.covariance().has_value());
  BOOST_CHECK_NE(*pars.covariance(), cov);
  // check Jacobian. should be identity since there was no propagation yet
  CHECK_CLOSE_ABS(jac, Jacobian(Jacobian::Identity()), eps);
  // check propagation length
  CHECK_CLOSE_ABS(pathLength, 0., eps);
}

// test building a curvilinear state object from the stepper state
BOOST_AUTO_TEST_CASE(BuildCurvilinear) {
  CurvilinearTrackParameters cp(pos4, unitDir, charge / absMom, cov,
                                particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options, cp);

  auto&& [pars, jac, pathLength] = stepper.curvilinearState(state);
  // check parameters
  CHECK_CLOSE_ABS(pars.position(geoCtx), pos, eps);
  CHECK_CLOSE_ABS(pars.time(), time, eps);
  CHECK_CLOSE_ABS(pars.momentum(), absMom * unitDir, eps);
  BOOST_CHECK_EQUAL(pars.charge(), charge);
  BOOST_CHECK(pars.covariance().has_value());
  BOOST_CHECK_NE(*pars.covariance(), cov);
  // check Jacobian. should be identity since there was no propagation yet
  CHECK_CLOSE_ABS(jac, Jacobian(Jacobian::Identity()), eps);
  // check propagation length
  CHECK_CLOSE_ABS(pathLength, 0., eps);
}

// test step method without covariance transport
BOOST_AUTO_TEST_CASE(Step) {
  CurvilinearTrackParameters cp(pos4, unitDir, charge / absMom, cov,
                                particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  auto state = stepper.makeState(options, cp);
  state.covTransport = false;

  // ensure step does not result in an error
  auto res = stepper.step(state, Direction::Backward(), nullptr);
  BOOST_CHECK(res.ok());

  // extract the actual step size
  auto h = res.value();
  BOOST_CHECK_EQUAL(state.stepSize.value(), stepSize);
  BOOST_CHECK_EQUAL(state.stepSize.value(), h * navDir);

  // check that the position has moved
  auto deltaPos = (stepper.position(state) - pos).eval();
  BOOST_CHECK_LT(0, deltaPos.norm());
  // check that time has changed
  auto deltaTime = stepper.time(state) - time;
  BOOST_CHECK_LT(0, std::abs(deltaTime));
  // check that the direction was rotated
  auto projDir = stepper.direction(state).dot(unitDir);
  BOOST_CHECK_LT(projDir, 1);

  // momentum and charge should not change
  CHECK_CLOSE_ABS(stepper.absoluteMomentum(state), absMom, eps);
  BOOST_CHECK_EQUAL(stepper.charge(state), charge);
}

// test step method with covariance transport
BOOST_AUTO_TEST_CASE(StepWithCovariance) {
  CurvilinearTrackParameters cp(pos4, unitDir, charge / absMom, cov,
                                particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  auto state = stepper.makeState(options, cp);
  state.covTransport = true;

  // ensure step does not result in an error
  auto res = stepper.step(state, Direction::Backward(), nullptr);
  BOOST_CHECK(res.ok());

  // extract the actual step size
  auto h = res.value();
  BOOST_CHECK_EQUAL(state.stepSize.value(), stepSize);
  BOOST_CHECK_EQUAL(state.stepSize.value(), h * navDir);

  // check that the position has moved
  auto deltaPos = (stepper.position(state) - pos).eval();
  BOOST_CHECK_LT(0, deltaPos.norm());
  // check that time has changed
  auto deltaTime = stepper.time(state) - time;
  BOOST_CHECK_LT(0, std::abs(deltaTime));
  // check that the direction was rotated
  auto projDir = stepper.direction(state).dot(unitDir);
  BOOST_CHECK_LT(projDir, 1);

  // momentum and charge should not change
  CHECK_CLOSE_ABS(stepper.absoluteMomentum(state), absMom, eps);
  BOOST_CHECK_EQUAL(stepper.charge(state), charge);

  stepper.transportCovarianceToCurvilinear(state);
  BOOST_CHECK_NE(state.cov, cov);
}

// test state reset method
BOOST_AUTO_TEST_CASE(Reset) {
  CurvilinearTrackParameters cp(pos4, unitDir, charge / absMom, cov,
                                particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  auto state = stepper.makeState(options, cp);
  state.covTransport = true;

  // ensure step does not result in an error
  stepper.step(state, Direction::Backward(), nullptr);

  // Construct the parameters
  Vector3 newPos(1.5, -2.5, 3.5);
  auto newAbsMom = 4.2 * absMom;
  double newTime = 7.5;
  double newCharge = 1.;
  BoundSquareMatrix newCov = 8.5 * Covariance::Identity();
  cp = CurvilinearTrackParameters(makeVector4(newPos, newTime), unitDir,
                                  newCharge / newAbsMom, newCov,
                                  particleHypothesis);
  FreeVector freeParams = transformBoundToFreeParameters(
      cp.referenceSurface(), geoCtx, cp.parameters());
  Direction navDir = Direction::Forward();
  double stepSize = -256.;

  auto copyState = [&](auto& field, const auto& other) {
    using field_t = std::decay_t<decltype(field)>;
    std::decay_t<decltype(other)> copy = stepper.makeState(options, cp);

    copy.state_ready = other.state_ready;
    copy.useJacobian = other.useJacobian;
    copy.step = other.step;
    copy.maxPathLength = other.maxPathLength;
    copy.mcondition = other.mcondition;
    copy.needgradient = other.needgradient;
    copy.newfield = other.newfield;
    copy.field = other.field;
    copy.pVector = other.pVector;
    std::copy(std::begin(other.parameters), std::end(other.parameters),
              std::begin(copy.parameters));
    copy.covariance = other.covariance;
    copy.covTransport = other.covTransport;
    std::copy(std::begin(other.jacobian), std::end(other.jacobian),
              std::begin(copy.jacobian));
    copy.pathAccumulated = other.pathAccumulated;
    copy.stepSize = other.stepSize;
    copy.previousStepSize = other.previousStepSize;

    copy.fieldCache = MagneticFieldProvider::Cache(
        std::in_place_type<typename field_t::Cache>,
        other.fieldCache.template as<typename field_t::Cache>());

    copy.debug = other.debug;
    copy.debugString = other.debugString;
    copy.debugPfxWidth = other.debugPfxWidth;
    copy.debugMsgWidth = other.debugMsgWidth;

    return copy;
  };

  // Reset all possible parameters
  Stepper::State stateCopy = copyState(*magneticField, state);
  BOOST_CHECK(cp.covariance().has_value());
  stepper.resetState(stateCopy, cp.parameters(), *cp.covariance(),
                     cp.referenceSurface(), stepSize);
  // Test all components
  BOOST_CHECK(stateCopy.covTransport);
  BOOST_CHECK_EQUAL(*stateCopy.covariance, newCov);
  BOOST_CHECK_EQUAL(stepper.position(stateCopy),
                    freeParams.template segment<3>(eFreePos0));
  BOOST_CHECK_EQUAL(stepper.direction(stateCopy),
                    freeParams.template segment<3>(eFreeDir0).normalized());
  BOOST_CHECK_EQUAL(stepper.absoluteMomentum(stateCopy),
                    std::abs(1. / freeParams[eFreeQOverP]));
  BOOST_CHECK_EQUAL(stepper.charge(stateCopy), stepper.charge(state));
  BOOST_CHECK_EQUAL(stepper.time(stateCopy), freeParams[eFreeTime]);
  BOOST_CHECK_EQUAL(stateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(stateCopy.stepSize.value(), navDir * stepSize);
  BOOST_CHECK_EQUAL(stateCopy.previousStepSize, state.previousStepSize);

  // Reset all possible parameters except the step size
  stateCopy = copyState(*magneticField, state);
  stepper.resetState(stateCopy, cp.parameters(), *cp.covariance(),
                     cp.referenceSurface());
  // Test all components
  BOOST_CHECK(stateCopy.covTransport);
  BOOST_CHECK_EQUAL(*stateCopy.covariance, newCov);
  BOOST_CHECK_EQUAL(stepper.position(stateCopy),
                    freeParams.template segment<3>(eFreePos0));
  BOOST_CHECK_EQUAL(stepper.direction(stateCopy),
                    freeParams.template segment<3>(eFreeDir0).normalized());
  BOOST_CHECK_EQUAL(stepper.absoluteMomentum(stateCopy),
                    std::abs(1. / freeParams[eFreeQOverP]));
  BOOST_CHECK_EQUAL(stepper.charge(stateCopy), stepper.charge(state));
  BOOST_CHECK_EQUAL(stepper.time(stateCopy), freeParams[eFreeTime]);
  BOOST_CHECK_EQUAL(stateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(stateCopy.stepSize.value(),
                    std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(stateCopy.previousStepSize, state.previousStepSize);

  // Reset the least amount of parameters
  stateCopy = copyState(*magneticField, state);
  stepper.resetState(stateCopy, cp.parameters(), *cp.covariance(),
                     cp.referenceSurface());
  // Test all components
  BOOST_CHECK(stateCopy.covTransport);
  BOOST_CHECK_EQUAL(*stateCopy.covariance, newCov);
  BOOST_CHECK_EQUAL(stepper.position(stateCopy),
                    freeParams.template segment<3>(eFreePos0));
  BOOST_CHECK_EQUAL(stepper.direction(stateCopy),
                    freeParams.template segment<3>(eFreeDir0).normalized());
  BOOST_CHECK_EQUAL(stepper.absoluteMomentum(stateCopy),
                    std::abs(1. / freeParams[eFreeQOverP]));
  BOOST_CHECK_EQUAL(stepper.charge(stateCopy), stepper.charge(state));
  BOOST_CHECK_EQUAL(stepper.time(stateCopy), freeParams[eFreeTime]);
  BOOST_CHECK_EQUAL(stateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(stateCopy.stepSize.value(),
                    std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(stateCopy.previousStepSize, state.previousStepSize);

  // Reset using different surface shapes
  // 1) Disc surface
  // Setting some parameters
  newPos << 0.5, -1.5, 0.;
  newAbsMom *= 1.23;
  newTime = 8.4;
  newCharge = -1.;
  newCov = 10.9 * Covariance::Identity();
  Transform3 trafo = Transform3::Identity();
  auto disc = Surface::makeShared<DiscSurface>(trafo);
  auto boundDisc = BoundTrackParameters::create(
                       disc, geoCtx, makeVector4(newPos, newTime), unitDir,
                       newCharge / newAbsMom, newCov, particleHypothesis)
                       .value();

  // Reset the state and test
  Stepper::State stateDisc = copyState(*magneticField, state);
  BOOST_CHECK(boundDisc.covariance().has_value());
  stepper.resetState(stateDisc, boundDisc.parameters(), *boundDisc.covariance(),
                     boundDisc.referenceSurface());

  CHECK_NE_COLLECTIONS(stateDisc.pVector, stateCopy.pVector);
  CHECK_NE_COLLECTIONS(stateDisc.pVector, state.pVector);

  // 2) Perigee surface
  // Setting some parameters
  newPos << -2.06155, -2.06155, 3.5;
  newAbsMom *= 0.45;
  newTime = 2.3;
  newCharge = 1.;
  newCov = 8.7 * Covariance::Identity();
  auto perigee = Surface::makeShared<PerigeeSurface>(trafo);
  auto boundPerigee =
      BoundTrackParameters::create(
          perigee, geoCtx, makeVector4(newPos, newTime), unitDir,
          newCharge / newAbsMom, newCov, particleHypothesis)
          .value();

  // Reset the state and test
  Stepper::State statePerigee = copyState(*magneticField, state);
  BOOST_CHECK(boundPerigee.covariance().has_value());
  stepper.resetState(statePerigee, boundPerigee.parameters(),
                     *boundPerigee.covariance(),
                     boundPerigee.referenceSurface());
  CHECK_NE_COLLECTIONS(statePerigee.pVector, stateCopy.pVector);
  CHECK_NE_COLLECTIONS(statePerigee.pVector, state.pVector);
  CHECK_NE_COLLECTIONS(statePerigee.pVector, stateDisc.pVector);

  // 3) Straw surface
  // Use the same parameters as for previous Perigee surface
  auto straw = Surface::makeShared<StrawSurface>(trafo);
  auto boundStraw = BoundTrackParameters::create(
                        straw, geoCtx, makeVector4(newPos, newTime), unitDir,
                        newCharge / newAbsMom, newCov, particleHypothesis)
                        .value();

  // Reset the state and test
  Stepper::State stateStraw = copyState(*magneticField, state);
  BOOST_CHECK(boundStraw.covariance().has_value());
  stepper.resetState(stateStraw, boundStraw.parameters(),
                     *boundStraw.covariance(), boundStraw.referenceSurface());
  CHECK_NE_COLLECTIONS(stateStraw.pVector, stateCopy.pVector);
  CHECK_NE_COLLECTIONS(stateStraw.pVector, state.pVector);
  CHECK_NE_COLLECTIONS(stateStraw.pVector, stateDisc.pVector);
  BOOST_CHECK_EQUAL_COLLECTIONS(
      stateStraw.pVector.begin(), stateStraw.pVector.end(),
      statePerigee.pVector.begin(), statePerigee.pVector.end());
}

BOOST_AUTO_TEST_CASE(StepSize) {
  CurvilinearTrackParameters cp(pos4, unitDir, charge / absMom, cov,
                                particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options, cp);

  stepper.updateStepSize(state, -5_cm, ConstrainedStep::Type::Navigator);
  BOOST_CHECK_EQUAL(state.previousStepSize, stepSize);
  BOOST_CHECK_EQUAL(state.stepSize.value(), -5_cm);

  stepper.releaseStepSize(state, ConstrainedStep::Type::Navigator);
  BOOST_CHECK_EQUAL(state.stepSize.value(), stepSize);
}

// test step size modification with target surfaces
BOOST_AUTO_TEST_CASE(StepSizeSurface) {
  CurvilinearTrackParameters cp(pos4, unitDir, charge / absMom, cov,
                                particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options, cp);

  auto distance = 10_mm;
  auto target = CurvilinearSurface(pos + navDir * distance * unitDir, unitDir)
                    .planeSurface();

  stepper.updateSurfaceStatus(
      state, *target, 0, navDir, BoundaryTolerance::Infinite(),
      s_onSurfaceTolerance, ConstrainedStep::Type::Navigator);
  BOOST_CHECK_EQUAL(state.stepSize.value(ConstrainedStep::Type::Navigator),
                    distance);

  // test the step size modification in the context of a surface
  stepper.updateStepSize(state,
                         target
                             ->intersect(geoCtx, stepper.position(state),
                                         navDir * stepper.direction(state),
                                         BoundaryTolerance::Infinite())
                             .closest(),
                         navDir, ConstrainedStep::Type::Navigator);
  BOOST_CHECK_EQUAL(state.stepSize.value(), distance);

  // start with a different step size
  state.stepSize.setUser(navDir * stepSize);
  stepper.releaseStepSize(state, ConstrainedStep::Type::Navigator);
  stepper.updateStepSize(state,
                         target
                             ->intersect(geoCtx, stepper.position(state),
                                         navDir * stepper.direction(state),
                                         BoundaryTolerance::Infinite())
                             .closest(),
                         navDir, ConstrainedStep::Type::Navigator);
  BOOST_CHECK_EQUAL(state.stepSize.value(), navDir * stepSize);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
