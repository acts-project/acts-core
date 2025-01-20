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
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/GenericCurvilinearTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Result.hpp"

#include <limits>
#include <memory>
#include <optional>
#include <string>

using Acts::VectorHelpers::makeVector4;

namespace Acts::Test {

using Covariance = BoundSquareMatrix;

static constexpr auto eps = 2 * std::numeric_limits<double>::epsilon();

/// These tests are aiming to test whether the state setup is working properly
BOOST_AUTO_TEST_CASE(straight_line_stepper_state_test) {
  // Set up some variables
  GeometryContext tgContext = GeometryContext();
  MagneticFieldContext mfContext = MagneticFieldContext();
  double stepSize = 123.;

  Vector3 pos(1., 2., 3.);
  Vector3 dir(4., 5., 6.);
  double time = 7.;
  double absMom = 8.;
  double charge = -1.;

  StraightLineStepper::Options slsOptions(tgContext, mfContext);
  slsOptions.maxStepSize = stepSize;

  // Test charged parameters without covariance matrix
  CurvilinearTrackParameters cp(makeVector4(pos, time), dir, charge / absMom,
                                std::nullopt, ParticleHypothesis::pion());

  StraightLineStepper sls;
  StraightLineStepper::State slsState = sls.makeState(slsOptions, cp);

  // Test the result & compare with the input/test for reasonable members
  BOOST_CHECK_EQUAL(slsState.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_EQUAL(slsState.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(slsState.derivative, FreeVector::Zero());
  BOOST_CHECK(!slsState.covTransport);
  BOOST_CHECK_EQUAL(slsState.cov, Covariance::Zero());
  CHECK_CLOSE_OR_SMALL(sls.position(slsState), pos, eps, eps);
  CHECK_CLOSE_OR_SMALL(sls.direction(slsState), dir.normalized(), eps, eps);
  CHECK_CLOSE_REL(sls.absoluteMomentum(slsState), absMom, eps);
  BOOST_CHECK_EQUAL(sls.charge(slsState), charge);
  CHECK_CLOSE_OR_SMALL(sls.time(slsState), time, eps, eps);
  BOOST_CHECK_EQUAL(slsState.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(slsState.stepSize.value(), stepSize);
  BOOST_CHECK_EQUAL(slsState.previousStepSize, 0.);

  // Test without charge and covariance matrix
  CurvilinearTrackParameters ncp(makeVector4(pos, time), dir, 1 / absMom,
                                 std::nullopt, ParticleHypothesis::pion0());
  slsState = sls.makeState(slsOptions, ncp);

  // Test with covariance matrix
  Covariance cov = 8. * Covariance::Identity();
  ncp = CurvilinearTrackParameters(makeVector4(pos, time), dir, 1 / absMom, cov,
                                   ParticleHypothesis::pion0());
  slsState = sls.makeState(slsOptions, ncp);
  BOOST_CHECK_NE(slsState.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK(slsState.covTransport);
  BOOST_CHECK_EQUAL(slsState.cov, cov);
}

/// These tests are aiming to test the functions of the StraightLineStepper
/// The numerical correctness of the stepper is tested in the integration tests
BOOST_AUTO_TEST_CASE(straight_line_stepper_test) {
  // Set up some variables for the state
  GeometryContext tgContext = GeometryContext();
  MagneticFieldContext mfContext = MagneticFieldContext();
  Direction navDir = Direction::Backward();
  double stepSize = 123.;

  // Construct the parameters
  Vector3 pos(1., 2., 3.);
  Vector3 dir = Vector3(4., 5., 6.).normalized();
  double time = 7.;
  double absMom = 8.;
  double charge = -1.;
  Covariance cov = 8. * Covariance::Identity();
  CurvilinearTrackParameters cp(makeVector4(pos, time), dir, charge / absMom,
                                cov, ParticleHypothesis::pion());

  StraightLineStepper::Options options(tgContext, mfContext);
  options.maxStepSize = stepSize;

  // Build the stepper and the state
  StraightLineStepper sls;
  StraightLineStepper::State slsState = sls.makeState(options, cp);

  // Test the getters
  CHECK_CLOSE_ABS(sls.position(slsState), pos, 1e-6);
  CHECK_CLOSE_ABS(sls.direction(slsState), dir, 1e-6);
  CHECK_CLOSE_ABS(sls.absoluteMomentum(slsState), absMom, 1e-6);
  BOOST_CHECK_EQUAL(sls.charge(slsState), charge);
  BOOST_CHECK_EQUAL(sls.time(slsState), time);

  // Step size modifies
  const std::string originalStepSize = slsState.stepSize.toString();

  sls.updateStepSize(slsState, -1337., ConstrainedStep::Type::Navigator);
  BOOST_CHECK_EQUAL(slsState.previousStepSize, stepSize);
  BOOST_CHECK_EQUAL(slsState.stepSize.value(), -1337.);

  sls.releaseStepSize(slsState, ConstrainedStep::Type::Navigator);
  BOOST_CHECK_EQUAL(slsState.stepSize.value(), stepSize);
  BOOST_CHECK_EQUAL(sls.outputStepSize(slsState), originalStepSize);

  // Test the curvilinear state construction
  auto curvState = sls.curvilinearState(slsState);
  auto curvPars = std::get<0>(curvState);
  CHECK_CLOSE_ABS(curvPars.position(tgContext), cp.position(tgContext), 1e-6);
  CHECK_CLOSE_ABS(curvPars.absoluteMomentum(), cp.absoluteMomentum(), 1e-6);
  CHECK_CLOSE_ABS(curvPars.charge(), cp.charge(), 1e-6);
  CHECK_CLOSE_ABS(curvPars.time(), cp.time(), 1e-6);
  BOOST_CHECK(curvPars.covariance().has_value());
  BOOST_CHECK_NE(*curvPars.covariance(), cov);
  CHECK_CLOSE_COVARIANCE(std::get<1>(curvState),
                         BoundMatrix(BoundMatrix::Identity()), 1e-6);
  CHECK_CLOSE_ABS(std::get<2>(curvState), 0., 1e-6);

  // Test the update method
  Vector3 newPos(2., 4., 8.);
  Vector3 newMom(3., 9., 27.);
  double newTime(321.);
  sls.update(slsState, newPos, newMom.normalized(), charge / newMom.norm(),
             newTime);
  CHECK_CLOSE_ABS(sls.position(slsState), newPos, 1e-6);
  CHECK_CLOSE_ABS(sls.direction(slsState), newMom.normalized(), 1e-6);
  CHECK_CLOSE_ABS(sls.absoluteMomentum(slsState), newMom.norm(), 1e-6);
  BOOST_CHECK_EQUAL(sls.charge(slsState), charge);
  BOOST_CHECK_EQUAL(sls.time(slsState), newTime);

  // The covariance transport
  slsState.cov = cov;
  sls.transportCovarianceToCurvilinear(slsState);
  BOOST_CHECK_NE(slsState.cov, cov);
  BOOST_CHECK_NE(slsState.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_EQUAL(slsState.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(slsState.derivative, FreeVector::Zero());

  // Perform a step without and with covariance transport
  slsState.cov = cov;

  slsState.covTransport = false;
  double h = sls.step(slsState, navDir, nullptr).value();
  BOOST_CHECK_EQUAL(slsState.stepSize.value(), stepSize);
  BOOST_CHECK_EQUAL(slsState.stepSize.value(), h * navDir);
  CHECK_CLOSE_COVARIANCE(slsState.cov, cov, 1e-6);
  BOOST_CHECK_GT(sls.position(slsState).norm(), newPos.norm());
  CHECK_CLOSE_ABS(sls.direction(slsState), newMom.normalized(), 1e-6);
  CHECK_CLOSE_ABS(sls.absoluteMomentum(slsState), newMom.norm(), 1e-6);
  CHECK_CLOSE_ABS(sls.charge(slsState), charge, 1e-6);
  BOOST_CHECK_LT(sls.time(slsState), newTime);
  BOOST_CHECK_EQUAL(slsState.derivative, FreeVector::Zero());
  BOOST_CHECK_EQUAL(slsState.jacTransport, FreeMatrix::Identity());

  slsState.covTransport = true;
  double h2 = sls.step(slsState, navDir, nullptr).value();
  BOOST_CHECK_EQUAL(slsState.stepSize.value(), stepSize);
  BOOST_CHECK_EQUAL(h2, h);
  CHECK_CLOSE_COVARIANCE(slsState.cov, cov, 1e-6);
  BOOST_CHECK_GT(sls.position(slsState).norm(), newPos.norm());
  CHECK_CLOSE_ABS(sls.direction(slsState), newMom.normalized(), 1e-6);
  CHECK_CLOSE_ABS(sls.absoluteMomentum(slsState), newMom.norm(), 1e-6);
  CHECK_CLOSE_ABS(sls.charge(slsState), charge, 1e-6);
  BOOST_CHECK_LT(sls.time(slsState), newTime);
  BOOST_CHECK_NE(slsState.derivative, FreeVector::Zero());
  BOOST_CHECK_NE(slsState.jacTransport, FreeMatrix::Identity());

  /// Test the state reset
  // Construct the parameters
  Vector3 pos2(1.5, -2.5, 3.5);
  Vector3 dir2 = Vector3(4.5, -5.5, 6.5).normalized();
  double time2 = 7.5;
  double absMom2 = 8.5;
  double charge2 = 1.;
  BoundSquareMatrix cov2 = 8.5 * Covariance::Identity();
  CurvilinearTrackParameters cp2(makeVector4(pos2, time2), dir2,
                                 charge2 / absMom2, cov2,
                                 ParticleHypothesis::pion());
  BOOST_CHECK(cp2.covariance().has_value());
  FreeVector freeParams = transformBoundToFreeParameters(
      cp2.referenceSurface(), tgContext, cp2.parameters());
  navDir = Direction::Forward();
  double stepSize2 = -2. * stepSize;

  // Reset all possible parameters
  StraightLineStepper::State slsStateCopy = slsState;
  sls.resetState(slsStateCopy, cp2.parameters(), *cp2.covariance(),
                 cp2.referenceSurface(), stepSize2);
  // Test all components
  BOOST_CHECK_NE(slsStateCopy.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_NE(slsStateCopy.jacToGlobal, slsState.jacToGlobal);
  BOOST_CHECK_EQUAL(slsStateCopy.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(slsStateCopy.derivative, FreeVector::Zero());
  BOOST_CHECK(slsStateCopy.covTransport);
  BOOST_CHECK_EQUAL(slsStateCopy.cov, cov2);
  CHECK_CLOSE_ABS(sls.position(slsStateCopy),
                  freeParams.template segment<3>(eFreePos0), 1e-6);
  CHECK_CLOSE_ABS(sls.direction(slsStateCopy),
                  freeParams.template segment<3>(eFreeDir0).normalized(), 1e-6);
  CHECK_CLOSE_ABS(sls.absoluteMomentum(slsStateCopy),
                  std::abs(1. / freeParams[eFreeQOverP]), 1e-6);
  CHECK_CLOSE_ABS(sls.charge(slsStateCopy), -sls.charge(slsState), 1e-6);
  CHECK_CLOSE_ABS(sls.time(slsStateCopy), freeParams[eFreeTime], 1e-6);
  BOOST_CHECK_EQUAL(slsStateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(slsStateCopy.stepSize.value(), navDir * stepSize2);
  BOOST_CHECK_EQUAL(slsStateCopy.previousStepSize, slsState.previousStepSize);

  // Reset all possible parameters except the step size
  slsStateCopy = slsState;
  sls.resetState(slsStateCopy, cp2.parameters(), *cp2.covariance(),
                 cp2.referenceSurface());
  // Test all components
  BOOST_CHECK_NE(slsStateCopy.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_NE(slsStateCopy.jacToGlobal, slsState.jacToGlobal);
  BOOST_CHECK_EQUAL(slsStateCopy.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(slsStateCopy.derivative, FreeVector::Zero());
  BOOST_CHECK(slsStateCopy.covTransport);
  BOOST_CHECK_EQUAL(slsStateCopy.cov, cov2);
  CHECK_CLOSE_ABS(sls.position(slsStateCopy),
                  freeParams.template segment<3>(eFreePos0), 1e-6);
  CHECK_CLOSE_ABS(sls.direction(slsStateCopy),
                  freeParams.template segment<3>(eFreeDir0), 1e-6);
  CHECK_CLOSE_ABS(sls.absoluteMomentum(slsStateCopy),
                  std::abs(1. / freeParams[eFreeQOverP]), 1e-6);
  CHECK_CLOSE_ABS(sls.charge(slsStateCopy), -sls.charge(slsState), 1e-6);
  CHECK_CLOSE_ABS(sls.time(slsStateCopy), freeParams[eFreeTime], 1e-6);
  BOOST_CHECK_EQUAL(slsStateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(slsStateCopy.stepSize.value(),
                    std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(slsStateCopy.previousStepSize, slsState.previousStepSize);

  // Reset the least amount of parameters
  slsStateCopy = slsState;
  sls.resetState(slsStateCopy, cp2.parameters(), *cp2.covariance(),
                 cp2.referenceSurface());
  // Test all components
  BOOST_CHECK_NE(slsStateCopy.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_NE(slsStateCopy.jacToGlobal, slsState.jacToGlobal);
  BOOST_CHECK_EQUAL(slsStateCopy.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(slsStateCopy.derivative, FreeVector::Zero());
  BOOST_CHECK(slsStateCopy.covTransport);
  BOOST_CHECK_EQUAL(slsStateCopy.cov, cov2);
  CHECK_CLOSE_ABS(sls.position(slsStateCopy),
                  freeParams.template segment<3>(eFreePos0), 1e-6);
  CHECK_CLOSE_ABS(sls.direction(slsStateCopy),
                  freeParams.template segment<3>(eFreeDir0).normalized(), 1e-6);
  CHECK_CLOSE_ABS(sls.absoluteMomentum(slsStateCopy),
                  std::abs(1. / freeParams[eFreeQOverP]), 1e-6);
  CHECK_CLOSE_ABS(sls.charge(slsStateCopy), -sls.charge(slsState), 1e-6);
  CHECK_CLOSE_ABS(sls.time(slsStateCopy), freeParams[eFreeTime], 1e-6);
  BOOST_CHECK_EQUAL(slsStateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(slsStateCopy.stepSize.value(),
                    std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(slsStateCopy.previousStepSize, slsState.previousStepSize);

  /// Repeat with surface related methods
  auto plane = CurvilinearSurface(pos, dir).planeSurface();
  auto bp = BoundTrackParameters::create(
                plane, tgContext, makeVector4(pos, time), dir, charge / absMom,
                cov, ParticleHypothesis::pion())
                .value();
  slsState = sls.makeState(options, bp);

  // Test the intersection in the context of a surface
  auto targetSurface =
      CurvilinearSurface(pos + navDir * 2. * dir, dir).planeSurface();
  sls.updateSurfaceStatus(slsState, *targetSurface, 0, navDir,
                          BoundaryTolerance::Infinite(), s_onSurfaceTolerance,
                          ConstrainedStep::Type::Navigator);
  CHECK_CLOSE_ABS(slsState.stepSize.value(ConstrainedStep::Type::Navigator),
                  navDir * 2., 1e-6);

  // Test the step size modification in the context of a surface
  sls.updateStepSize(slsState,
                     targetSurface
                         ->intersect(tgContext, sls.position(slsState),
                                     navDir * sls.direction(slsState),
                                     BoundaryTolerance::Infinite())
                         .closest(),
                     navDir, ConstrainedStep::Type::Navigator);
  CHECK_CLOSE_ABS(slsState.stepSize.value(), 2, 1e-6);
  slsState.stepSize.setUser(navDir * stepSize);
  sls.releaseStepSize(slsState, ConstrainedStep::Type::Navigator);
  sls.updateStepSize(slsState,
                     targetSurface
                         ->intersect(tgContext, sls.position(slsState),
                                     navDir * sls.direction(slsState),
                                     BoundaryTolerance::Infinite())
                         .closest(),
                     navDir, ConstrainedStep::Type::Navigator);
  CHECK_CLOSE_ABS(slsState.stepSize.value(), 2, 1e-6);

  // Test the bound state construction
  FreeToBoundCorrection freeToBoundCorrection(false);
  auto boundState =
      sls.boundState(slsState, *plane, true, freeToBoundCorrection).value();
  auto boundPars = std::get<0>(boundState);
  CHECK_CLOSE_ABS(boundPars.position(tgContext), bp.position(tgContext), 1e-6);
  CHECK_CLOSE_ABS(boundPars.momentum(), bp.momentum(), 1e-6);
  CHECK_CLOSE_ABS(boundPars.charge(), bp.charge(), 1e-6);
  CHECK_CLOSE_ABS(boundPars.time(), bp.time(), 1e-6);
  BOOST_CHECK(boundPars.covariance().has_value());
  BOOST_CHECK_NE(*boundPars.covariance(), cov);
  CHECK_CLOSE_COVARIANCE(std::get<1>(boundState),
                         BoundMatrix(BoundMatrix::Identity()), 1e-6);
  CHECK_CLOSE_ABS(std::get<2>(boundState), 0., 1e-6);

  // Transport the covariance in the context of a surface
  sls.transportCovarianceToBound(slsState, *plane, freeToBoundCorrection);
  BOOST_CHECK_NE(slsState.cov, cov);
  BOOST_CHECK_NE(slsState.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_EQUAL(slsState.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(slsState.derivative, FreeVector::Zero());

  // Update in context of a surface
  freeParams = transformBoundToFreeParameters(bp.referenceSurface(), tgContext,
                                              bp.parameters());

  BOOST_CHECK(bp.covariance().has_value());
  sls.update(slsState, freeParams, bp.parameters(), 2 * (*bp.covariance()),
             *plane);
  CHECK_CLOSE_OR_SMALL(sls.position(slsState), pos, eps, eps);
  BOOST_CHECK_EQUAL(sls.charge(slsState), charge);
  CHECK_CLOSE_OR_SMALL(sls.time(slsState), time, eps, eps);
  CHECK_CLOSE_COVARIANCE(slsState.cov, Covariance(2. * cov), 1e-6);
}

}  // namespace Acts::Test
