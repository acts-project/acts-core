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
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/GenericCurvilinearTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/EigenStepperDenseExtension.hpp"
#include "Acts/Propagator/EigenStepperError.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <cmath>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace Acts::UnitLiterals;
using Acts::VectorHelpers::makeVector4;

namespace Acts::Test {

using Covariance = BoundSquareMatrix;

static constexpr auto eps = 3 * std::numeric_limits<double>::epsilon();

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

/// @brief Aborter for the case that a particle leaves the detector or reaches
/// a custom made threshold.
///
struct EndOfWorld {
  /// Maximum value in x-direction of the detector
  double maxX = 1_m;

  /// @brief Main call operator for the abort operation
  ///
  /// @tparam propagator_state_t State of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in] state State of the propagation
  /// @param [in] stepper Stepper of the propagation
  /// @param [in] navigator Navigator of the propagation
  ///
  /// @return Boolean statement if the particle is still in the detector
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool checkAbort(propagator_state_t& state, const stepper_t& stepper,
                  const navigator_t& /*navigator*/,
                  const Logger& /*logger*/) const {
    const double tolerance = state.options.surfaceTolerance;
    if (maxX - std::abs(stepper.position(state.stepping).x()) <= tolerance ||
        std::abs(stepper.position(state.stepping).y()) >= 0.5_m ||
        std::abs(stepper.position(state.stepping).z()) >= 0.5_m) {
      return true;
    }
    return false;
  }
};

///
/// @brief Data collector while propagation
///
struct StepCollector {
  ///
  /// @brief Data container for result analysis
  ///
  struct this_result {
    // Position of the propagator after each step
    std::vector<Vector3> position;
    // Momentum of the propagator after each step
    std::vector<Vector3> momentum;
  };

  using result_type = this_result;

  /// @brief Main call operator for the action list. It stores the data for
  /// analysis afterwards
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in] state State of the propagator
  /// @param [in] stepper Stepper of the propagation
  /// @param [in] navigator Navigator of the propagation
  /// @param [out] result Struct which is filled with the data
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  void act(propagator_state_t& state, const stepper_t& stepper,
           const navigator_t& /*navigator*/, result_type& result,
           const Logger& /*logger*/) const {
    result.position.push_back(stepper.position(state.stepping));
    result.momentum.push_back(stepper.momentum(state.stepping));
  }

  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool checkAbort(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                  const navigator_t& /*navigator*/, result_type& /*result*/,
                  const Logger& /*logger*/) const {
    return false;
  }
};

/// These tests are aiming to test whether the state setup is working properly
BOOST_AUTO_TEST_CASE(eigen_stepper_state_test) {
  // Set up some variables
  auto bField = std::make_shared<ConstantBField>(Vector3(1., 2.5, 33.33));

  Vector3 pos(1., 2., 3.);
  Vector3 dir(4., 5., 6.);
  double time = 7.;
  double absMom = 8.;
  double charge = -1.;

  EigenStepper<>::Options esOptions(tgContext, mfContext);

  // Test charged parameters without covariance matrix
  CurvilinearTrackParameters cp(makeVector4(pos, time), dir, charge / absMom,
                                std::nullopt, ParticleHypothesis::pion());
  EigenStepper<> es(bField);
  EigenStepper<>::State esState = es.makeState(esOptions, cp);

  // Test the result & compare with the input/test for reasonable members
  BOOST_CHECK_EQUAL(esState.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_EQUAL(esState.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(esState.derivative, FreeVector::Zero());
  BOOST_CHECK(!esState.covTransport);
  BOOST_CHECK_EQUAL(esState.cov, Covariance::Zero());
  BOOST_CHECK_EQUAL(esState.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(esState.previousStepSize, 0.);

  // Test without charge and covariance matrix
  CurvilinearTrackParameters ncp(makeVector4(pos, time), dir, 1 / absMom,
                                 std::nullopt, ParticleHypothesis::pion0());
  esOptions = EigenStepper<>::Options(tgContext, mfContext);
  esState = es.makeState(esOptions, ncp);
  BOOST_CHECK_EQUAL(es.charge(esState), 0.);

  // Test with covariance matrix
  Covariance cov = 8. * Covariance::Identity();
  ncp = CurvilinearTrackParameters(makeVector4(pos, time), dir, 1 / absMom, cov,
                                   ParticleHypothesis::pion0());
  esState = es.makeState(esOptions, ncp);
  BOOST_CHECK_NE(esState.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK(esState.covTransport);
  BOOST_CHECK_EQUAL(esState.cov, cov);
}

/// These tests are aiming to test the functions of the EigenStepper
/// The numerical correctness of the stepper is tested in the integration tests
BOOST_AUTO_TEST_CASE(eigen_stepper_test) {
  // Set up some variables for the state
  Direction navDir = Direction::Backward();
  double stepSize = 123.;
  auto bField = std::make_shared<ConstantBField>(Vector3(1., 2.5, 33.33));
  auto bCache = bField->makeCache(mfContext);

  // Construct the parameters
  Vector3 pos(1., 2., 3.);
  Vector3 dir = Vector3(4., 5., 6.).normalized();
  double time = 7.;
  double absMom = 8.;
  double charge = -1.;
  Covariance cov = 8. * Covariance::Identity();
  CurvilinearTrackParameters cp(makeVector4(pos, time), dir, charge / absMom,
                                cov, ParticleHypothesis::pion());

  EigenStepper<>::Options esOptions(tgContext, mfContext);
  esOptions.maxStepSize = stepSize;

  // Build the stepper and the state
  EigenStepper<> es(bField);
  EigenStepper<>::State esState = es.makeState(esOptions, cp);

  // Test the getters
  CHECK_CLOSE_ABS(es.position(esState), pos, eps);
  CHECK_CLOSE_ABS(es.direction(esState), dir, eps);
  CHECK_CLOSE_ABS(es.absoluteMomentum(esState), absMom, eps);
  CHECK_CLOSE_ABS(es.charge(esState), charge, eps);
  CHECK_CLOSE_ABS(es.time(esState), time, eps);
  BOOST_CHECK_EQUAL(es.getField(esState, pos).value(),
                    bField->getField(pos, bCache).value());

  // Step size modifies
  const std::string originalStepSize = esState.stepSize.toString();

  es.updateStepSize(esState, -1337., ConstrainedStep::Type::Navigator);
  BOOST_CHECK_EQUAL(esState.previousStepSize, stepSize);
  BOOST_CHECK_EQUAL(esState.stepSize.value(), -1337.);

  es.releaseStepSize(esState, ConstrainedStep::Type::Navigator);
  BOOST_CHECK_EQUAL(esState.stepSize.value(), stepSize);
  BOOST_CHECK_EQUAL(es.outputStepSize(esState), originalStepSize);

  // Test the curvilinear state construction
  auto curvState = es.curvilinearState(esState);
  auto curvPars = std::get<0>(curvState);
  CHECK_CLOSE_ABS(curvPars.position(tgContext), cp.position(tgContext), eps);
  CHECK_CLOSE_ABS(curvPars.momentum(), cp.momentum(), 10e-6);
  CHECK_CLOSE_ABS(curvPars.charge(), cp.charge(), eps);
  CHECK_CLOSE_ABS(curvPars.time(), cp.time(), eps);
  BOOST_CHECK(curvPars.covariance().has_value());
  BOOST_CHECK_NE(*curvPars.covariance(), cov);
  CHECK_CLOSE_COVARIANCE(std::get<1>(curvState),
                         BoundMatrix(BoundMatrix::Identity()), eps);
  CHECK_CLOSE_ABS(std::get<2>(curvState), 0., eps);

  // Test the update method
  Vector3 newPos(2., 4., 8.);
  Vector3 newMom(3., 9., 27.);
  double newTime(321.);
  es.update(esState, newPos, newMom.normalized(), charge / newMom.norm(),
            newTime);
  BOOST_CHECK_EQUAL(es.position(esState), newPos);
  BOOST_CHECK_EQUAL(es.direction(esState), newMom.normalized());
  BOOST_CHECK_EQUAL(es.absoluteMomentum(esState), newMom.norm());
  BOOST_CHECK_EQUAL(es.charge(esState), charge);
  BOOST_CHECK_EQUAL(es.time(esState), newTime);

  // The covariance transport
  esState.cov = cov;
  es.transportCovarianceToCurvilinear(esState);
  BOOST_CHECK_NE(esState.cov, cov);
  BOOST_CHECK_NE(esState.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_EQUAL(esState.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(esState.derivative, FreeVector::Zero());

  // Perform a step without and with covariance transport
  esState.cov = cov;

  esState.covTransport = false;
  es.step(esState, navDir, nullptr).value();
  CHECK_CLOSE_COVARIANCE(esState.cov, cov, eps);
  BOOST_CHECK_NE(es.position(esState).norm(), newPos.norm());
  BOOST_CHECK_NE(es.direction(esState), newMom.normalized());
  BOOST_CHECK_EQUAL(es.charge(esState), charge);
  BOOST_CHECK_LT(es.time(esState), newTime);
  BOOST_CHECK_EQUAL(esState.derivative, FreeVector::Zero());
  BOOST_CHECK_EQUAL(esState.jacTransport, FreeMatrix::Identity());

  esState.covTransport = true;
  es.step(esState, navDir, nullptr).value();
  CHECK_CLOSE_COVARIANCE(esState.cov, cov, eps);
  BOOST_CHECK_NE(es.position(esState).norm(), newPos.norm());
  BOOST_CHECK_NE(es.direction(esState), newMom.normalized());
  BOOST_CHECK_EQUAL(es.charge(esState), charge);
  BOOST_CHECK_LT(es.time(esState), newTime);
  BOOST_CHECK_NE(esState.derivative, FreeVector::Zero());
  BOOST_CHECK_NE(esState.jacTransport, FreeMatrix::Identity());

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
  FreeVector freeParams = transformBoundToFreeParameters(
      cp2.referenceSurface(), tgContext, cp2.parameters());
  navDir = Direction::Forward();
  double stepSize2 = -2. * stepSize;

  auto copyState = [&](auto& field, const auto& state) {
    using field_t = std::decay_t<decltype(field)>;
    std::decay_t<decltype(state)> copy = es.makeState(esOptions, cp);
    copy.pars = state.pars;
    copy.covTransport = state.covTransport;
    copy.cov = state.cov;
    copy.jacobian = state.jacobian;
    copy.jacToGlobal = state.jacToGlobal;
    copy.jacTransport = state.jacTransport;
    copy.derivative = state.derivative;
    copy.pathAccumulated = state.pathAccumulated;
    copy.stepSize = state.stepSize;
    copy.previousStepSize = state.previousStepSize;

    copy.fieldCache = MagneticFieldProvider::Cache(
        std::in_place_type<typename field_t::Cache>,
        state.fieldCache.template as<typename field_t::Cache>());

    copy.extension = state.extension;
    copy.stepData = state.stepData;

    return copy;
  };

  // Reset all possible parameters
  EigenStepper<>::State esStateCopy = copyState(*bField, esState);
  BOOST_CHECK(cp2.covariance().has_value());
  es.resetState(esStateCopy, cp2.parameters(), *cp2.covariance(),
                cp2.referenceSurface(), stepSize2);
  // Test all components
  BOOST_CHECK_NE(esStateCopy.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_NE(esStateCopy.jacToGlobal, esState.jacToGlobal);
  BOOST_CHECK_EQUAL(esStateCopy.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(esStateCopy.derivative, FreeVector::Zero());
  BOOST_CHECK(esStateCopy.covTransport);
  BOOST_CHECK_EQUAL(esStateCopy.cov, cov2);
  BOOST_CHECK_EQUAL(es.position(esStateCopy),
                    freeParams.template segment<3>(eFreePos0));
  BOOST_CHECK_EQUAL(es.direction(esStateCopy),
                    freeParams.template segment<3>(eFreeDir0).normalized());
  BOOST_CHECK_EQUAL(es.absoluteMomentum(esStateCopy),
                    std::abs(1. / freeParams[eFreeQOverP]));
  BOOST_CHECK_EQUAL(es.charge(esStateCopy), -es.charge(esState));
  BOOST_CHECK_EQUAL(es.time(esStateCopy), freeParams[eFreeTime]);
  BOOST_CHECK_EQUAL(esStateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(esStateCopy.stepSize.value(), navDir * stepSize2);
  BOOST_CHECK_EQUAL(esStateCopy.previousStepSize, esState.previousStepSize);

  // Reset all possible parameters except the step size
  esStateCopy = copyState(*bField, esState);
  es.resetState(esStateCopy, cp2.parameters(), *cp2.covariance(),
                cp2.referenceSurface());
  // Test all components
  BOOST_CHECK_NE(esStateCopy.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_NE(esStateCopy.jacToGlobal, esState.jacToGlobal);
  BOOST_CHECK_EQUAL(esStateCopy.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(esStateCopy.derivative, FreeVector::Zero());
  BOOST_CHECK(esStateCopy.covTransport);
  BOOST_CHECK_EQUAL(esStateCopy.cov, cov2);
  BOOST_CHECK_EQUAL(es.position(esStateCopy),
                    freeParams.template segment<3>(eFreePos0));
  BOOST_CHECK_EQUAL(es.direction(esStateCopy),
                    freeParams.template segment<3>(eFreeDir0).normalized());
  BOOST_CHECK_EQUAL(es.absoluteMomentum(esStateCopy),
                    std::abs(1. / freeParams[eFreeQOverP]));
  BOOST_CHECK_EQUAL(es.charge(esStateCopy), -es.charge(esState));
  BOOST_CHECK_EQUAL(es.time(esStateCopy), freeParams[eFreeTime]);
  BOOST_CHECK_EQUAL(esStateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(esStateCopy.stepSize.value(),
                    std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(esStateCopy.previousStepSize, esState.previousStepSize);

  // Reset the least amount of parameters
  esStateCopy = copyState(*bField, esState);
  es.resetState(esStateCopy, cp2.parameters(), *cp2.covariance(),
                cp2.referenceSurface());
  // Test all components
  BOOST_CHECK_NE(esStateCopy.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_NE(esStateCopy.jacToGlobal, esState.jacToGlobal);
  BOOST_CHECK_EQUAL(esStateCopy.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(esStateCopy.derivative, FreeVector::Zero());
  BOOST_CHECK(esStateCopy.covTransport);
  BOOST_CHECK_EQUAL(esStateCopy.cov, cov2);
  BOOST_CHECK_EQUAL(es.position(esStateCopy),
                    freeParams.template segment<3>(eFreePos0));
  BOOST_CHECK_EQUAL(es.direction(esStateCopy),
                    freeParams.template segment<3>(eFreeDir0).normalized());
  BOOST_CHECK_EQUAL(es.absoluteMomentum(esStateCopy),
                    std::abs(1. / freeParams[eFreeQOverP]));
  BOOST_CHECK_EQUAL(es.charge(esStateCopy), -es.charge(esState));
  BOOST_CHECK_EQUAL(es.time(esStateCopy), freeParams[eFreeTime]);
  BOOST_CHECK_EQUAL(esStateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(esStateCopy.stepSize.value(),
                    std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(esStateCopy.previousStepSize, esState.previousStepSize);

  /// Repeat with surface related methods
  auto plane = CurvilinearSurface(pos, dir.normalized()).planeSurface();
  auto bp = BoundTrackParameters::create(
                plane, tgContext, makeVector4(pos, time), dir, charge / absMom,
                cov, ParticleHypothesis::pion())
                .value();
  esState = es.makeState(esOptions, bp);

  // Test the intersection in the context of a surface
  auto targetSurface =
      CurvilinearSurface(pos + navDir * 2. * dir, dir).planeSurface();
  es.updateSurfaceStatus(esState, *targetSurface, 0, navDir,
                         BoundaryTolerance::Infinite(), s_onSurfaceTolerance,
                         ConstrainedStep::Type::Navigator);
  CHECK_CLOSE_ABS(esState.stepSize.value(ConstrainedStep::Type::Navigator),
                  navDir * 2., eps);

  // Test the step size modification in the context of a surface
  es.updateStepSize(esState,
                    targetSurface
                        ->intersect(tgContext, es.position(esState),
                                    navDir * es.direction(esState),
                                    BoundaryTolerance::Infinite())
                        .closest(),
                    navDir, ConstrainedStep::Type::Navigator);
  CHECK_CLOSE_ABS(esState.stepSize.value(), 2., eps);
  esState.stepSize.setUser(navDir * stepSize);
  es.releaseStepSize(esState, ConstrainedStep::Type::Navigator);
  es.updateStepSize(esState,
                    targetSurface
                        ->intersect(tgContext, es.position(esState),
                                    navDir * es.direction(esState),
                                    BoundaryTolerance::Infinite())
                        .closest(),
                    navDir, ConstrainedStep::Type::Navigator);
  CHECK_CLOSE_ABS(esState.stepSize.value(), 2., eps);

  // Test the bound state construction
  auto boundState = es.boundState(esState, *plane).value();
  auto boundPars = std::get<0>(boundState);
  CHECK_CLOSE_ABS(boundPars.position(tgContext), bp.position(tgContext), eps);
  CHECK_CLOSE_ABS(boundPars.momentum(), bp.momentum(), 1e-7);
  CHECK_CLOSE_ABS(boundPars.charge(), bp.charge(), eps);
  CHECK_CLOSE_ABS(boundPars.time(), bp.time(), eps);
  BOOST_CHECK(boundPars.covariance().has_value());
  BOOST_CHECK_NE(*boundPars.covariance(), cov);
  CHECK_CLOSE_COVARIANCE(std::get<1>(boundState),
                         BoundMatrix(BoundMatrix::Identity()), eps);
  CHECK_CLOSE_ABS(std::get<2>(boundState), 0., eps);

  // Transport the covariance in the context of a surface
  es.transportCovarianceToBound(esState, *plane);
  BOOST_CHECK_NE(esState.cov, cov);
  BOOST_CHECK_NE(esState.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_EQUAL(esState.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(esState.derivative, FreeVector::Zero());

  // Update in context of a surface
  freeParams = transformBoundToFreeParameters(bp.referenceSurface(), tgContext,
                                              bp.parameters());

  es.update(esState, freeParams, bp.parameters(), 2 * (*bp.covariance()),
            *plane);
  CHECK_CLOSE_OR_SMALL(es.position(esState), pos, eps, eps);
  CHECK_CLOSE_OR_SMALL(es.direction(esState), dir, eps, eps);
  CHECK_CLOSE_REL(es.absoluteMomentum(esState), absMom, eps);
  BOOST_CHECK_EQUAL(es.charge(esState), charge);
  CHECK_CLOSE_OR_SMALL(es.time(esState), time, eps, eps);
  CHECK_CLOSE_COVARIANCE(esState.cov, Covariance(2. * cov), eps);

  // Test a case where no step size adjustment is required
  esState.options.stepTolerance = 2. * 4.4258e+09;
  double h0 = esState.stepSize.value();
  es.step(esState, navDir, nullptr);
  CHECK_CLOSE_ABS(h0, esState.stepSize.value(), eps);

  // Produce some errors
  auto nBfield = std::make_shared<NullBField>();
  EigenStepper<> nes(nBfield);
  EigenStepper<>::Options nesOptions(tgContext, mfContext);
  EigenStepper<>::State nesState = nes.makeState(nesOptions, cp);
  auto nEsState = copyState(*nBfield, nesState);
  // Test that we can reach the minimum step size
  nEsState.options.stepTolerance = 1e-21;
  nEsState.options.stepSizeCutOff = 1e20;
  auto res = nes.step(nEsState, navDir, nullptr);
  BOOST_CHECK(!res.ok());
  BOOST_CHECK_EQUAL(res.error(), EigenStepperError::StepSizeStalled);

  // Test that the number of trials exceeds
  nEsState.options.stepSizeCutOff = 0.;
  nEsState.options.maxRungeKuttaStepTrials = 0.;
  res = nes.step(nEsState, navDir, nullptr);
  BOOST_CHECK(!res.ok());
  BOOST_CHECK_EQUAL(res.error(), EigenStepperError::StepSizeAdjustmentFailed);
}

/// @brief This function tests the EigenStepper with the EigenStepperDefaultExtension and
/// the EigenStepperDenseExtension. The focus of this tests lies in
/// the choosing of the right extension for the individual use case. This is
/// performed with three different detectors:
/// a) Pure vacuum -> DefaultExtension needs to act
/// b) Pure Be -> DenseEnvironmentExtension needs to act
/// c) Vacuum - Be - Vacuum -> Both should act and switch during the propagation

// Test case a). The DenseEnvironmentExtension should state that it is not
// valid in this case.
BOOST_AUTO_TEST_CASE(step_extension_vacuum_test) {
  CuboidVolumeBuilder cvb;
  CuboidVolumeBuilder::VolumeConfig vConf;
  vConf.position = {0.5_m, 0., 0.};
  vConf.length = {1_m, 1_m, 1_m};
  CuboidVolumeBuilder::Config conf;
  conf.volumeCfg.push_back(vConf);
  conf.position = {0.5_m, 0., 0.};
  conf.length = {1_m, 1_m, 1_m};

  // Build detector
  cvb.setConfig(conf);
  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto& vb) {
        return cvb.trackingVolume(context, inner, vb);
      });
  TrackingGeometryBuilder tgb(tgbCfg);
  std::shared_ptr<const TrackingGeometry> vacuum =
      tgb.trackingGeometry(tgContext);

  // Build navigator
  Navigator naviVac({vacuum, true, true, true});

  // Set initial parameters for the particle track
  Covariance cov = Covariance::Identity();
  const Vector3 startDir = makeDirectionFromPhiTheta(0_degree, 90_degree);
  const Vector3 startMom = 1_GeV * startDir;
  const CurvilinearTrackParameters sbtp(Vector4::Zero(), startDir, 1_e / 1_GeV,
                                        cov, ParticleHypothesis::pion());

  using Stepper = EigenStepper<EigenStepperDenseExtension>;
  using Propagator = Propagator<Stepper, Navigator>;
  using PropagatorOptions =
      Propagator::Options<ActorList<StepCollector, EndOfWorld>>;

  // Set options for propagator
  PropagatorOptions propOpts(tgContext, mfContext);
  propOpts.maxSteps = 100;
  propOpts.stepping.maxStepSize = 1.5_m;

  // Build stepper and propagator
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 0.));
  Stepper es(bField);
  Propagator prop(es, naviVac);

  // Launch and collect results
  const auto& result = prop.propagate(sbtp, propOpts).value();
  const StepCollector::this_result& stepResult =
      result.get<typename StepCollector::result_type>();

  // Check that the propagation happened without interactions
  for (const auto& pos : stepResult.position) {
    CHECK_SMALL(pos.y(), 1_um);
    CHECK_SMALL(pos.z(), 1_um);
    if (pos == stepResult.position.back()) {
      CHECK_CLOSE_ABS(pos.x(), 1_m, 1_um);
    }
  }
  for (const auto& mom : stepResult.momentum) {
    CHECK_CLOSE_ABS(mom, startMom, 1_keV);
  }

  using DefStepper = EigenStepper<EigenStepperDenseExtension>;
  using DefPropagator = Acts::Propagator<DefStepper, Navigator>;
  using DefPropagatorOptions =
      DefPropagator::Options<ActorList<StepCollector, EndOfWorld>>;

  // Set options for propagator
  DefPropagatorOptions propOptsDef(tgContext, mfContext);
  propOptsDef.maxSteps = 100;
  propOptsDef.stepping.maxStepSize = 1.5_m;

  DefStepper esDef(bField);
  DefPropagator propDef(esDef, naviVac);

  // Launch and collect results
  const auto& resultDef = propDef.propagate(sbtp, propOptsDef).value();
  const StepCollector::this_result& stepResultDef =
      resultDef.get<typename StepCollector::result_type>();

  // Check that the right extension was chosen
  // If chosen correctly, the number of elements should be identical
  BOOST_CHECK_EQUAL(stepResult.position.size(), stepResultDef.position.size());
  for (unsigned int i = 0; i < stepResult.position.size(); i++) {
    CHECK_CLOSE_ABS(stepResult.position[i], stepResultDef.position[i], 1_um);
  }
  BOOST_CHECK_EQUAL(stepResult.momentum.size(), stepResultDef.momentum.size());
  for (unsigned int i = 0; i < stepResult.momentum.size(); i++) {
    CHECK_CLOSE_ABS(stepResult.momentum[i], stepResultDef.momentum[i], 1_keV);
  }
}
// Test case b). The DefaultExtension should state that it is invalid here.
BOOST_AUTO_TEST_CASE(step_extension_material_test) {
  CuboidVolumeBuilder cvb;
  CuboidVolumeBuilder::VolumeConfig vConf;
  vConf.position = {0.5_m, 0., 0.};
  vConf.length = {1_m, 1_m, 1_m};
  vConf.volumeMaterial =
      std::make_shared<HomogeneousVolumeMaterial>(makeBeryllium());
  CuboidVolumeBuilder::Config conf;
  conf.volumeCfg.push_back(vConf);
  conf.position = {0.5_m, 0., 0.};
  conf.length = {1_m, 1_m, 1_m};

  // Build detector
  cvb.setConfig(conf);
  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto& vb) {
        return cvb.trackingVolume(context, inner, vb);
      });
  TrackingGeometryBuilder tgb(tgbCfg);
  std::shared_ptr<const TrackingGeometry> material =
      tgb.trackingGeometry(tgContext);

  // Build navigator
  Navigator naviMat({material, true, true, true});

  // Set initial parameters for the particle track
  Covariance cov = Covariance::Identity();
  const Vector3 startDir = makeDirectionFromPhiTheta(0_degree, 90_degree);
  const Vector3 startMom = 5_GeV * startDir;
  const CurvilinearTrackParameters sbtp(Vector4::Zero(), startDir, 1_e / 5_GeV,
                                        cov, ParticleHypothesis::pion());

  using Stepper = EigenStepper<EigenStepperDenseExtension>;
  using Propagator = Propagator<Stepper, Navigator>;
  using PropagatorOptions =
      Propagator::Options<ActorList<StepCollector, EndOfWorld>>;

  // Set options for propagator
  PropagatorOptions propOpts(tgContext, mfContext);
  propOpts.maxSteps = 10000;
  propOpts.stepping.maxStepSize = 1.5_m;

  // Build stepper and propagator
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 0.));
  Stepper es(bField);
  Propagator prop(es, naviMat,
                  Acts::getDefaultLogger("Propagator", Acts::Logging::VERBOSE));

  // Launch and collect results
  const auto& result = prop.propagate(sbtp, propOpts).value();
  const StepCollector::this_result& stepResult =
      result.get<typename StepCollector::result_type>();

  // Check that there occurred interaction
  for (const auto& pos : stepResult.position) {
    CHECK_SMALL(pos.y(), 1_um);
    CHECK_SMALL(pos.z(), 1_um);
    if (pos == stepResult.position.front()) {
      CHECK_SMALL(pos.x(), 1_um);
    } else {
      BOOST_CHECK_GT(std::abs(pos.x()), 1_um);
    }
  }
  for (const auto& mom : stepResult.momentum) {
    CHECK_SMALL(mom.y(), 1_keV);
    CHECK_SMALL(mom.z(), 1_keV);
    if (mom == stepResult.momentum.front()) {
      CHECK_CLOSE_ABS(mom.x(), 5_GeV, 1_keV);
    } else {
      BOOST_CHECK_LT(mom.x(), 5_GeV);
    }
  }

  using DenseStepper = EigenStepper<EigenStepperDenseExtension>;
  using DensePropagator = Acts::Propagator<DenseStepper, Navigator>;
  using DensePropagatorOptions =
      DensePropagator::Options<ActorList<StepCollector, EndOfWorld>>;

  // Rebuild and check the choice of extension
  // Set options for propagator
  DensePropagatorOptions propOptsDense(tgContext, mfContext);
  propOptsDense.maxSteps = 1000;
  propOptsDense.stepping.maxStepSize = 1.5_m;

  // Build stepper and propagator
  DenseStepper esDense(bField);
  DensePropagator propDense(esDense, naviMat);

  // Launch and collect results
  const auto& resultDense = propDense.propagate(sbtp, propOptsDense).value();
  const StepCollector::this_result& stepResultDense =
      resultDense.get<typename StepCollector::result_type>();

  // Check that the right extension was chosen
  // If chosen correctly, the number of elements should be identical
  BOOST_CHECK_EQUAL(stepResult.position.size(),
                    stepResultDense.position.size());
  for (unsigned int i = 0; i < stepResult.position.size(); i++) {
    CHECK_CLOSE_ABS(stepResult.position[i], stepResultDense.position[i], 1_um);
  }
  BOOST_CHECK_EQUAL(stepResult.momentum.size(),
                    stepResultDense.momentum.size());
  for (unsigned int i = 0; i < stepResult.momentum.size(); i++) {
    CHECK_CLOSE_ABS(stepResult.momentum[i], stepResultDense.momentum[i], 1_keV);
  }

  ////////////////////////////////////////////////////////////////////

  // Re-launch the configuration with magnetic field
  bField->setField(Vector3{0., 1_T, 0.});
  Stepper esB(bField);
  Propagator propB(esB, naviMat);

  const auto& resultB = propB.propagate(sbtp, propOptsDense).value();
  const StepCollector::this_result& stepResultB =
      resultB.get<typename StepCollector::result_type>();

  // Check that there occurred interaction
  for (const auto& pos : stepResultB.position) {
    if (pos == stepResultB.position.front()) {
      CHECK_SMALL(pos, 1_um);
    } else {
      BOOST_CHECK_GT(std::abs(pos.x()), 1_um);
      CHECK_SMALL(pos.y(), 1_um);
      BOOST_CHECK_GT(std::abs(pos.z()), 0.125_um);
    }
  }
  for (const auto& mom : stepResultB.momentum) {
    if (mom == stepResultB.momentum.front()) {
      CHECK_CLOSE_ABS(mom, startMom, 1_keV);
    } else {
      BOOST_CHECK_NE(mom.x(), 5_GeV);
      CHECK_SMALL(mom.y(), 1_keV);
      BOOST_CHECK_NE(mom.z(), 0.);
    }
  }
}
// Test case c). Both should be involved in their part of the detector
BOOST_AUTO_TEST_CASE(step_extension_vacmatvac_test) {
  CuboidVolumeBuilder cvb;
  CuboidVolumeBuilder::VolumeConfig vConfVac1;
  vConfVac1.position = {0.5_m, 0., 0.};
  vConfVac1.length = {1_m, 1_m, 1_m};
  vConfVac1.name = "First vacuum volume";
  CuboidVolumeBuilder::VolumeConfig vConfMat;
  vConfMat.position = {1.5_m, 0., 0.};
  vConfMat.length = {1_m, 1_m, 1_m};
  vConfMat.volumeMaterial =
      std::make_shared<const HomogeneousVolumeMaterial>(makeBeryllium());
  vConfMat.name = "Material volume";
  CuboidVolumeBuilder::VolumeConfig vConfVac2;
  vConfVac2.position = {2.5_m, 0., 0.};
  vConfVac2.length = {1_m, 1_m, 1_m};
  vConfVac2.name = "Second vacuum volume";
  CuboidVolumeBuilder::Config conf;
  conf.volumeCfg = {vConfVac1, vConfMat, vConfVac2};
  conf.position = {1.5_m, 0., 0.};
  conf.length = {3_m, 1_m, 1_m};

  // Build detector
  cvb.setConfig(conf);
  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto& vb) {
        return cvb.trackingVolume(context, inner, vb);
      });
  TrackingGeometryBuilder tgb(tgbCfg);
  std::shared_ptr<const TrackingGeometry> det = tgb.trackingGeometry(tgContext);

  // Build navigator
  Navigator naviDet({det, true, true, true});

  // Set initial parameters for the particle track
  CurvilinearTrackParameters sbtp(Vector4::Zero(), 0_degree, 90_degree,
                                  1_e / 5_GeV, Covariance::Identity(),
                                  ParticleHypothesis::pion());

  using Stepper = EigenStepper<EigenStepperDenseExtension>;
  using Propagator = Acts::Propagator<Stepper, Navigator>;
  using PropagatorOptions =
      Propagator::Options<ActorList<StepCollector, EndOfWorld>>;

  // Set options for propagator
  PropagatorOptions propOpts(tgContext, mfContext);
  propOpts.actorList.get<EndOfWorld>().maxX = 3_m;
  propOpts.maxSteps = 1000;
  propOpts.stepping.maxStepSize = 1.5_m;

  // Build stepper and propagator
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 1_T, 0.));
  Stepper es(bField);
  Propagator prop(es, naviDet);

  // Launch and collect results
  const auto& result = prop.propagate(sbtp, propOpts).value();
  const StepCollector::this_result& stepResult =
      result.get<typename StepCollector::result_type>();

  // Manually set the extensions for each step and propagate through each
  // volume by propagation to the boundaries
  // Collect boundaries
  std::vector<Surface const*> surs;
  std::vector<std::shared_ptr<const BoundarySurfaceT<TrackingVolume>>>
      boundaries = det->lowestTrackingVolume(tgContext, {0.5_m, 0., 0.})
                       ->boundarySurfaces();
  for (auto& b : boundaries) {
    if (b->surfaceRepresentation().center(tgContext).x() == 1_m) {
      surs.push_back(&(b->surfaceRepresentation()));
      break;
    }
  }
  boundaries =
      det->lowestTrackingVolume(tgContext, {1.5_m, 0., 0.})->boundarySurfaces();
  for (auto& b : boundaries) {
    if (b->surfaceRepresentation().center(tgContext).x() == 2_m) {
      surs.push_back(&(b->surfaceRepresentation()));
      break;
    }
  }
  boundaries =
      det->lowestTrackingVolume(tgContext, {2.5_m, 0., 0.})->boundarySurfaces();
  for (auto& b : boundaries) {
    if (b->surfaceRepresentation().center(tgContext).x() == 3_m) {
      surs.push_back(&(b->surfaceRepresentation()));
      break;
    }
  }

  // Build launcher through vacuum
  // Set options for propagator

  using DefStepper = EigenStepper<EigenStepperDenseExtension>;
  using DefPropagator = Acts::Propagator<DefStepper, Navigator>;
  using DefPropagatorOptions =
      DefPropagator::Options<ActorList<StepCollector, EndOfWorld>>;

  DefPropagatorOptions propOptsDef(tgContext, mfContext);
  propOptsDef.actorList.get<EndOfWorld>().maxX = 3_m;
  propOptsDef.maxSteps = 1000;
  propOptsDef.stepping.maxStepSize = 1.5_m;

  // Build stepper and propagator
  DefStepper esDef(bField);
  DefPropagator propDef(esDef, naviDet);

  // Launch and collect results
  const auto& resultDef = propDef.propagate(sbtp, propOptsDef).value();
  const StepCollector::this_result& stepResultDef =
      resultDef.get<typename StepCollector::result_type>();

  // Check the exit situation of the first volume
  std::pair<Vector3, Vector3> endParams, endParamsControl;
  for (unsigned int i = 0; i < stepResultDef.position.size(); i++) {
    if (1_m - stepResultDef.position[i].x() < 1e-4) {
      endParams =
          std::make_pair(stepResultDef.position[i], stepResultDef.momentum[i]);
      break;
    }
  }
  for (unsigned int i = 0; i < stepResult.position.size(); i++) {
    if (1_m - stepResult.position[i].x() < 1e-4) {
      endParamsControl =
          std::make_pair(stepResult.position[i], stepResult.momentum[i]);
      break;
    }
  }

  CHECK_CLOSE_ABS(endParams.first, endParamsControl.first, 1_um);
  CHECK_CLOSE_ABS(endParams.second, endParamsControl.second, 1_um);

  CHECK_CLOSE_ABS(endParams.first.x(), endParamsControl.first.x(), 1e-5);
  CHECK_CLOSE_ABS(endParams.first.y(), endParamsControl.first.y(), 1e-5);
  CHECK_CLOSE_ABS(endParams.first.z(), endParamsControl.first.z(), 1e-5);
  CHECK_CLOSE_ABS(endParams.second.x(), endParamsControl.second.x(), 1e-5);
  CHECK_CLOSE_ABS(endParams.second.y(), endParamsControl.second.y(), 1e-5);
  CHECK_CLOSE_ABS(endParams.second.z(), endParamsControl.second.z(), 1e-5);

  // Build launcher through material
  // Set initial parameters for the particle track by using the result of the
  // first volume

  using DenseStepper = EigenStepper<EigenStepperDenseExtension>;
  using DensePropagator = Acts::Propagator<DenseStepper, Navigator>;
  using DensePropagatorOptions =
      DensePropagator::Options<ActorList<StepCollector, EndOfWorld>>;

  // Set options for propagator
  DensePropagatorOptions propOptsDense(tgContext, mfContext);
  propOptsDense.actorList.get<EndOfWorld>().maxX = 3_m;
  propOptsDense.maxSteps = 1000;
  propOptsDense.stepping.maxStepSize = 1.5_m;

  // Build stepper and propagator
  DenseStepper esDense(bField);
  DensePropagator propDense(esDense, naviDet);

  // Launch and collect results
  const auto& resultDense = propDense.propagate(sbtp, propOptsDense).value();
  const StepCollector::this_result& stepResultDense =
      resultDense.get<typename StepCollector::result_type>();

  // Check the exit situation of the second volume
  for (unsigned int i = 0; i < stepResultDense.position.size(); i++) {
    if (1_m - stepResultDense.position[i].x() < 1e-4) {
      endParams = std::make_pair(stepResultDense.position[i],
                                 stepResultDense.momentum[i]);
      break;
    }
  }
  for (unsigned int i = 0; i < stepResult.position.size(); i++) {
    if (1_m - stepResult.position[i].x() < 1e-4) {
      endParamsControl =
          std::make_pair(stepResult.position[i], stepResult.momentum[i]);
      break;
    }
  }

  CHECK_CLOSE_ABS(endParams.first, endParamsControl.first, 1_um);
  CHECK_CLOSE_ABS(endParams.second, endParamsControl.second, 1_um);
}

// Test case a). The DenseEnvironmentExtension should state that it is not
// valid in this case.
BOOST_AUTO_TEST_CASE(step_extension_trackercalomdt_test) {
  double rotationAngle = std::numbers::pi / 2.;
  Vector3 xPos(cos(rotationAngle), 0., sin(rotationAngle));
  Vector3 yPos(0., 1., 0.);
  Vector3 zPos(-sin(rotationAngle), 0., cos(rotationAngle));
  MaterialSlab matProp(makeBeryllium(), 0.5_mm);

  CuboidVolumeBuilder cvb;
  CuboidVolumeBuilder::SurfaceConfig sConf1;
  sConf1.position = Vector3(0.3_m, 0., 0.);
  sConf1.rotation.col(0) = xPos;
  sConf1.rotation.col(1) = yPos;
  sConf1.rotation.col(2) = zPos;
  sConf1.rBounds =
      std::make_shared<const RectangleBounds>(RectangleBounds(0.5_m, 0.5_m));
  sConf1.surMat = std::shared_ptr<const ISurfaceMaterial>(
      new HomogeneousSurfaceMaterial(matProp));
  sConf1.thickness = 1._mm;
  CuboidVolumeBuilder::LayerConfig lConf1;
  lConf1.surfaceCfg = {sConf1};

  CuboidVolumeBuilder::SurfaceConfig sConf2;
  sConf2.position = Vector3(0.6_m, 0., 0.);
  sConf2.rotation.col(0) = xPos;
  sConf2.rotation.col(1) = yPos;
  sConf2.rotation.col(2) = zPos;
  sConf2.rBounds =
      std::make_shared<const RectangleBounds>(RectangleBounds(0.5_m, 0.5_m));
  sConf2.surMat = std::shared_ptr<const ISurfaceMaterial>(
      new HomogeneousSurfaceMaterial(matProp));
  sConf2.thickness = 1._mm;
  CuboidVolumeBuilder::LayerConfig lConf2;
  lConf2.surfaceCfg = {sConf2};

  CuboidVolumeBuilder::VolumeConfig muConf1;
  muConf1.position = {2.3_m, 0., 0.};
  muConf1.length = {20._cm, 20._cm, 20._cm};
  muConf1.volumeMaterial =
      std::make_shared<HomogeneousVolumeMaterial>(makeBeryllium());
  muConf1.name = "MDT1";
  CuboidVolumeBuilder::VolumeConfig muConf2;
  muConf2.position = {2.7_m, 0., 0.};
  muConf2.length = {20._cm, 20._cm, 20._cm};
  muConf2.volumeMaterial =
      std::make_shared<HomogeneousVolumeMaterial>(makeBeryllium());
  muConf2.name = "MDT2";

  CuboidVolumeBuilder::VolumeConfig vConf1;
  vConf1.position = {0.5_m, 0., 0.};
  vConf1.length = {1._m, 1._m, 1._m};
  vConf1.layerCfg = {lConf1, lConf2};
  vConf1.name = "Tracker";
  CuboidVolumeBuilder::VolumeConfig vConf2;
  vConf2.position = {1.5_m, 0., 0.};
  vConf2.length = {1._m, 1._m, 1._m};
  vConf2.volumeMaterial =
      std::make_shared<HomogeneousVolumeMaterial>(makeBeryllium());
  vConf2.name = "Calorimeter";
  CuboidVolumeBuilder::VolumeConfig vConf3;
  vConf3.position = {2.5_m, 0., 0.};
  vConf3.length = {1._m, 1._m, 1._m};
  vConf3.volumeCfg = {muConf1, muConf2};
  vConf3.name = "Muon system";
  CuboidVolumeBuilder::Config conf;
  conf.volumeCfg = {vConf1, vConf2, vConf3};
  conf.position = {1.5_m, 0., 0.};
  conf.length = {3._m, 1._m, 1._m};

  // Build detector
  cvb.setConfig(conf);
  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto& vb) {
        return cvb.trackingVolume(context, inner, vb);
      });
  TrackingGeometryBuilder tgb(tgbCfg);
  std::shared_ptr<const TrackingGeometry> detector =
      tgb.trackingGeometry(tgContext);

  // Build navigator
  Navigator naviVac({detector, true, true, true});

  // Set initial parameters for the particle track
  CurvilinearTrackParameters sbtp(Vector4::Zero(), 0_degree, 90_degree,
                                  1_e / 1_GeV, Covariance::Identity(),
                                  ParticleHypothesis::pion());

  using Stepper = EigenStepper<EigenStepperDenseExtension>;
  using Propagator = Acts::Propagator<Stepper, Navigator>;
  using PropagatorOptions = Propagator::Options<
      ActorList<StepCollector, MaterialInteractor, EndOfWorld>>;

  // Set options for propagator
  PropagatorOptions propOpts(tgContext, mfContext);
  propOpts.actorList.get<EndOfWorld>().maxX = 3._m;
  propOpts.maxSteps = 10000;

  // Build stepper and propagator
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 0.));
  Stepper es(bField);
  Propagator prop(es, naviVac);

  // Launch and collect results
  const auto& result = prop.propagate(sbtp, propOpts).value();
  const StepCollector::this_result& stepResult =
      result.get<typename StepCollector::result_type>();

  // Test that momentum changes only occurred at the right detector parts
  double lastMomentum = stepResult.momentum[0].x();
  for (unsigned int i = 0; i < stepResult.position.size(); i++) {
    // Test for changes
    if ((stepResult.position[i].x() > 0.3_m &&
         stepResult.position[i].x() < 0.6_m) ||
        (stepResult.position[i].x() > 0.6_m &&
         stepResult.position[i].x() <= 1._m) ||
        (stepResult.position[i].x() > 1._m &&
         stepResult.position[i].x() <= 2._m) ||
        (stepResult.position[i].x() > 2.2_m &&
         stepResult.position[i].x() <= 2.4_m) ||
        (stepResult.position[i].x() > 2.6_m &&
         stepResult.position[i].x() <= 2.8_m)) {
      BOOST_CHECK_LE(stepResult.momentum[i].x(), lastMomentum);
      lastMomentum = stepResult.momentum[i].x();
    } else
    // Test the absence of momentum loss
    {
      if (stepResult.position[i].x() < 0.3_m ||
          (stepResult.position[i].x() > 2._m &&
           stepResult.position[i].x() <= 2.2_m) ||
          (stepResult.position[i].x() > 2.4_m &&
           stepResult.position[i].x() <= 2.6_m) ||
          (stepResult.position[i].x() > 2.8_m &&
           stepResult.position[i].x() <= 3._m)) {
        BOOST_CHECK_EQUAL(stepResult.momentum[i].x(), lastMomentum);
      }
    }
  }
}

}  // namespace Acts::Test
