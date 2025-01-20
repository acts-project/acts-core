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
#include "Acts/EventData/MultiComponentTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/EigenStepperDefaultExtension.hpp"
#include "Acts/Propagator/MultiEigenStepperLoop.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <numbers>
#include <optional>
#include <random>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace Acts {
namespace Concepts {
template <typename T>
concept has_components = requires { typename T::components; };
}  // namespace Concepts

struct MultiStepperSurfaceReached;
}  // namespace Acts

using namespace Acts;
using namespace Acts::VectorHelpers;

/////////////////////////////////////////////////////
// Some useful global objects, typedefs and structs
/////////////////////////////////////////////////////
const MagneticFieldContext magCtx;
const GeometryContext geoCtx;

using MultiStepperLoop = MultiEigenStepperLoop<EigenStepperDefaultExtension>;
using SingleStepper = EigenStepper<EigenStepperDefaultExtension>;

const double defaultStepSize = 123.;
const auto defaultNDir = Direction::Backward();

const auto defaultBField =
    std::make_shared<ConstantBField>(Vector3(1., 2.5, 33.33));
const auto defaultNullBField = std::make_shared<NullBField>();

const auto particleHypothesis = ParticleHypothesis::pion();

// Makes random bound parameters and covariance and a plane surface at {0,0,0}
// with normal {1,0,0}. Optionally some external fixed bound parameters can be
// supplied
auto makeDefaultBoundPars(bool cov = true, std::size_t n = 4,
                          std::optional<BoundVector> ext_pars = std::nullopt) {
  std::vector<std::tuple<double, BoundVector, std::optional<BoundSquareMatrix>>>
      cmps;
  using Opt = std::optional<BoundSquareMatrix>;

  auto make_random_sym_matrix = []() {
    auto c = BoundSquareMatrix::Random().eval();
    c *= c.transpose();
    return c;
  };

  // note that we are using the default random device
  std::mt19937 gen;
  std::uniform_real_distribution<> locDis(-10., 10.);
  std::uniform_real_distribution<> phiDis(-std::numbers::pi, std::numbers::pi);
  std::uniform_real_distribution<> thetaDis(0., std::numbers::pi);
  std::uniform_real_distribution<> qOverPDis(-10., 10.);
  std::uniform_real_distribution<> timeDis(0., 100.);

  for (auto i = 0ul; i < n; ++i) {
    BoundVector params = BoundVector::Zero();

    if (ext_pars) {
      params = *ext_pars;
    } else {
      params[eBoundLoc0] = locDis(gen);
      params[eBoundLoc1] = locDis(gen);
      params[eBoundPhi] = phiDis(gen);
      params[eBoundTheta] = thetaDis(gen);
      params[eBoundQOverP] = qOverPDis(gen);
      params[eBoundTime] = timeDis(gen);
    }

    cmps.push_back(
        {1. / n, params, cov ? Opt{make_random_sym_matrix()} : Opt{}});
  }

  auto surface = Acts::CurvilinearSurface(Vector3::Zero(), Vector3{1., 0., 0.})
                     .planeSurface();

  return MultiComponentBoundTrackParameters(surface, cmps, particleHypothesis);
}

//////////////////////
/// Test the reducers
//////////////////////
BOOST_AUTO_TEST_CASE(test_max_weight_reducer) {
  // Can use this multistepper since we only care about the state which is
  // invariant
  using MultiOptions = MultiStepperLoop::Options;
  using MultiState = MultiStepperLoop::State;

  MultiOptions options(geoCtx, magCtx);
  options.maxStepSize = defaultStepSize;

  SingleStepper singleStepper(defaultBField);
  MultiStepperLoop multiStepper(defaultBField);

  constexpr std::size_t N = 4;
  const auto multi_pars = makeDefaultBoundPars(false, N);
  MultiState state = multiStepper.makeState(options, multi_pars);

  double w = 0.1;
  double wSum = 0.0;
  for (auto &[sstate, weight, _] : state.components) {
    weight = w;
    wSum += w;
    w += 0.1;
  }
  BOOST_CHECK_EQUAL(wSum, 1.0);
  BOOST_CHECK_EQUAL(state.components.back().weight, 0.4);

  MaxWeightReducerLoop reducer{};
  BOOST_CHECK_EQUAL(reducer.position(state),
                    singleStepper.position(state.components.back().state));
  BOOST_CHECK_EQUAL(reducer.direction(state),
                    singleStepper.direction(state.components.back().state));
}

BOOST_AUTO_TEST_CASE(test_max_momentum_reducer) {
  // Can use this multistepper since we only care about the state which is
  // invariant
  using MultiOptions = MultiStepperLoop::Options;
  using MultiState = typename MultiStepperLoop::State;

  MultiOptions options(geoCtx, magCtx);
  options.maxStepSize = defaultStepSize;

  SingleStepper singleStepper(defaultBField);
  MultiStepperLoop multiStepper(defaultBField);

  constexpr std::size_t N = 4;
  const auto multi_pars = makeDefaultBoundPars(false, N);
  MultiState state = multiStepper.makeState(options, multi_pars);

  double p = 1.0;
  double q = 1.0;
  for (auto &[sstate, weight, _] : state.components) {
    sstate.pars[eFreeQOverP] = q / p;
    p *= 2.0;
  }
  BOOST_CHECK_EQUAL(state.components.back().state.pars[eFreeQOverP], q / 8.0);

  MaxMomentumReducerLoop reducer{};
  BOOST_CHECK_EQUAL(reducer.position(state),
                    singleStepper.position(state.components.back().state));
  BOOST_CHECK_EQUAL(reducer.direction(state),
                    singleStepper.direction(state.components.back().state));
}

//////////////////////////////////////////////////////
/// Test the construction of the MultiStepper::State
//////////////////////////////////////////////////////
template <typename multi_stepper_t, bool Cov>
void test_multi_stepper_state() {
  using MultiOptions = typename multi_stepper_t::Options;
  using MultiState = typename multi_stepper_t::State;
  using MultiStepper = multi_stepper_t;

  MultiOptions options(geoCtx, magCtx);
  options.maxStepSize = defaultStepSize;

  SingleStepper singleStepper(defaultBField);
  MultiStepper multiStepper(defaultBField);

  constexpr std::size_t N = 4;
  const auto multi_pars = makeDefaultBoundPars(Cov, N, BoundVector::Ones());

  MultiState state = multiStepper.makeState(options, multi_pars);

  BOOST_CHECK_EQUAL(N, multiStepper.numberComponents(state));

  // Test the result & compare with the input/test for reasonable members
  auto const_iterable = multiStepper.constComponentIterable(state);
  for (const auto cmp : const_iterable) {
    BOOST_CHECK_EQUAL(cmp.jacTransport(), FreeMatrix::Identity());
    BOOST_CHECK_EQUAL(cmp.derivative(), FreeVector::Zero());
    if constexpr (!Cov) {
      BOOST_CHECK_EQUAL(cmp.jacToGlobal(), BoundToFreeMatrix::Zero());
      BOOST_CHECK_EQUAL(cmp.cov(), BoundSquareMatrix::Zero());
    }
  }

  BOOST_CHECK_EQUAL(state.pathAccumulated, 0.);
  for (const auto cmp : const_iterable) {
    BOOST_CHECK_EQUAL(cmp.pathAccumulated(), 0.);
  }

  // covTransport in the MultiEigenStepperLoop is redundant and
  // thus not part of the interface. However, we want to check them for
  // consistency.
  if constexpr (Acts::Concepts::has_components<MultiState>) {
    BOOST_CHECK(!state.covTransport);
    for (const auto &cmp : state.components) {
      BOOST_CHECK_EQUAL(cmp.state.covTransport, Cov);
    }
  }
}

BOOST_AUTO_TEST_CASE(multi_stepper_state_no_cov) {
  test_multi_stepper_state<MultiStepperLoop, false>();
}

template <typename multi_stepper_t>
void test_multi_stepper_state_invalid() {
  using MultiOptions = typename multi_stepper_t::Options;

  MultiOptions options(geoCtx, magCtx);
  options.maxStepSize = defaultStepSize;

  MultiStepperLoop multi_stepper(defaultBField);

  // Empty component vector
  const auto multi_pars = makeDefaultBoundPars(false, 0);

  BOOST_CHECK_THROW(multi_stepper.makeState(options, multi_pars),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(multi_eigen_stepper_state_invalid) {
  test_multi_stepper_state_invalid<MultiStepperLoop>();
}

////////////////////////////////////////////////////////////////////////
// Compare the Multi-Stepper against the Eigen-Stepper for consistency
////////////////////////////////////////////////////////////////////////
template <typename multi_stepper_t>
void test_multi_stepper_vs_eigen_stepper() {
  using MultiOptions = typename multi_stepper_t::Options;
  using MultiState = typename multi_stepper_t::State;
  using MultiStepper = multi_stepper_t;

  MultiOptions options(geoCtx, magCtx);
  options.maxStepSize = defaultStepSize;

  MultiStepper multi_stepper(defaultBField);
  SingleStepper single_stepper(defaultBField);

  const BoundVector pars = BoundVector::Ones();
  const BoundSquareMatrix cov = BoundSquareMatrix::Identity();

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSquareMatrix>>>
      cmps(4, {0.25, pars, cov});

  auto surface =
      Acts::CurvilinearSurface(Vector3::Zero(), Vector3::Ones().normalized())
          .planeSurface();

  MultiComponentBoundTrackParameters multi_pars(surface, cmps,
                                                particleHypothesis);
  BoundTrackParameters single_pars(surface, pars, cov, particleHypothesis);

  MultiState multi_state = multi_stepper.makeState(options, multi_pars);
  SingleStepper::State single_state =
      single_stepper.makeState(options, single_pars);

  for (auto cmp : multi_stepper.componentIterable(multi_state)) {
    cmp.status() = Acts::IntersectionStatus::reachable;
  }

  // Do some steps and check that the results match
  for (int i = 0; i < 10; ++i) {
    // Single stepper
    auto single_result =
        single_stepper.step(single_state, defaultNDir, nullptr);
    single_stepper.transportCovarianceToCurvilinear(single_state);

    // Multi stepper;
    auto multi_result = multi_stepper.step(multi_state, defaultNDir, nullptr);
    multi_stepper.transportCovarianceToCurvilinear(multi_state);

    // Check equality
    BOOST_REQUIRE(multi_result.ok());
    BOOST_REQUIRE_EQUAL(multi_result.ok(), single_result.ok());

    BOOST_CHECK_EQUAL(*single_result, *multi_result);

    for (const auto cmp : multi_stepper.constComponentIterable(multi_state)) {
      BOOST_CHECK_EQUAL(cmp.pars(), single_state.pars);
      BOOST_CHECK_EQUAL(cmp.cov(), single_state.cov);
      BOOST_CHECK_EQUAL(cmp.jacTransport(), single_state.jacTransport);
      BOOST_CHECK_EQUAL(cmp.jacToGlobal(), single_state.jacToGlobal);
      BOOST_CHECK_EQUAL(cmp.derivative(), single_state.derivative);
      BOOST_CHECK_EQUAL(cmp.pathAccumulated(), single_state.pathAccumulated);
    }
  }
}

BOOST_AUTO_TEST_CASE(multi_eigen_vs_single_eigen) {
  test_multi_stepper_vs_eigen_stepper<MultiStepperLoop>();
}

/////////////////////////////
// Test stepsize accessors
/////////////////////////////

// TODO do this later, when we introduce the MultiEigenStepperSIMD, which there
// needs new interfaces...

////////////////////////////////////////////////////
// Test the modifying accessors to the components
////////////////////////////////////////////////////
template <typename multi_stepper_t>
void test_components_modifying_accessors() {
  using MultiOptions = typename multi_stepper_t::Options;
  using MultiState = typename multi_stepper_t::State;
  using MultiStepper = multi_stepper_t;

  MultiOptions options(geoCtx, magCtx);
  options.maxStepSize = defaultStepSize;

  const auto multi_pars = makeDefaultBoundPars();

  MultiStepper multi_stepper(defaultBField);

  MultiState mutable_multi_state = multi_stepper.makeState(options, multi_pars);
  const MultiState const_multi_state =
      multi_stepper.makeState(options, multi_pars);

  auto modify = [&](const auto &projector) {
    // Here test the mutable overloads of the mutable iterable
    for (auto cmp : multi_stepper.componentIterable(mutable_multi_state)) {
      using type = std::decay_t<decltype(projector(cmp))>;
      if constexpr (std::is_enum_v<type>) {
        projector(cmp) =
            static_cast<type>(static_cast<int>(projector(cmp)) + 1);
      } else {
        projector(cmp) *= 2.0;
      }
    }
  };

  auto check = [&](const auto &projector) {
    // Here test the const-member functions of the mutable iterable
    auto mutable_state_iterable =
        multi_stepper.componentIterable(mutable_multi_state);
    // Here test the const iterable
    auto const_state_iterable =
        multi_stepper.constComponentIterable(const_multi_state);

    auto mstate_it = mutable_state_iterable.begin();
    auto cstate_it = const_state_iterable.begin();
    for (; cstate_it != const_state_iterable.end(); ++mstate_it, ++cstate_it) {
      const auto mstate_cmp = *mstate_it;
      auto cstate_cmp = *cstate_it;

      using type = std::decay_t<decltype(projector(mstate_cmp))>;

      if constexpr (std::is_arithmetic_v<type>) {
        BOOST_CHECK_CLOSE(projector(mstate_cmp), 2.0 * projector(cstate_cmp),
                          1.e-8);
      } else if constexpr (std::is_enum_v<type>) {
        BOOST_CHECK_EQUAL(static_cast<int>(projector(mstate_cmp)),
                          1 + static_cast<int>(projector(cstate_cmp)));
      } else {
        BOOST_CHECK(
            projector(mstate_cmp).isApprox(2.0 * projector(cstate_cmp), 1.e-8));
      }
    }
  };

  const auto projectors = std::make_tuple(
      [](auto &cmp) -> decltype(auto) { return cmp.status(); },
      [](auto &cmp) -> decltype(auto) { return cmp.pathAccumulated(); },
      [](auto &cmp) -> decltype(auto) { return cmp.weight(); },
      [](auto &cmp) -> decltype(auto) { return cmp.pars(); },
      [](auto &cmp) -> decltype(auto) { return cmp.cov(); },
      [](auto &cmp) -> decltype(auto) { return cmp.jacTransport(); },
      [](auto &cmp) -> decltype(auto) { return cmp.derivative(); },
      [](auto &cmp) -> decltype(auto) { return cmp.jacobian(); },
      [](auto &cmp) -> decltype(auto) { return cmp.jacToGlobal(); });

  std::apply(
      [&](const auto &...projs) {
        // clang-format off
        ( [&]() { modify(projs); check(projs); }(), ...);
        // clang-format on
      },
      projectors);
}

BOOST_AUTO_TEST_CASE(multi_eigen_component_iterable_with_modification) {
  test_components_modifying_accessors<MultiStepperLoop>();
}

/////////////////////////////////////////////
// Test if the surface status update works
/////////////////////////////////////////////
template <typename multi_stepper_t>
void test_multi_stepper_surface_status_update() {
  using MultiOptions = typename multi_stepper_t::Options;
  using MultiState = typename multi_stepper_t::State;
  using MultiStepper = multi_stepper_t;

  MultiOptions options(geoCtx, magCtx);
  options.maxStepSize = defaultStepSize;

  MultiStepper multi_stepper(defaultNullBField);
  SingleStepper single_stepper(defaultNullBField);

  auto start_surface =
      Acts::CurvilinearSurface(Vector3::Zero(), Vector3{1.0, 0.0, 0.0})
          .planeSurface();

  auto right_surface =
      Acts::CurvilinearSurface(Vector3{1.0, 0.0, 0.0}, Vector3{1.0, 0.0, 0.0})
          .planeSurface();

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSquareMatrix>>>
      cmps(2, {0.5, BoundVector::Zero(), std::nullopt});
  std::get<BoundVector>(cmps[0])[eBoundTheta] = std::numbers::pi / 2.;
  std::get<BoundVector>(cmps[1])[eBoundPhi] = std::numbers::pi;
  std::get<BoundVector>(cmps[1])[eBoundTheta] = std::numbers::pi / 2.;
  std::get<BoundVector>(cmps[0])[eBoundQOverP] = 1.0;
  std::get<BoundVector>(cmps[1])[eBoundQOverP] = 1.0;

  MultiComponentBoundTrackParameters multi_pars(start_surface, cmps,
                                                particleHypothesis);

  BOOST_REQUIRE(std::get<1>(multi_pars[0])
                    .direction()
                    .isApprox(Vector3{1.0, 0.0, 0.0}, 1.e-10));
  BOOST_REQUIRE(std::get<1>(multi_pars[1])
                    .direction()
                    .isApprox(Vector3{-1.0, 0.0, 0.0}, 1.e-10));

  MultiState multi_state = multi_stepper.makeState(options, multi_pars);
  SingleStepper::State single_state =
      single_stepper.makeState(options, std::get<1>(multi_pars[0]));

  // Update surface status and check
  {
    auto status = multi_stepper.updateSurfaceStatus(
        multi_state, *right_surface, 0, Direction::Forward(),
        BoundaryTolerance::Infinite(), s_onSurfaceTolerance,
        ConstrainedStep::Type::Navigator);

    BOOST_CHECK_EQUAL(status, IntersectionStatus::reachable);

    auto cmp_iterable = multi_stepper.constComponentIterable(multi_state);
    auto cmp_1 = *cmp_iterable.begin();
    auto cmp_2 = *(++cmp_iterable.begin());

    BOOST_CHECK_EQUAL(cmp_1.status(), IntersectionStatus::reachable);
    BOOST_CHECK_EQUAL(cmp_2.status(), IntersectionStatus::reachable);

    BOOST_CHECK_EQUAL(cmp_1.cmp.state.stepSize.value(), 1.0);
    BOOST_CHECK_EQUAL(cmp_2.cmp.state.stepSize.value(), -1.0);
  }

  // Step forward now
  {
    multi_stepper.step(multi_state, Direction::Forward(), nullptr);

    // Single stepper
    single_stepper.step(single_state, Direction::Forward(), nullptr);
  }

  // Update surface status and check again
  {
    auto status = multi_stepper.updateSurfaceStatus(
        multi_state, *right_surface, 0, Direction::Forward(),
        BoundaryTolerance::Infinite(), s_onSurfaceTolerance,
        ConstrainedStep::Type::Navigator);

    BOOST_CHECK_EQUAL(status, IntersectionStatus::onSurface);

    auto cmp_iterable = multi_stepper.constComponentIterable(multi_state);
    auto cmp_1 = *cmp_iterable.begin();
    auto cmp_2 = *(++cmp_iterable.begin());

    BOOST_CHECK_EQUAL(cmp_1.status(), IntersectionStatus::onSurface);
    BOOST_CHECK_EQUAL(cmp_2.status(), IntersectionStatus::onSurface);
  }

  // Start surface should be reachable
  {
    auto status = multi_stepper.updateSurfaceStatus(
        multi_state, *start_surface, 0, Direction::Forward(),
        BoundaryTolerance::Infinite(), s_onSurfaceTolerance,
        ConstrainedStep::Type::Navigator);

    BOOST_CHECK_EQUAL(status, IntersectionStatus::reachable);

    auto cmp_iterable = multi_stepper.constComponentIterable(multi_state);
    auto cmp_1 = *cmp_iterable.begin();
    auto cmp_2 = *(++cmp_iterable.begin());

    BOOST_CHECK_EQUAL(cmp_1.status(), IntersectionStatus::reachable);
    BOOST_CHECK_EQUAL(cmp_2.status(), IntersectionStatus::reachable);

    BOOST_CHECK_EQUAL(cmp_1.cmp.state.stepSize.value(), -1.0);
    BOOST_CHECK_EQUAL(cmp_2.cmp.state.stepSize.value(), 1.0);
  }
}

BOOST_AUTO_TEST_CASE(test_surface_status_and_cmpwise_bound_state) {
  test_multi_stepper_surface_status_update<MultiStepperLoop>();
}

//////////////////////////////////
// Test Bound state computations
//////////////////////////////////
template <typename multi_stepper_t>
void test_component_bound_state() {
  using MultiOptions = typename multi_stepper_t::Options;
  using MultiState = typename multi_stepper_t::State;
  using MultiStepper = multi_stepper_t;

  MultiOptions options(geoCtx, magCtx);
  options.maxStepSize = defaultStepSize;

  MultiStepper multi_stepper(defaultNullBField);
  SingleStepper single_stepper(defaultNullBField);

  auto start_surface =
      Acts::CurvilinearSurface(Vector3::Zero(), Vector3{1.0, 0.0, 0.0})
          .planeSurface();

  auto right_surface =
      Acts::CurvilinearSurface(Vector3{1.0, 0.0, 0.0}, Vector3{1.0, 0.0, 0.0})
          .planeSurface();

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSquareMatrix>>>
      cmps(2, {0.5, BoundVector::Zero(), std::nullopt});
  std::get<BoundVector>(cmps[0])[eBoundTheta] = std::numbers::pi / 2.;
  std::get<BoundVector>(cmps[1])[eBoundPhi] = std::numbers::pi;
  std::get<BoundVector>(cmps[1])[eBoundTheta] = std::numbers::pi / 2.;
  std::get<BoundVector>(cmps[0])[eBoundQOverP] = 1.0;
  std::get<BoundVector>(cmps[1])[eBoundQOverP] = 1.0;

  MultiComponentBoundTrackParameters multi_pars(start_surface, cmps,
                                                particleHypothesis);

  BOOST_REQUIRE(std::get<1>(multi_pars[0])
                    .direction()
                    .isApprox(Vector3{1.0, 0.0, 0.0}, 1.e-10));
  BOOST_REQUIRE(std::get<1>(multi_pars[1])
                    .direction()
                    .isApprox(Vector3{-1.0, 0.0, 0.0}, 1.e-10));

  MultiState multi_state = multi_stepper.makeState(options, multi_pars);
  SingleStepper::State single_state =
      single_stepper.makeState(options, std::get<1>(multi_pars[0]));

  // Step forward now
  {
    multi_stepper.updateSurfaceStatus(
        multi_state, *right_surface, 0, Direction::Forward(),
        BoundaryTolerance::Infinite(), s_onSurfaceTolerance,
        ConstrainedStep::Type::Navigator);
    multi_stepper.step(multi_state, Direction::Forward(), nullptr);

    // Single stepper
    single_stepper.updateSurfaceStatus(
        single_state, *right_surface, 0, Direction::Forward(),
        BoundaryTolerance::Infinite(), s_onSurfaceTolerance,
        ConstrainedStep::Type::Navigator);
    single_stepper.step(single_state, Direction::Forward(), nullptr);
  }

  // Check component-wise bound-state
  {
    auto single_bound_state = single_stepper.boundState(
        single_state, *right_surface, true, FreeToBoundCorrection(false));
    BOOST_REQUIRE(single_bound_state.ok());

    auto cmp_iterable = multi_stepper.componentIterable(multi_state);
    auto cmp_1 = *cmp_iterable.begin();
    auto cmp_2 = *(++cmp_iterable.begin());

    auto bound_state_1 =
        cmp_1.boundState(*right_surface, true, FreeToBoundCorrection(false));
    BOOST_REQUIRE(bound_state_1.ok());
    BOOST_CHECK(*single_bound_state == *bound_state_1);

    auto bound_state_2 =
        cmp_2.boundState(*right_surface, true, FreeToBoundCorrection(false));
    BOOST_CHECK(bound_state_2.ok());
  }
}

BOOST_AUTO_TEST_CASE(test_component_wise_bound_state) {
  test_component_bound_state<MultiStepperLoop>();
}

template <typename multi_stepper_t>
void test_combined_bound_state_function() {
  using MultiOptions = typename multi_stepper_t::Options;
  using MultiState = typename multi_stepper_t::State;
  using MultiStepper = multi_stepper_t;

  MultiOptions options(geoCtx, magCtx);
  options.maxStepSize = defaultStepSize;

  MultiStepper multi_stepper(defaultBField);

  auto surface =
      Acts::CurvilinearSurface(Vector3::Zero(), Vector3{1.0, 0.0, 0.0})
          .planeSurface();

  // Use Ones() here, so that the angles are in correct range
  const auto pars = BoundVector::Ones().eval();
  const auto cov = []() {
    auto c = BoundSquareMatrix::Random().eval();
    c *= c.transpose();
    return c;
  }();

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSquareMatrix>>>
      cmps(4, {0.25, pars, cov});

  MultiComponentBoundTrackParameters multi_pars(surface, cmps,
                                                particleHypothesis);
  MultiState multi_state = multi_stepper.makeState(options, multi_pars);

  auto res = multi_stepper.boundState(multi_state, *surface, true,
                                      FreeToBoundCorrection(false));

  BOOST_REQUIRE(res.ok());

  const auto [bound_pars, jacobian, pathLength] = *res;

  BOOST_CHECK_EQUAL(jacobian, decltype(jacobian)::Zero());
  BOOST_CHECK_EQUAL(pathLength, 0.0);
  BOOST_CHECK(bound_pars.parameters().isApprox(pars, 1.e-8));
  BOOST_CHECK(bound_pars.covariance()->isApprox(cov, 1.e-8));
}

BOOST_AUTO_TEST_CASE(test_combined_bound_state) {
  test_combined_bound_state_function<MultiStepperLoop>();
}

//////////////////////////////////////////////////
// Test the combined curvilinear state function
//////////////////////////////////////////////////
template <typename multi_stepper_t>
void test_combined_curvilinear_state_function() {
  using MultiOptions = typename multi_stepper_t::Options;
  using MultiState = typename multi_stepper_t::State;
  using MultiStepper = multi_stepper_t;

  MultiOptions options(geoCtx, magCtx);
  options.maxStepSize = defaultStepSize;

  MultiStepper multi_stepper(defaultBField);

  auto surface =
      Acts::CurvilinearSurface(Vector3::Zero(), Vector3{1.0, 0.0, 0.0})
          .planeSurface();

  // Use Ones() here, so that the angles are in correct range
  const auto pars = BoundVector::Ones().eval();
  const auto cov = []() {
    auto c = BoundSquareMatrix::Random().eval();
    c *= c.transpose();
    return c;
  }();

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSquareMatrix>>>
      cmps(4, {0.25, pars, cov});
  BoundTrackParameters check_pars(surface, pars, cov, particleHypothesis);

  MultiComponentBoundTrackParameters multi_pars(surface, cmps,
                                                particleHypothesis);
  MultiState multi_state = multi_stepper.makeState(options, multi_pars);

  const auto [curv_pars, jac, pathLength] =
      multi_stepper.curvilinearState(multi_state);

  BOOST_CHECK(curv_pars.fourPosition(geoCtx).isApprox(
      check_pars.fourPosition(geoCtx), 1.e-8));
  BOOST_CHECK(curv_pars.direction().isApprox(check_pars.direction(), 1.e-8));
  BOOST_CHECK_CLOSE(curv_pars.absoluteMomentum(), check_pars.absoluteMomentum(),
                    1.e-8);
  BOOST_CHECK_CLOSE(curv_pars.charge(), check_pars.charge(), 1.e-8);
}

BOOST_AUTO_TEST_CASE(test_curvilinear_state) {
  test_combined_curvilinear_state_function<MultiStepperLoop>();
}

////////////////////////////////////
// Test single component interface
////////////////////////////////////

template <typename multi_stepper_t>
void test_single_component_interface_function() {
  using MultiOptions = typename multi_stepper_t::Options;
  using MultiState = typename multi_stepper_t::State;
  using MultiStepper = multi_stepper_t;

  MultiOptions options(geoCtx, magCtx);
  options.maxStepSize = defaultStepSize;

  MultiStepper multi_stepper(defaultBField);

  MultiComponentBoundTrackParameters multi_pars = makeDefaultBoundPars(true, 4);

  MultiState multi_state = multi_stepper.makeState(options, multi_pars);

  // Check at least some properties at the moment
  auto check = [&](auto cmp) {
    auto sstepper = cmp.singleStepper(multi_stepper);
    auto &sstepping = cmp.state();

    BOOST_CHECK_EQUAL(sstepper.position(sstepping),
                      cmp.pars().template segment<3>(eFreePos0));
    BOOST_CHECK_EQUAL(sstepper.direction(sstepping),
                      cmp.pars().template segment<3>(eFreeDir0));
    BOOST_CHECK_EQUAL(sstepper.time(sstepping), cmp.pars()[eFreeTime]);
    BOOST_CHECK_CLOSE(sstepper.qOverP(sstepping), cmp.pars()[eFreeQOverP],
                      1.e-8);
  };

  for (const auto cmp : multi_stepper.constComponentIterable(multi_state)) {
    check(cmp);
  }

  for (auto cmp : multi_stepper.componentIterable(multi_state)) {
    check(cmp);
  }
}

BOOST_AUTO_TEST_CASE(test_single_component_interface) {
  test_single_component_interface_function<MultiStepperLoop>();
}

//////////////////////////////
// Remove and add components
//////////////////////////////

template <typename multi_stepper_t>
void remove_add_components_function() {
  using MultiOptions = typename multi_stepper_t::Options;
  using MultiState = typename multi_stepper_t::State;
  using MultiStepper = multi_stepper_t;

  MultiOptions options(geoCtx, magCtx);

  MultiStepper multi_stepper(defaultBField);

  const auto multi_pars = makeDefaultBoundPars(4);

  MultiState multi_state = multi_stepper.makeState(options, multi_pars);

  {
    BoundTrackParameters pars(multi_pars.referenceSurface().getSharedPtr(),
                              BoundVector::Ones(), std::nullopt,
                              particleHypothesis);
    multi_stepper.addComponent(multi_state, pars, 0.0);
  }

  BOOST_CHECK_EQUAL(multi_stepper.numberComponents(multi_state),
                    multi_pars.components().size() + 1);

  multi_stepper.clearComponents(multi_state);

  BOOST_CHECK_EQUAL(multi_stepper.numberComponents(multi_state), 0);
}

BOOST_AUTO_TEST_CASE(remove_add_components_test) {
  remove_add_components_function<MultiStepperLoop>();
}

//////////////////////////////////////////////////
// Instantiate a Propagator with the MultiStepper
//////////////////////////////////////////////////

template <typename multi_stepper_t>
void propagator_instatiation_test_function() {
  auto bField = std::make_shared<NullBField>();
  multi_stepper_t multi_stepper(bField);
  Propagator<multi_stepper_t, Navigator> propagator(
      std::move(multi_stepper), Navigator{Navigator::Config{}});

  auto surface =
      Acts::CurvilinearSurface(Vector3::Zero(), Vector3{1.0, 0.0, 0.0})
          .planeSurface();
  using PropagatorOptions =
      typename Propagator<multi_stepper_t, Navigator>::template Options<>;
  PropagatorOptions options(geoCtx, magCtx);

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSquareMatrix>>>
      cmps(4, {0.25, BoundVector::Ones().eval(),
               BoundSquareMatrix::Identity().eval()});
  MultiComponentBoundTrackParameters pars(surface, cmps, particleHypothesis);

  // This only checks that this compiles, not that it runs without errors
  // @TODO: Add test that checks the target aborter works correctly

  // Instantiate with target
  using type_a =
      decltype(propagator.template propagate<decltype(pars), decltype(options),
                                             MultiStepperSurfaceReached>(
          pars, *surface, options));
  static_assert(!std::is_same_v<type_a, void>);

  // Instantiate without target
  using type_b = decltype(propagator.propagate(pars, options));
  static_assert(!std::is_same_v<type_b, void>);
}

BOOST_AUTO_TEST_CASE(propagator_instatiation_test) {
  propagator_instatiation_test_function<MultiStepperLoop>();
}
