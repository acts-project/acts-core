// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/PropagatorTraits.hpp"
#include "Acts/Propagator/StepperOptions.hpp"
#include "Acts/Propagator/StepperStatistics.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cmath>
#include <limits>
#include <string>
#include <tuple>

namespace Acts {

class IVolumeMaterial;

/// @brief straight line stepper based on Surface intersection
///
/// The straight line stepper is a simple navigation stepper
/// to be used to navigate through the tracking geometry. It can be
/// used for simple material mapping, navigation validation
class StraightLineStepper {
 public:
  using Jacobian = BoundMatrix;
  using Covariance = BoundSquareMatrix;
  using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;
  using CurvilinearState =
      std::tuple<CurvilinearTrackParameters, Jacobian, double>;
  using BField = NullBField;

  struct Config {};

  struct Options : public StepperPlainOptions {
    Options(const GeometryContext& gctx, const MagneticFieldContext& mctx)
        : StepperPlainOptions(gctx, mctx) {}

    void setPlainOptions(const StepperPlainOptions& options) {
      static_cast<StepperPlainOptions&>(*this) = options;
    }
  };

  /// State for track parameter propagation
  ///
  struct State {
    /// Constructor from the initial bound track parameters
    ///
    /// @param [in] optionsIn The options for the stepper
    ///
    /// @note the covariance matrix is copied when needed
    explicit State(const Options& optionsIn) : options(optionsIn) {}

    Options options;

    /// Jacobian from local to the global frame
    BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();

    /// Pure transport jacobian part from runge kutta integration
    FreeMatrix jacTransport = FreeMatrix::Identity();

    /// The full jacobian of the transport entire transport
    Jacobian jacobian = Jacobian::Identity();

    /// The propagation derivative
    FreeVector derivative = FreeVector::Zero();

    /// Internal free vector parameters
    FreeVector pars = FreeVector::Zero();

    /// Particle hypothesis
    ParticleHypothesis particleHypothesis = ParticleHypothesis::pion();

    /// Boolean to indicate if you need covariance transport
    bool covTransport = false;
    Covariance cov = Covariance::Zero();

    /// accummulated path length state
    double pathAccumulated = 0.;

    /// Total number of performed steps
    std::size_t nSteps = 0;

    /// Totoal number of attempted steps
    std::size_t nStepTrials = 0;

    /// adaptive step size of the runge-kutta integration
    ConstrainedStep stepSize;

    // Previous step size for overstep estimation (ignored for SL stepper)
    double previousStepSize = 0.;

    /// Statistics of the stepper
    StepperStatistics statistics;
  };

  StraightLineStepper() = default;

  State makeState(const Options& options,
                  const BoundTrackParameters& par) const {
    State state{options};

    state.particleHypothesis = par.particleHypothesis();

    Vector3 position = par.position(options.geoContext);
    Vector3 direction = par.direction();
    state.pars.template segment<3>(eFreePos0) = position;
    state.pars.template segment<3>(eFreeDir0) = direction;
    state.pars[eFreeTime] = par.time();
    state.pars[eFreeQOverP] = par.parameters()[eBoundQOverP];

    // Init the jacobian matrix if needed
    if (par.covariance()) {
      // Get the reference surface for navigation
      const auto& surface = par.referenceSurface();
      // set the covariance transport flag to true and copy
      state.covTransport = true;
      state.cov = BoundSquareMatrix(*par.covariance());
      state.jacToGlobal =
          surface.boundToFreeJacobian(options.geoContext, position, direction);
    }

    state.stepSize = ConstrainedStep(options.maxStepSize);

    return state;
  }

  /// @brief Resets the state
  ///
  /// @param [in, out] state State of the stepper
  /// @param [in] boundParams Parameters in bound parametrisation
  /// @param [in] cov Covariance matrix
  /// @param [in] surface The reset @c State will be on this surface
  /// @param [in] stepSize Step size
  void resetState(
      State& state, const BoundVector& boundParams,
      const BoundSquareMatrix& cov, const Surface& surface,
      const double stepSize = std::numeric_limits<double>::max()) const;

  /// Get the field for the stepping, this gives back a zero field
  Result<Vector3> getField(State& /*state*/, const Vector3& /*pos*/) const {
    // get the field from the cell
    return Result<Vector3>::success({0., 0., 0.});
  }

  /// Global particle position accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  Vector3 position(const State& state) const {
    return state.pars.template segment<3>(eFreePos0);
  }

  /// Momentum direction accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  Vector3 direction(const State& state) const {
    return state.pars.template segment<3>(eFreeDir0);
  }

  /// QoP direction accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double qOverP(const State& state) const { return state.pars[eFreeQOverP]; }

  /// Absolute momentum accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double absoluteMomentum(const State& state) const {
    return particleHypothesis(state).extractMomentum(qOverP(state));
  }

  /// Momentum accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  Vector3 momentum(const State& state) const {
    return absoluteMomentum(state) * direction(state);
  }

  /// Charge access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double charge(const State& state) const {
    return particleHypothesis(state).extractCharge(qOverP(state));
  }

  /// Particle hypothesis
  ///
  /// @param state [in] The stepping state (thread-local cache)
  const ParticleHypothesis& particleHypothesis(const State& state) const {
    return state.particleHypothesis;
  }

  /// Time access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double time(const State& state) const { return state.pars[eFreeTime]; }

  /// Update surface status
  ///
  /// This method intersects the provided surface and update the navigation
  /// step estimation accordingly (hence it changes the state). It also
  /// returns the status of the intersection to trigger onSurface in case
  /// the surface is reached.
  ///
  /// @param [in,out] state The stepping state (thread-local cache)
  /// @param [in] surface The surface provided
  /// @param [in] index The surface intersection index
  /// @param [in] navDir The navigation direction
  /// @param [in] boundaryTolerance The boundary check for this status update
  /// @param [in] surfaceTolerance Surface tolerance used for intersection
  /// @param [in] stype The step size type to be set
  /// @param [in] logger A logger instance
  IntersectionStatus updateSurfaceStatus(
      State& state, const Surface& surface, std::uint8_t index,
      Direction navDir, const BoundaryTolerance& boundaryTolerance,
      double surfaceTolerance, ConstrainedStep::Type stype,
      const Logger& logger = getDummyLogger()) const {
    return detail::updateSingleSurfaceStatus<StraightLineStepper>(
        *this, state, surface, index, navDir, boundaryTolerance,
        surfaceTolerance, stype, logger);
  }

  /// Update step size
  ///
  /// It checks the status to the reference surface & updates
  /// the step size accordingly
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param oIntersection [in] The ObjectIntersection to layer, boundary, etc
  /// @param direction [in] The propagation direction
  /// @param stype [in] The step size type to be set
  template <typename object_intersection_t>
  void updateStepSize(State& state, const object_intersection_t& oIntersection,
                      Direction direction, ConstrainedStep::Type stype) const {
    (void)direction;
    double stepSize = oIntersection.pathLength();
    updateStepSize(state, stepSize, stype);
  }

  /// Update step size - explicitly with a double
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param stepSize [in] The step size value
  /// @param stype [in] The step size type to be set
  void updateStepSize(State& state, double stepSize,
                      ConstrainedStep::Type stype) const {
    state.previousStepSize = state.stepSize.value();
    state.stepSize.update(stepSize, stype);
  }

  /// Get the step size
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @param stype [in] The step size type to be returned
  double getStepSize(const State& state, ConstrainedStep::Type stype) const {
    return state.stepSize.value(stype);
  }

  /// Release the Step size
  ///
  /// @param [in,out] state The stepping state (thread-local cache)
  /// @param [in] stype The step size type to be released
  void releaseStepSize(State& state, ConstrainedStep::Type stype) const {
    state.stepSize.release(stype);
  }

  /// Output the Step Size - single component
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  std::string outputStepSize(const State& state) const {
    return state.stepSize.toString();
  }

  /// Create and return the bound state at the current position
  ///
  /// @brief It does not check if the transported state is at the surface, this
  /// needs to be guaranteed by the propagator
  ///
  /// @param [in] state State that will be presented as @c BoundState
  /// @param [in] surface The surface to which we bind the state
  /// @param [in] transportCov Flag steering covariance transport
  /// @param [in] freeToBoundCorrection Correction for non-linearity effect during transform from free to bound
  ///
  /// @return A bound state:
  ///   - the parameters at the surface
  ///   - the stepwise jacobian towards it (from last bound)
  ///   - and the path length (from start - for ordering)
  Result<BoundState> boundState(
      State& state, const Surface& surface, bool transportCov = true,
      const FreeToBoundCorrection& freeToBoundCorrection =
          FreeToBoundCorrection(false)) const;

  /// @brief If necessary fill additional members needed for curvilinearState
  ///
  /// Compute path length derivatives in case they have not been computed
  /// yet, which is the case if no step has been executed yet.
  ///
  /// @param [in, out] state The stepping state (thread-local cache)
  /// @return true if nothing is missing after this call, false otherwise.
  bool prepareCurvilinearState(State& state) const {
    // test whether the accumulated path has still its initial value.
    if (state.pathAccumulated != 0) {
      return true;
    }

    // dr/ds :
    state.derivative.template head<3>() = direction(state);
    // dt / ds
    state.derivative(eFreeTime) = fastHypot(
        1., state.particleHypothesis.mass() / absoluteMomentum(state));
    // d (dr/ds) / ds : == 0
    state.derivative.template segment<3>(4) = Acts::Vector3::Zero().transpose();
    // d qop / ds  == 0
    state.derivative(eFreeQOverP) = 0.;

    return true;
  }

  /// Create and return a curvilinear state at the current position
  ///
  /// @brief This creates a curvilinear state.
  ///
  /// @param [in] state State that will be presented as @c CurvilinearState
  /// @param [in] transportCov Flag steering covariance transport
  ///
  /// @return A curvilinear state:
  ///   - the curvilinear parameters at given position
  ///   - the stepweise jacobian towards it (from last bound)
  ///   - and the path length (from start - for ordering)
  CurvilinearState curvilinearState(State& state,
                                    bool transportCov = true) const;

  /// Method to update a stepper state to the some parameters
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] freeParams Free parameters that will be written into @p state
  /// @param [in] boundParams Corresponding bound parameters used to update jacToGlobal in @p state
  /// @param [in] covariance Covariance that will be written into @p state
  /// @param [in] surface The surface used to update the jacToGlobal
  void update(State& state, const FreeVector& freeParams,
              const BoundVector& boundParams, const Covariance& covariance,
              const Surface& surface) const;

  /// Method to update the stepper state
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] uposition the updated position
  /// @param [in] udirection the updated direction
  /// @param [in] qop the updated qop value
  /// @param [in] time the updated time value
  void update(State& state, const Vector3& uposition, const Vector3& udirection,
              double qop, double time) const;

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state - for the moment a dummy method
  ///
  /// @param [in,out] state State of the stepper
  void transportCovarianceToCurvilinear(State& state) const;

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state - for the moment a dummy method
  ///
  /// @tparam surface_t the surface type - ignored here
  ///
  /// @param [in,out] state The stepper state
  /// @param [in] surface is the surface to which the covariance is
  ///        forwarded to
  /// @note no check is done if the position is actually on the surface
  /// @param [in] freeToBoundCorrection Correction for non-linearity effect during transform from free to bound
  ///
  void transportCovarianceToBound(
      State& state, const Surface& surface,
      const FreeToBoundCorrection& freeToBoundCorrection =
          FreeToBoundCorrection(false)) const;

  /// Perform a straight line propagation step
  ///
  /// @param [in,out] state The state of the stepper
  /// @param propDir is the direction of propagation
  /// @param material is the material properties
  /// @return the result of the step
  ///
  /// The state contains the desired step size,
  /// it can be negative during backwards track propagation,
  /// and since we're using an adaptive algorithm, it can
  /// be modified by the stepper class during propagation.
  Result<double> step(State& state, Direction propDir,
                      const IVolumeMaterial* material) const {
    (void)material;

    // use the adjusted step size
    const auto h = state.stepSize.value() * propDir;
    const auto m = state.particleHypothesis.mass();
    const auto p = absoluteMomentum(state);
    // time propagates along distance as 1/b = sqrt(1 + m²/p²)
    const auto dtds = fastHypot(1., m / p);
    // Update the track parameters according to the equations of motion
    Vector3 dir = direction(state);
    state.pars.template segment<3>(eFreePos0) += h * dir;
    state.pars[eFreeTime] += h * dtds;

    // Propagate the jacobian
    if (state.covTransport) {
      // The step transport matrix in global coordinates
      FreeMatrix D = FreeMatrix::Identity();
      D.block<3, 3>(0, 4) = ActsSquareMatrix<3>::Identity() * h;
      // Extend the calculation by the time propagation
      // Evaluate dt/dlambda
      D(3, 7) = h * m * m * state.pars[eFreeQOverP] / dtds;
      // Set the derivative factor the time
      state.derivative(3) = dtds;
      // Update jacobian and derivative
      state.jacTransport = D * state.jacTransport;
      state.derivative.template head<3>() = dir;
    }

    // state the path length
    state.pathAccumulated += h;
    ++state.nSteps;
    ++state.nStepTrials;

    ++state.statistics.nAttemptedSteps;
    ++state.statistics.nSuccessfulSteps;
    if (propDir != Direction::fromScalarZeroAsPositive(h)) {
      ++state.statistics.nReverseSteps;
    }
    state.statistics.pathLength += h;
    state.statistics.absolutePathLength += std::abs(h);

    return h;
  }
};

template <>
struct SupportsBoundParameters<StraightLineStepper> : public std::true_type {};

}  // namespace Acts
