// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/PropagatorTraits.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

namespace Acts {

class OptimizedStepper {
 public:
  /// Jacobian, Covariance and State definitions
  using Jacobian = BoundMatrix;
  using Covariance = BoundSquareMatrix;
  using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;
  using CurvilinearState =
      std::tuple<CurvilinearTrackParameters, Jacobian, double>;

  /// @brief State for track parameter propagation
  ///
  /// It contains the stepping information and is provided thread local
  /// by the propagator
  struct State {
    State() = delete;

    /// Constructor from the initial bound track parameters
    ///
    /// @param [in] gctx is the context object for the geometry
    /// @param [in] fieldCacheIn is the cache object for the magnetic field
    /// @param [in] par The track parameters at start
    /// @param [in] ssize is the maximum step size
    ///
    /// @note the covariance matrix is copied when needed
    explicit State(const GeometryContext& gctx,
                   MagneticFieldProvider::Cache fieldCacheIn,
                   const BoundTrackParameters& par,
                   double ssize = std::numeric_limits<double>::max())
        : particleHypothesis(par.particleHypothesis()),
          stepSize(ssize),
          fieldCache(std::move(fieldCacheIn)),
          geoContext(gctx) {
      Vector3 position = par.position(gctx);
      Vector3 direction = par.direction();
      pars.template segment<3>(eFreePos0) = position;
      pars.template segment<3>(eFreeDir0) = direction;
      pars[eFreeTime] = par.time();
      pars[eFreeQOverP] = par.parameters()[eBoundQOverP];

      // Init the jacobian matrix if needed
      if (par.covariance()) {
        // Get the reference surface for navigation
        const auto& surface = par.referenceSurface();
        // set the covariance transport flag to true and copy
        covTransport = true;
        cov = BoundSquareMatrix(*par.covariance());
        jacToGlobal = surface.boundToFreeJacobian(gctx, position, direction);
      }
    }

    /// Internal free vector parameters
    FreeVector pars = FreeVector::Zero();

    /// Particle hypothesis
    ParticleHypothesis particleHypothesis = ParticleHypothesis::pion();

    /// Covariance matrix (and indicator)
    /// associated with the initial error on track parameters
    bool covTransport = false;
    Covariance cov = Covariance::Zero();

    /// The full jacobian of the transport entire transport
    Jacobian jacobian = Jacobian::Identity();

    /// Jacobian from local to the global frame
    BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();

    /// Pure transport jacobian part from runge kutta integration
    FreeMatrix jacTransport = FreeMatrix::Identity();

    /// The propagation derivative
    FreeVector derivative = FreeVector::Zero();

    /// Accummulated path length state
    double pathAccumulated = 0.;

    /// Adaptive step size of the runge-kutta integration
    ConstrainedStep stepSize;

    /// Last performed step (for overstep limit calculation)
    double previousStepSize = 0.;

    /// This caches the current magnetic field cell and stays
    /// (and interpolates) within it as long as this is valid.
    /// See step() code for details.
    MagneticFieldProvider::Cache fieldCache;

    /// The geometry context
    std::reference_wrapper<const GeometryContext> geoContext;

    /// @brief Storage of magnetic field and the sub steps during a RKN4 step
    struct {
      /// Magnetic field evaulations
      Vector3 B_first, B_middle, B_last;
      /// k_i of the RKN4 algorithm
      Vector3 k1, k2, k3, k4;
      /// k_i elements of the momenta
      std::array<double, 4> kQoP{};
    } stepData;
  };

  /// Constructor requires knowledge of the detector's magnetic field
  OptimizedStepper(std::shared_ptr<const MagneticFieldProvider> bField,
                   double overstepLimit = 100 * UnitConstants::um);

  State makeState(std::reference_wrapper<const GeometryContext> gctx,
                  std::reference_wrapper<const MagneticFieldContext> mctx,
                  const BoundTrackParameters& par,
                  double ssize = std::numeric_limits<double>::max()) const;

  /// @brief Resets the state
  ///
  /// @param [in, out] state State of the stepper
  /// @param [in] boundParams Parameters in bound parametrisation
  /// @param [in] cov Covariance matrix
  /// @param [in] surface The reference surface of the bound parameters
  /// @param [in] stepSize Step size
  void resetState(
      State& state, const BoundVector& boundParams,
      const BoundSquareMatrix& cov, const Surface& surface,
      const double stepSize = std::numeric_limits<double>::max()) const;

  /// Get the field for the stepping, it checks first if the access is still
  /// within the Cell, and updates the cell if necessary.
  ///
  /// @param [in,out] state is the propagation state associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  Result<Vector3> getField(State& state, const Vector3& pos) const {
    // get the field from the cell
    return m_bField->getField(pos, state.fieldCache);
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
  /// It checks the status to the reference surface & updates
  /// the step size accordingly
  ///
  /// @param [in,out] state The stepping state (thread-local cache)
  /// @param [in] surface The surface provided
  /// @param [in] index The surface intersection index
  /// @param [in] navDir The navigation direction
  /// @param [in] bcheck The boundary check for this status update
  /// @param [in] surfaceTolerance Surface tolerance used for intersection
  /// @param [in] logger A @c Logger instance
  Intersection3D::Status updateSurfaceStatus(
      State& state, const Surface& surface, std::uint8_t index,
      Direction navDir, const BoundaryCheck& bcheck,
      ActsScalar surfaceTolerance = s_onSurfaceTolerance,
      const Logger& logger = getDummyLogger()) const {
    return detail::updateSingleSurfaceStatus<OptimizedStepper>(
        *this, state, surface, index, navDir, bcheck, surfaceTolerance, logger);
  }

  /// Update step size
  ///
  /// This method intersects the provided surface and update the navigation
  /// step estimation accordingly (hence it changes the state). It also
  /// returns the status of the intersection to trigger onSurface in case
  /// the surface is reached.
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param oIntersection [in] The ObjectIntersection to layer, boundary, etc
  /// @param release [in] boolean to trigger step size release
  template <typename object_intersection_t>
  void updateStepSize(State& state, const object_intersection_t& oIntersection,
                      Direction /*direction*/, bool release = true) const {
    detail::updateSingleStepSize<OptimizedStepper>(state, oIntersection,
                                                   release);
  }

  /// Update step size - explicitly with a double
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param stepSize [in] The step size value
  /// @param stype [in] The step size type to be set
  /// @param release [in] Do we release the step size?
  void updateStepSize(State& state, double stepSize,
                      ConstrainedStep::Type stype, bool release = true) const {
    state.previousStepSize = state.stepSize.value();
    state.stepSize.update(stepSize, stype, release);
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
  /// @param state [in,out] The stepping state (thread-local cache)
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

  /// Overstep limit
  double overstepLimit(const State& /*state*/) const {
    // A dynamic overstep limit could sit here
    return -m_overstepLimit;
  }

  /// Create and return the bound state at the current position
  ///
  /// @brief This transports (if necessary) the covariance
  /// to the surface and creates a bound state. It does not check
  /// if the transported state is at the surface, this needs to
  /// be guaranteed by the propagator
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
  /// @param [in, out] prop_state State that will be presented as @c BoundState
  /// @param [in] navigator the navigator of the propagation
  /// @return true if nothing is missing after this call, false otherwise.
  template <typename propagator_state_t, typename navigator_t>
  bool prepareCurvilinearState(propagator_state_t& prop_state,
                               const navigator_t& navigator) const;

  /// Create and return a curvilinear state at the current position
  ///
  /// @brief This transports (if necessary) the covariance
  /// to the current position and creates a curvilinear state.
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
  /// @param [in] covariance The covariance that will be written into @p state
  /// @param [in] surface The surface used to update the jacToGlobal
  void update(State& state, const FreeVector& freeParams,
              const BoundVector& boundParams, const Covariance& covariance,
              const Surface& surface) const;

  /// Method to update the stepper state
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] uposition the updated position
  /// @param [in] udirection the updated direction
  /// @param [in] qOverP the updated qOverP value
  /// @param [in] time the updated time value
  void update(State& state, const Vector3& uposition, const Vector3& udirection,
              double qOverP, double time) const;

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state
  ///
  /// @param [in,out] state State of the stepper
  void transportCovarianceToCurvilinear(State& state) const;

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current position,
  /// or direction of the state
  ///
  /// @tparam surface_t the Surface type
  ///
  /// @param [in,out] state State of the stepper
  /// @param [in] surface is the surface to which the covariance is forwarded to
  /// @param [in] freeToBoundCorrection Correction for non-linearity effect during transform from free to bound
  /// @note no check is done if the position is actually on the surface
  void transportCovarianceToBound(
      State& state, const Surface& surface,
      const FreeToBoundCorrection& freeToBoundCorrection =
          FreeToBoundCorrection(false)) const;

  /// Perform a Runge-Kutta track parameter propagation step
  ///
  /// @param [in,out] state the propagation state
  /// @param [in] navigator the navigator of the propagation
  /// @note The state contains the desired step size.  It can be negative during
  ///       backwards track propagation, and since we're using an adaptive
  ///       algorithm, it can be modified by the stepper class during
  ///       propagation.
  template <typename propagator_state_t, typename navigator_t>
  Result<double> step(propagator_state_t& state,
                      const navigator_t& navigator) const;

  /// Method that reset the Jacobian to the Identity for when no bound state are
  /// available
  ///
  /// @param [in,out] state State of the stepper
  void setIdentityJacobian(State& state) const;

 protected:
  /// Magnetic field inside of the detector
  std::shared_ptr<const MagneticFieldProvider> m_bField;

  /// Overstep limit
  double m_overstepLimit;
};

template <typename navigator_t>
struct SupportsBoundParameters<OptimizedStepper, navigator_t>
    : public std::true_type {};

}  // namespace Acts
