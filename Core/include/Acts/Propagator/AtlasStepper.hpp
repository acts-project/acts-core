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
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/StepperOptions.hpp"
#include "Acts/Propagator/StepperStatistics.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cmath>

namespace Acts {

class IVolumeMaterial;

/// @brief the AtlasStepper implementation for the
///
/// This is based original stepper code from the ATLAS RungeKuttaPropagator
class AtlasStepper {
 public:
  using Jacobian = BoundMatrix;
  using Covariance = BoundSquareMatrix;
  using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;
  using CurvilinearState =
      std::tuple<CurvilinearTrackParameters, Jacobian, double>;

  struct Config {
    std::shared_ptr<const MagneticFieldProvider> bField;
  };

  struct Options : public StepperPlainOptions {
    Options(const GeometryContext& gctx, const MagneticFieldContext& mctx)
        : StepperPlainOptions(gctx, mctx) {}

    void setPlainOptions(const StepperPlainOptions& options) {
      static_cast<StepperPlainOptions&>(*this) = options;
    }
  };

  /// @brief Nested State struct for the local caching
  struct State {
    /// Constructor
    ///
    /// @tparam Type of TrackParameters
    ///
    /// @param [in] optionsIn The options for the stepper
    /// @param [in] fieldCacheIn The magnetic field cache
    State(const Options& optionsIn, MagneticFieldProvider::Cache fieldCacheIn)
        : options(optionsIn), fieldCache(std::move(fieldCacheIn)) {}

    Options options;

    ParticleHypothesis particleHypothesis = ParticleHypothesis::pion();

    // optimisation that init is not called twice
    bool state_ready = false;
    // configuration
    bool useJacobian = false;
    double step = 0;
    double maxPathLength = 0;
    bool mcondition = false;
    bool needgradient = false;
    bool newfield = true;
    // internal parameters to be used
    Vector3 field = Vector3::Zero();
    std::array<double, 60> pVector{};

    /// Storage pattern of pVector
    ///                   /dL0    /dL1    /dPhi   /the   /dCM   /dT
    /// X  ->P[0]  dX /   P[ 8]   P[16]   P[24]   P[32]   P[40]  P[48]
    /// Y  ->P[1]  dY /   P[ 9]   P[17]   P[25]   P[33]   P[41]  P[49]
    /// Z  ->P[2]  dZ /   P[10]   P[18]   P[26]   P[34]   P[42]  P[50]
    /// T  ->P[3]  dT/    P[11]   P[19]   P[27]   P[35]   P[43]  P[51]
    /// Ax ->P[4]  dAx/   P[12]   P[20]   P[28]   P[36]   P[44]  P[52]
    /// Ay ->P[5]  dAy/   P[13]   P[21]   P[29]   P[37]   P[45]  P[53]
    /// Az ->P[6]  dAz/   P[14]   P[22]   P[30]   P[38]   P[46]  P[54]
    /// CM ->P[7]  dCM/   P[15]   P[23]   P[31]   P[39]   P[47]  P[55]
    /// Cache: P[56] - P[59]

    // result
    double parameters[eBoundSize] = {0., 0., 0., 0., 0., 0.};
    const Covariance* covariance = nullptr;
    Covariance cov = Covariance::Zero();
    bool covTransport = false;
    double jacobian[eBoundSize * eBoundSize] = {};

    // accummulated path length cache
    double pathAccumulated = 0.;

    /// Total number of performed steps
    std::size_t nSteps = 0;

    /// Totoal number of attempted steps
    std::size_t nStepTrials = 0;

    // Adaptive step size of the runge-kutta integration
    ConstrainedStep stepSize;

    // Previous step size for overstep estimation
    double previousStepSize = 0.;

    /// It caches the current magnetic field cell and stays (and interpolates)
    ///  within as long as this is valid. See step() code for details.
    MagneticFieldProvider::Cache fieldCache;

    /// Debug output
    /// the string where debug messages are stored (optionally)
    bool debug = false;
    std::string debugString = "";
    /// buffer & formatting for consistent output
    std::size_t debugPfxWidth = 30;
    std::size_t debugMsgWidth = 50;

    /// The statistics of the stepper
    StepperStatistics statistics;
  };

  explicit AtlasStepper(std::shared_ptr<const MagneticFieldProvider> bField)
      : m_bField(std::move(bField)) {}

  explicit AtlasStepper(const Config& config) : m_bField(config.bField) {}

  State makeState(const Options& options,
                  const BoundTrackParameters& par) const {
    State state{options, m_bField->makeCache(options.magFieldContext)};

    state.particleHypothesis = par.particleHypothesis();

    // The rest of this constructor is copy&paste of AtlasStepper::update() -
    // this is a nasty but working solution for the stepper state without
    // functions

    const auto pos = par.position(options.geoContext);
    const auto Vp = par.parameters();

    double Sf = std::sin(Vp[eBoundPhi]);
    double Cf = std::cos(Vp[eBoundPhi]);
    double Se = std::sin(Vp[eBoundTheta]);
    double Ce = std::cos(Vp[eBoundTheta]);

    double* pVector = state.pVector.data();

    pVector[0] = pos[ePos0];
    pVector[1] = pos[ePos1];
    pVector[2] = pos[ePos2];
    pVector[3] = par.time();
    pVector[4] = Cf * Se;
    pVector[5] = Sf * Se;
    pVector[6] = Ce;
    pVector[7] = Vp[eBoundQOverP];

    // @todo: remove magic numbers - is that the charge ?
    if (std::abs(pVector[7]) < .000000000000001) {
      pVector[7] < 0. ? pVector[7] = -.000000000000001
                      : pVector[7] = .000000000000001;
    }

    // prepare the jacobian if we have a covariance
    if (par.covariance()) {
      // copy the covariance matrix
      state.covariance = new BoundSquareMatrix(*par.covariance());
      state.covTransport = true;
      state.useJacobian = true;
      const auto transform = par.referenceSurface().referenceFrame(
          options.geoContext, pos, par.direction());

      pVector[8] = transform(0, eBoundLoc0);
      pVector[16] = transform(0, eBoundLoc1);
      pVector[24] = 0.;
      pVector[32] = 0.;
      pVector[40] = 0.;
      pVector[48] = 0.;  // dX /

      pVector[9] = transform(1, eBoundLoc0);
      pVector[17] = transform(1, eBoundLoc1);
      pVector[25] = 0.;
      pVector[33] = 0.;
      pVector[41] = 0.;
      pVector[49] = 0.;  // dY /

      pVector[10] = transform(2, eBoundLoc0);
      pVector[18] = transform(2, eBoundLoc1);
      pVector[26] = 0.;
      pVector[34] = 0.;
      pVector[42] = 0.;
      pVector[50] = 0.;  // dZ /

      pVector[11] = 0.;
      pVector[19] = 0.;
      pVector[27] = 0.;
      pVector[35] = 0.;
      pVector[43] = 0.;
      pVector[51] = 1.;  // dT/

      pVector[12] = 0.;
      pVector[20] = 0.;
      pVector[28] = -Sf * Se;  // - sin(phi) * cos(theta)
      pVector[36] = Cf * Ce;   // cos(phi) * cos(theta)
      pVector[44] = 0.;
      pVector[52] = 0.;  // dAx/

      pVector[13] = 0.;
      pVector[21] = 0.;
      pVector[29] = Cf * Se;  // cos(phi) * sin(theta)
      pVector[37] = Sf * Ce;  // sin(phi) * cos(theta)
      pVector[45] = 0.;
      pVector[53] = 0.;  // dAy/

      pVector[14] = 0.;
      pVector[22] = 0.;
      pVector[30] = 0.;
      pVector[38] = -Se;  // - sin(theta)
      pVector[46] = 0.;
      pVector[54] = 0.;  // dAz/

      pVector[15] = 0.;
      pVector[23] = 0.;
      pVector[31] = 0.;
      pVector[39] = 0.;
      pVector[47] = 1.;
      pVector[55] = 0.;  // dCM/

      pVector[56] = 0.;
      pVector[57] = 0.;
      pVector[58] = 0.;
      pVector[59] = 0.;

      // special treatment for surface types
      const auto& surface = par.referenceSurface();
      // the disc needs polar coordinate adaptations
      if (surface.type() == Surface::Disc) {
        double lCf = std::cos(Vp[1]);
        double lSf = std::sin(Vp[1]);
        double Ax[3] = {transform(0, 0), transform(1, 0), transform(2, 0)};
        double Ay[3] = {transform(0, 1), transform(1, 1), transform(2, 1)};
        double d0 = lCf * Ax[0] + lSf * Ay[0];
        double d1 = lCf * Ax[1] + lSf * Ay[1];
        double d2 = lCf * Ax[2] + lSf * Ay[2];
        pVector[8] = d0;
        pVector[9] = d1;
        pVector[10] = d2;
        pVector[16] = Vp[0] * (lCf * Ay[0] - lSf * Ax[0]);
        pVector[17] = Vp[0] * (lCf * Ay[1] - lSf * Ax[1]);
        pVector[18] = Vp[0] * (lCf * Ay[2] - lSf * Ax[2]);
      }
      // the line needs components that relate direction change
      // with global frame change
      if (surface.type() == Surface::Perigee ||
          surface.type() == Surface::Straw) {
        // sticking to the nomenclature of the original RkPropagator
        // - axis pointing along the drift/transverse direction
        double B[3] = {transform(0, 0), transform(1, 0), transform(2, 0)};
        // - axis along the straw
        double A[3] = {transform(0, 1), transform(1, 1), transform(2, 1)};
        // - normal vector of the reference frame
        double C[3] = {transform(0, 2), transform(1, 2), transform(2, 2)};

        // projection of direction onto normal vector of reference frame
        double PC = pVector[4] * C[0] + pVector[5] * C[1] + pVector[6] * C[2];
        double Bn = 1. / PC;

        double Bx2 = -A[2] * pVector[29];
        double Bx3 = A[1] * pVector[38] - A[2] * pVector[37];

        double By2 = A[2] * pVector[28];
        double By3 = A[2] * pVector[36] - A[0] * pVector[38];

        double Bz2 = A[0] * pVector[29] - A[1] * pVector[28];
        double Bz3 = A[0] * pVector[37] - A[1] * pVector[36];

        double B2 = B[0] * Bx2 + B[1] * By2 + B[2] * Bz2;
        double B3 = B[0] * Bx3 + B[1] * By3 + B[2] * Bz3;

        Bx2 = (Bx2 - B[0] * B2) * Bn;
        Bx3 = (Bx3 - B[0] * B3) * Bn;
        By2 = (By2 - B[1] * B2) * Bn;
        By3 = (By3 - B[1] * B3) * Bn;
        Bz2 = (Bz2 - B[2] * B2) * Bn;
        Bz3 = (Bz3 - B[2] * B3) * Bn;

        //  /dPhi      |     /dThe       |
        pVector[24] = Bx2 * Vp[0];
        pVector[32] = Bx3 * Vp[0];  // dX/
        pVector[25] = By2 * Vp[0];
        pVector[33] = By3 * Vp[0];  // dY/
        pVector[26] = Bz2 * Vp[0];
        pVector[34] = Bz3 * Vp[0];  // dZ/
      }
    }

    state.stepSize = ConstrainedStep(options.maxStepSize);

    // now declare the state as ready
    state.state_ready = true;

    return state;
  }

  /// @brief Resets the state
  ///
  /// @param [in, out] state State of the stepper
  /// @param [in] boundParams Parameters in bound parametrisation
  /// @param [in] cov Covariance matrix
  /// @param [in] surface Reset state will be on this surface
  /// @param [in] stepSize Step size
  void resetState(
      State& state, const BoundVector& boundParams,
      const BoundSquareMatrix& cov, const Surface& surface,
      const double stepSize = std::numeric_limits<double>::max()) const {
    // Update the stepping state
    update(state,
           transformBoundToFreeParameters(surface, state.options.geoContext,
                                          boundParams),
           boundParams, cov, surface);
    state.stepSize = ConstrainedStep(stepSize);
    state.pathAccumulated = 0.;

    setIdentityJacobian(state);
  }

  /// Get the field for the stepping
  /// It checks first if the access is still within the Cell,
  /// and updates the cell if necessary, then it takes the field
  /// from the cell
  /// @param [in,out] state is the stepper state associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  Result<Vector3> getField(State& state, const Vector3& pos) const {
    // get the field from the cell
    auto res = m_bField->getField(pos, state.fieldCache);
    if (res.ok()) {
      state.field = *res;
    }
    return res;
  }

  Vector3 position(const State& state) const {
    return Vector3(state.pVector[0], state.pVector[1], state.pVector[2]);
  }

  Vector3 direction(const State& state) const {
    return Vector3(state.pVector[4], state.pVector[5], state.pVector[6]);
  }

  double qOverP(const State& state) const { return state.pVector[7]; }

  /// Absolute momentum accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double absoluteMomentum(const State& state) const {
    return particleHypothesis(state).extractMomentum(qOverP(state));
  }

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

  /// Overstep limit
  double overstepLimit(const State& /*state*/) const {
    return -m_overstepLimit;
  }

  /// Time access
  double time(const State& state) const { return state.pVector[3]; }

  /// Update surface status
  ///
  /// This method intersect the provided surface and update the navigation
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
  /// @param [in] logger Logger instance to use
  IntersectionStatus updateSurfaceStatus(
      State& state, const Surface& surface, std::uint8_t index,
      Direction navDir, const BoundaryTolerance& boundaryTolerance,
      double surfaceTolerance, ConstrainedStep::Type stype,
      const Logger& logger = getDummyLogger()) const {
    return detail::updateSingleSurfaceStatus<AtlasStepper>(
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
  /// @param [in,out] state The stepping state (thread-local cache)
  /// @param [in] stepSize The step size value
  /// @param [in] stype The step size type to be set
  void updateStepSize(State& state, double stepSize,
                      ConstrainedStep::Type stype) const {
    state.previousStepSize = state.stepSize.value();
    state.stepSize.update(stepSize, stype);
  }

  /// Release the Step size
  ///
  /// @param [in,out] state The stepping state (thread-local cache)
  /// @param [in] stype The step size type to be released
  void releaseStepSize(State& state, ConstrainedStep::Type stype) const {
    state.stepSize.release(stype);
  }

  /// Get the step size
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @param stype [in] The step size type to be returned
  double getStepSize(const State& state, ConstrainedStep::Type stype) const {
    return state.stepSize.value(stype);
  }

  /// Output the Step Size - single component
  ///
  /// @param [in,out] state The stepping state (thread-local cache)
  std::string outputStepSize(const State& state) const {
    return state.stepSize.toString();
  }

  /// Create and return the bound state at the current position
  ///
  ///
  /// @param [in] state State that will be presented as @c BoundState
  /// @param [in] surface The surface to which we bind the state
  /// @param [in] transportCov Flag steering covariance transport
  /// @param [in] freeToBoundCorrection Correction for non-linearity effect during transform from free to bound
  ///
  /// @return A bound state:
  ///   - the parameters at the surface
  ///   - the stepwise jacobian towards it
  ///   - and the path length (from start - for ordering)
  Result<BoundState> boundState(
      State& state, const Surface& surface, bool transportCov = true,
      const FreeToBoundCorrection& freeToBoundCorrection =
          FreeToBoundCorrection(false)) const {
    // the convert method invalidates the state (in case it's reused)
    state.state_ready = false;
    // extract state information
    Acts::Vector4 pos4;
    pos4[ePos0] = state.pVector[0];
    pos4[ePos1] = state.pVector[1];
    pos4[ePos2] = state.pVector[2];
    pos4[eTime] = state.pVector[3];
    Acts::Vector3 dir;
    dir[eMom0] = state.pVector[4];
    dir[eMom1] = state.pVector[5];
    dir[eMom2] = state.pVector[6];
    const auto qOverP = state.pVector[7];

    // The transport of the covariance
    std::optional<Covariance> covOpt = std::nullopt;
    if (state.covTransport && transportCov) {
      transportCovarianceToBound(state, surface, freeToBoundCorrection);
    }
    if (state.cov != Covariance::Zero()) {
      covOpt = state.cov;
    }

    // Fill the end parameters
    auto parameters = BoundTrackParameters::create(
        surface.getSharedPtr(), state.options.geoContext, pos4, dir, qOverP,
        std::move(covOpt), state.particleHypothesis);
    if (!parameters.ok()) {
      return parameters.error();
    }

    Jacobian jacobian(state.jacobian);

    return BoundState(std::move(*parameters), jacobian.transpose(),
                      state.pathAccumulated);
  }

  /// @brief If necessary fill additional members needed for curvilinearState
  ///
  /// Compute path length derivatives in case they have not been computed
  /// yet, which is the case if no step has been executed yet.
  ///
  /// @param [in, out] state The stepping state (thread-local cache)
  /// @return true if nothing is missing after this call, false otherwise.
  bool prepareCurvilinearState(State& state) const {
    (void)state;
    return true;
  }

  /// Create and return a curvilinear state at the current position
  ///
  ///
  /// @param [in] state State that will be presented as @c CurvilinearState
  /// @param [in] transportCov Flag steering covariance transport
  ///
  /// @return A curvilinear state:
  ///   - the curvilinear parameters at given position
  ///   - the stepweise jacobian towards it
  ///   - and the path length (from start - for ordering)
  CurvilinearState curvilinearState(State& state,
                                    bool transportCov = true) const {
    // the convert method invalidates the state (in case it's reused)
    state.state_ready = false;
    // extract state information
    Acts::Vector4 pos4;
    pos4[ePos0] = state.pVector[0];
    pos4[ePos1] = state.pVector[1];
    pos4[ePos2] = state.pVector[2];
    pos4[eTime] = state.pVector[3];
    Acts::Vector3 dir;
    dir[eMom0] = state.pVector[4];
    dir[eMom1] = state.pVector[5];
    dir[eMom2] = state.pVector[6];
    const auto qOverP = state.pVector[7];

    std::optional<Covariance> covOpt = std::nullopt;
    if (state.covTransport && transportCov) {
      transportCovarianceToCurvilinear(state);
    }
    if (state.cov != Covariance::Zero()) {
      covOpt = state.cov;
    }

    CurvilinearTrackParameters parameters(pos4, dir, qOverP, std::move(covOpt),
                                          state.particleHypothesis);

    Jacobian jacobian(state.jacobian);

    return CurvilinearState(std::move(parameters), jacobian.transpose(),
                            state.pathAccumulated);
  }

  /// The state update method
  ///
  /// @param [in,out] state The stepper state for
  /// @param [in] parameters The new free track parameters at start
  /// @param [in] boundParams Corresponding bound parameters
  /// @param [in] covariance The updated covariance matrix
  /// @param [in] surface The surface used to update the pVector
  void update(State& state, const FreeVector& parameters,
              const BoundVector& boundParams, const Covariance& covariance,
              const Surface& surface) const {
    Vector3 direction = parameters.template segment<3>(eFreeDir0).normalized();
    state.pVector[0] = parameters[eFreePos0];
    state.pVector[1] = parameters[eFreePos1];
    state.pVector[2] = parameters[eFreePos2];
    state.pVector[3] = parameters[eFreeTime];
    state.pVector[4] = direction.x();
    state.pVector[5] = direction.y();
    state.pVector[6] = direction.z();
    state.pVector[7] = std::copysign(parameters[eFreeQOverP], state.pVector[7]);

    // @todo: remove magic numbers - is that the charge ?
    if (std::abs(state.pVector[7]) < .000000000000001) {
      state.pVector[7] < 0. ? state.pVector[7] = -.000000000000001
                            : state.pVector[7] = .000000000000001;
    }

    // prepare the jacobian if we have a covariance
    // copy the covariance matrix

    Vector3 pos(state.pVector[0], state.pVector[1], state.pVector[2]);
    Vector3 mom(state.pVector[4], state.pVector[5], state.pVector[6]);

    double Sf = std::sin(boundParams[eBoundPhi]);
    double Cf = std::cos(boundParams[eBoundPhi]);
    double Se = std::sin(boundParams[eBoundTheta]);
    double Ce = std::cos(boundParams[eBoundTheta]);

    const auto transform =
        surface.referenceFrame(state.options.geoContext, pos, mom);

    state.pVector[8] = transform(0, eBoundLoc0);
    state.pVector[16] = transform(0, eBoundLoc1);
    state.pVector[24] = 0.;
    state.pVector[32] = 0.;
    state.pVector[40] = 0.;
    state.pVector[48] = 0.;  // dX /

    state.pVector[9] = transform(1, eBoundLoc0);
    state.pVector[17] = transform(1, eBoundLoc1);
    state.pVector[25] = 0.;
    state.pVector[33] = 0.;
    state.pVector[41] = 0.;
    state.pVector[49] = 0.;  // dY /

    state.pVector[10] = transform(2, eBoundLoc0);
    state.pVector[18] = transform(2, eBoundLoc1);
    state.pVector[26] = 0.;
    state.pVector[34] = 0.;
    state.pVector[42] = 0.;
    state.pVector[50] = 0.;  // dZ /

    state.pVector[11] = 0.;
    state.pVector[19] = 0.;
    state.pVector[27] = 0.;
    state.pVector[35] = 0.;
    state.pVector[43] = 0.;
    state.pVector[51] = 1.;  // dT/

    state.pVector[12] = 0.;
    state.pVector[20] = 0.;
    state.pVector[28] = -Sf * Se;  // - sin(phi) * cos(theta)
    state.pVector[36] = Cf * Ce;   // cos(phi) * cos(theta)
    state.pVector[44] = 0.;
    state.pVector[52] = 0.;  // dAx/

    state.pVector[13] = 0.;
    state.pVector[21] = 0.;
    state.pVector[29] = Cf * Se;  // cos(phi) * sin(theta)
    state.pVector[37] = Sf * Ce;  // sin(phi) * cos(theta)
    state.pVector[45] = 0.;
    state.pVector[53] = 0.;  // dAy/

    state.pVector[14] = 0.;
    state.pVector[22] = 0.;
    state.pVector[30] = 0.;
    state.pVector[38] = -Se;  // - sin(theta)
    state.pVector[46] = 0.;
    state.pVector[54] = 0.;  // dAz/

    state.pVector[15] = 0.;
    state.pVector[23] = 0.;
    state.pVector[31] = 0.;
    state.pVector[39] = 0.;
    state.pVector[47] = 1.;
    state.pVector[55] = 0.;  // dCM/

    state.pVector[56] = 0.;
    state.pVector[57] = 0.;
    state.pVector[58] = 0.;
    state.pVector[59] = 0.;

    // special treatment for surface types
    // the disc needs polar coordinate adaptations
    if (surface.type() == Surface::Disc) {
      double lCf = std::cos(boundParams[eBoundLoc1]);
      double lSf = std::sin(boundParams[eBoundLoc1]);
      double Ax[3] = {transform(0, 0), transform(1, 0), transform(2, 0)};
      double Ay[3] = {transform(0, 1), transform(1, 1), transform(2, 1)};
      double d0 = lCf * Ax[0] + lSf * Ay[0];
      double d1 = lCf * Ax[1] + lSf * Ay[1];
      double d2 = lCf * Ax[2] + lSf * Ay[2];
      state.pVector[8] = d0;
      state.pVector[9] = d1;
      state.pVector[10] = d2;
      state.pVector[16] = boundParams[eBoundLoc0] * (lCf * Ay[0] - lSf * Ax[0]);
      state.pVector[17] = boundParams[eBoundLoc0] * (lCf * Ay[1] - lSf * Ax[1]);
      state.pVector[18] = boundParams[eBoundLoc0] * (lCf * Ay[2] - lSf * Ax[2]);
    }
    // the line needs components that relate direction change
    // with global frame change
    if (surface.type() == Surface::Perigee ||
        surface.type() == Surface::Straw) {
      // sticking to the nomenclature of the original RkPropagator
      // - axis pointing along the drift/transverse direction
      double B[3] = {transform(0, 0), transform(1, 0), transform(2, 0)};
      // - axis along the straw
      double A[3] = {transform(0, 1), transform(1, 1), transform(2, 1)};
      // - normal vector of the reference frame
      double C[3] = {transform(0, 2), transform(1, 2), transform(2, 2)};

      // projection of direction onto normal vector of reference frame
      double PC = state.pVector[4] * C[0] + state.pVector[5] * C[1] +
                  state.pVector[6] * C[2];
      double Bn = 1. / PC;

      double Bx2 = -A[2] * state.pVector[29];
      double Bx3 = A[1] * state.pVector[38] - A[2] * state.pVector[37];

      double By2 = A[2] * state.pVector[28];
      double By3 = A[2] * state.pVector[36] - A[0] * state.pVector[38];

      double Bz2 = A[0] * state.pVector[29] - A[1] * state.pVector[28];
      double Bz3 = A[0] * state.pVector[37] - A[1] * state.pVector[36];

      double B2 = B[0] * Bx2 + B[1] * By2 + B[2] * Bz2;
      double B3 = B[0] * Bx3 + B[1] * By3 + B[2] * Bz3;

      Bx2 = (Bx2 - B[0] * B2) * Bn;
      Bx3 = (Bx3 - B[0] * B3) * Bn;
      By2 = (By2 - B[1] * B2) * Bn;
      By3 = (By3 - B[1] * B3) * Bn;
      Bz2 = (Bz2 - B[2] * B2) * Bn;
      Bz3 = (Bz3 - B[2] * B3) * Bn;

      //  /dPhi      |     /dThe       |
      state.pVector[24] = Bx2 * boundParams[eBoundLoc0];
      state.pVector[32] = Bx3 * boundParams[eBoundLoc0];  // dX/
      state.pVector[25] = By2 * boundParams[eBoundLoc0];
      state.pVector[33] = By3 * boundParams[eBoundLoc0];  // dY/
      state.pVector[26] = Bz2 * boundParams[eBoundLoc0];
      state.pVector[34] = Bz3 * boundParams[eBoundLoc0];  // dZ/
    }

    state.covariance = new BoundSquareMatrix(covariance);
    state.covTransport = true;
    state.useJacobian = true;

    // declare the state as ready
    state.state_ready = true;
  }

  /// Method to update momentum, direction and p
  ///
  /// @param state The state object
  /// @param uposition the updated position
  /// @param udirection the updated direction
  /// @param qop the updated momentum value
  /// @param time the update time
  void update(State& state, const Vector3& uposition, const Vector3& udirection,
              double qop, double time) const {
    // update the vector
    state.pVector[0] = uposition[0];
    state.pVector[1] = uposition[1];
    state.pVector[2] = uposition[2];
    state.pVector[3] = time;
    state.pVector[4] = udirection[0];
    state.pVector[5] = udirection[1];
    state.pVector[6] = udirection[2];
    state.pVector[7] = qop;
  }

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state
  ///
  /// @param [in,out] state State of the stepper
  void transportCovarianceToCurvilinear(State& state) const {
    double P[60];
    for (unsigned int i = 0; i < 60; ++i) {
      P[i] = state.pVector[i];
    }

    double p = 1. / P[7];
    P[40] *= p;
    P[41] *= p;
    P[42] *= p;
    P[44] *= p;
    P[45] *= p;
    P[46] *= p;

    double An = std::sqrt(P[4] * P[4] + P[5] * P[5]);
    double Ax[3];
    if (An != 0.) {
      Ax[0] = -P[5] / An;
      Ax[1] = P[4] / An;
      Ax[2] = 0.;
    } else {
      Ax[0] = 1.;
      Ax[1] = 0.;
      Ax[2] = 0.;
    }

    double Ay[3] = {-Ax[1] * P[6], Ax[0] * P[6], An};
    double S[3] = {P[4], P[5], P[6]};

    double A = P[4] * S[0] + P[5] * S[1] + P[6] * S[2];
    if (A != 0.) {
      A = 1. / A;
    }
    S[0] *= A;
    S[1] *= A;
    S[2] *= A;

    double s0 = P[8] * S[0] + P[9] * S[1] + P[10] * S[2];
    double s1 = P[16] * S[0] + P[17] * S[1] + P[18] * S[2];
    double s2 = P[24] * S[0] + P[25] * S[1] + P[26] * S[2];
    double s3 = P[32] * S[0] + P[33] * S[1] + P[34] * S[2];
    double s4 = P[40] * S[0] + P[41] * S[1] + P[42] * S[2];

    P[8] -= (s0 * P[4]);
    P[9] -= (s0 * P[5]);
    P[10] -= (s0 * P[6]);
    P[11] -= (s0 * P[59]);
    P[12] -= (s0 * P[56]);
    P[13] -= (s0 * P[57]);
    P[14] -= (s0 * P[58]);
    P[16] -= (s1 * P[4]);
    P[17] -= (s1 * P[5]);
    P[18] -= (s1 * P[6]);
    P[19] -= (s1 * P[59]);
    P[20] -= (s1 * P[56]);
    P[21] -= (s1 * P[57]);
    P[22] -= (s1 * P[58]);
    P[24] -= (s2 * P[4]);
    P[25] -= (s2 * P[5]);
    P[26] -= (s2 * P[6]);
    P[27] -= (s2 * P[59]);
    P[28] -= (s2 * P[56]);
    P[29] -= (s2 * P[57]);
    P[30] -= (s2 * P[58]);
    P[32] -= (s3 * P[4]);
    P[33] -= (s3 * P[5]);
    P[34] -= (s3 * P[6]);
    P[35] -= (s3 * P[59]);
    P[36] -= (s3 * P[56]);
    P[37] -= (s3 * P[57]);
    P[38] -= (s3 * P[58]);
    P[40] -= (s4 * P[4]);
    P[41] -= (s4 * P[5]);
    P[42] -= (s4 * P[6]);
    P[43] -= (s4 * P[59]);
    P[44] -= (s4 * P[56]);
    P[45] -= (s4 * P[57]);
    P[46] -= (s4 * P[58]);

    double P3 = 0, P4 = 0, C = P[4] * P[4] + P[5] * P[5];
    if (C > 1.e-20) {
      C = 1. / C;
      P3 = P[4] * C;
      P4 = P[5] * C;
      C = -sqrt(C);
    } else {
      C = -1.e10;
      P3 = 1.;
      P4 = 0.;
    }

    // Jacobian production
    //
    state.jacobian[0] = Ax[0] * P[8] + Ax[1] * P[9];    // dL0/dL0
    state.jacobian[1] = Ax[0] * P[16] + Ax[1] * P[17];  // dL0/dL1
    state.jacobian[2] = Ax[0] * P[24] + Ax[1] * P[25];  // dL0/dPhi
    state.jacobian[3] = Ax[0] * P[32] + Ax[1] * P[33];  // dL0/dThe
    state.jacobian[4] = Ax[0] * P[40] + Ax[1] * P[41];  // dL0/dCM
    state.jacobian[5] = 0.;                             // dL0/dT

    state.jacobian[6] = Ay[0] * P[8] + Ay[1] * P[9] + Ay[2] * P[10];  // dL1/dL0
    state.jacobian[7] =
        Ay[0] * P[16] + Ay[1] * P[17] + Ay[2] * P[18];  // dL1/dL1
    state.jacobian[8] =
        Ay[0] * P[24] + Ay[1] * P[25] + Ay[2] * P[26];  // dL1/dPhi
    state.jacobian[9] =
        Ay[0] * P[32] + Ay[1] * P[33] + Ay[2] * P[34];  // dL1/dThe
    state.jacobian[10] =
        Ay[0] * P[40] + Ay[1] * P[41] + Ay[2] * P[42];  // dL1/dCM
    state.jacobian[11] = 0.;                            // dL1/dT

    state.jacobian[12] = P3 * P[13] - P4 * P[12];  // dPhi/dL0
    state.jacobian[13] = P3 * P[21] - P4 * P[20];  // dPhi/dL1
    state.jacobian[14] = P3 * P[29] - P4 * P[28];  // dPhi/dPhi
    state.jacobian[15] = P3 * P[37] - P4 * P[36];  // dPhi/dThe
    state.jacobian[16] = P3 * P[45] - P4 * P[44];  // dPhi/dCM
    state.jacobian[17] = 0.;                       // dPhi/dT

    state.jacobian[18] = C * P[14];  // dThe/dL0
    state.jacobian[19] = C * P[22];  // dThe/dL1
    state.jacobian[20] = C * P[30];  // dThe/dPhi
    state.jacobian[21] = C * P[38];  // dThe/dThe
    state.jacobian[22] = C * P[46];  // dThe/dCM
    state.jacobian[23] = 0.;         // dThe/dT

    state.jacobian[24] = 0.;     // dCM /dL0
    state.jacobian[25] = 0.;     // dCM /dL1
    state.jacobian[26] = 0.;     // dCM /dPhi
    state.jacobian[27] = 0.;     // dCM /dTheta
    state.jacobian[28] = P[47];  // dCM /dCM
    state.jacobian[29] = 0.;     // dCM/dT

    state.jacobian[30] = P[11];  // dT/dL0
    state.jacobian[31] = P[19];  // dT/dL1
    state.jacobian[32] = P[27];  // dT/dPhi
    state.jacobian[33] = P[35];  // dT/dThe
    state.jacobian[34] = P[43];  // dT/dCM
    state.jacobian[35] = P[51];  // dT/dT

    Eigen::Map<Eigen::Matrix<double, eBoundSize, eBoundSize, Eigen::RowMajor>>
        J(state.jacobian);
    state.cov = J * (*state.covariance) * J.transpose();
  }

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current position,
  /// or direction of the state
  ///
  /// @param [in,out] state State of the stepper
  /// @param [in] surface is the surface to which the covariance is forwarded to
  void transportCovarianceToBound(
      State& state, const Surface& surface,
      const FreeToBoundCorrection& /*freeToBoundCorrection*/ =
          FreeToBoundCorrection(false)) const {
    Acts::Vector3 gp(state.pVector[0], state.pVector[1], state.pVector[2]);
    Acts::Vector3 mom(state.pVector[4], state.pVector[5], state.pVector[6]);

    double P[60];
    for (unsigned int i = 0; i < 60; ++i) {
      P[i] = state.pVector[i];
    }

    mom /= std::abs(state.pVector[7]);

    double p = 1. / state.pVector[7];
    P[40] *= p;
    P[41] *= p;
    P[42] *= p;
    P[44] *= p;
    P[45] *= p;
    P[46] *= p;

    const auto fFrame =
        surface.referenceFrame(state.options.geoContext, gp, mom);

    double Ax[3] = {fFrame(0, 0), fFrame(1, 0), fFrame(2, 0)};
    double Ay[3] = {fFrame(0, 1), fFrame(1, 1), fFrame(2, 1)};
    double S[3] = {fFrame(0, 2), fFrame(1, 2), fFrame(2, 2)};

    // this is the projection of direction onto the local normal vector
    double A = P[4] * S[0] + P[5] * S[1] + P[6] * S[2];

    if (A != 0.) {
      A = 1. / A;
    }

    S[0] *= A;
    S[1] *= A;
    S[2] *= A;

    double s0 = P[8] * S[0] + P[9] * S[1] + P[10] * S[2];
    double s1 = P[16] * S[0] + P[17] * S[1] + P[18] * S[2];
    double s2 = P[24] * S[0] + P[25] * S[1] + P[26] * S[2];
    double s3 = P[32] * S[0] + P[33] * S[1] + P[34] * S[2];
    double s4 = P[40] * S[0] + P[41] * S[1] + P[42] * S[2];

    // in case of line-type surfaces - we need to take into account that
    // the reference frame changes with variations of all local
    // parameters
    if (surface.type() == Surface::Straw ||
        surface.type() == Surface::Perigee) {
      // vector from position to center
      double x = P[0] - surface.center(state.options.geoContext).x();
      double y = P[1] - surface.center(state.options.geoContext).y();
      double z = P[2] - surface.center(state.options.geoContext).z();

      // this is the projection of the direction onto the local y axis
      double d = P[4] * Ay[0] + P[5] * Ay[1] + P[6] * Ay[2];

      // this is cos(beta)
      double a = (1. - d) * (1. + d);
      if (a != 0.) {
        a = 1. / a;  // i.e. 1./(1-d^2)
      }

      // that's the modified norm vector
      double X = d * Ay[0] - P[4];  //
      double Y = d * Ay[1] - P[5];  //
      double Z = d * Ay[2] - P[6];  //

      // d0 to d1
      double d0 = P[12] * Ay[0] + P[13] * Ay[1] + P[14] * Ay[2];
      double d1 = P[20] * Ay[0] + P[21] * Ay[1] + P[22] * Ay[2];
      double d2 = P[28] * Ay[0] + P[29] * Ay[1] + P[30] * Ay[2];
      double d3 = P[36] * Ay[0] + P[37] * Ay[1] + P[38] * Ay[2];
      double d4 = P[44] * Ay[0] + P[45] * Ay[1] + P[46] * Ay[2];

      s0 = (((P[8] * X + P[9] * Y + P[10] * Z) + x * (d0 * Ay[0] - P[12])) +
            (y * (d0 * Ay[1] - P[13]) + z * (d0 * Ay[2] - P[14]))) *
           (-a);

      s1 = (((P[16] * X + P[17] * Y + P[18] * Z) + x * (d1 * Ay[0] - P[20])) +
            (y * (d1 * Ay[1] - P[21]) + z * (d1 * Ay[2] - P[22]))) *
           (-a);
      s2 = (((P[24] * X + P[25] * Y + P[26] * Z) + x * (d2 * Ay[0] - P[28])) +
            (y * (d2 * Ay[1] - P[29]) + z * (d2 * Ay[2] - P[30]))) *
           (-a);
      s3 = (((P[32] * X + P[33] * Y + P[34] * Z) + x * (d3 * Ay[0] - P[36])) +
            (y * (d3 * Ay[1] - P[37]) + z * (d3 * Ay[2] - P[38]))) *
           (-a);
      s4 = (((P[40] * X + P[41] * Y + P[42] * Z) + x * (d4 * Ay[0] - P[44])) +
            (y * (d4 * Ay[1] - P[45]) + z * (d4 * Ay[2] - P[46]))) *
           (-a);
    }

    P[8] -= (s0 * P[4]);
    P[9] -= (s0 * P[5]);
    P[10] -= (s0 * P[6]);
    P[11] -= (s0 * P[59]);
    P[12] -= (s0 * P[56]);
    P[13] -= (s0 * P[57]);
    P[14] -= (s0 * P[58]);

    P[16] -= (s1 * P[4]);
    P[17] -= (s1 * P[5]);
    P[18] -= (s1 * P[6]);
    P[19] -= (s1 * P[59]);
    P[20] -= (s1 * P[56]);
    P[21] -= (s1 * P[57]);
    P[22] -= (s1 * P[58]);

    P[24] -= (s2 * P[4]);
    P[25] -= (s2 * P[5]);
    P[26] -= (s2 * P[6]);
    P[27] -= (s2 * P[59]);
    P[28] -= (s2 * P[56]);
    P[29] -= (s2 * P[57]);
    P[30] -= (s2 * P[58]);

    P[32] -= (s3 * P[4]);
    P[33] -= (s3 * P[5]);
    P[34] -= (s3 * P[6]);
    P[35] -= (s3 * P[59]);
    P[36] -= (s3 * P[56]);
    P[37] -= (s3 * P[57]);
    P[38] -= (s3 * P[58]);

    P[40] -= (s4 * P[4]);
    P[41] -= (s4 * P[5]);
    P[42] -= (s4 * P[6]);
    P[43] -= (s4 * P[59]);
    P[44] -= (s4 * P[56]);
    P[45] -= (s4 * P[57]);
    P[46] -= (s4 * P[58]);

    double P3 = 0, P4 = 0, C = P[4] * P[4] + P[5] * P[5];
    if (C > 1.e-20) {
      C = 1. / C;
      P3 = P[4] * C;
      P4 = P[5] * C;
      C = -sqrt(C);
    } else {
      C = -1.e10;
      P3 = 1.;
      P4 = 0.;
    }

    double MA[3] = {Ax[0], Ax[1], Ax[2]};
    double MB[3] = {Ay[0], Ay[1], Ay[2]};
    // Jacobian production of transport and to_local
    if (surface.type() == Surface::Disc) {
      // the vector from the disc surface to the p
      const auto& sfc = surface.center(state.options.geoContext);
      double d[3] = {P[0] - sfc(0), P[1] - sfc(1), P[2] - sfc(2)};
      // this needs the transformation to polar coordinates
      double RC = d[0] * Ax[0] + d[1] * Ax[1] + d[2] * Ax[2];
      double RS = d[0] * Ay[0] + d[1] * Ay[1] + d[2] * Ay[2];
      double R2 = RC * RC + RS * RS;

      // inverse radius
      double Ri = 1. / sqrt(R2);
      MA[0] = (RC * Ax[0] + RS * Ay[0]) * Ri;
      MA[1] = (RC * Ax[1] + RS * Ay[1]) * Ri;
      MA[2] = (RC * Ax[2] + RS * Ay[2]) * Ri;
      MB[0] = (RC * Ay[0] - RS * Ax[0]) * (Ri = 1. / R2);
      MB[1] = (RC * Ay[1] - RS * Ax[1]) * Ri;
      MB[2] = (RC * Ay[2] - RS * Ax[2]) * Ri;
    }

    state.jacobian[0] = MA[0] * P[8] + MA[1] * P[9] + MA[2] * P[10];  // dL0/dL0
    state.jacobian[1] =
        MA[0] * P[16] + MA[1] * P[17] + MA[2] * P[18];  // dL0/dL1
    state.jacobian[2] =
        MA[0] * P[24] + MA[1] * P[25] + MA[2] * P[26];  // dL0/dPhi
    state.jacobian[3] =
        MA[0] * P[32] + MA[1] * P[33] + MA[2] * P[34];  // dL0/dThe
    state.jacobian[4] =
        MA[0] * P[40] + MA[1] * P[41] + MA[2] * P[42];  // dL0/dCM
    state.jacobian[5] = 0.;                             // dL0/dT

    state.jacobian[6] = MB[0] * P[8] + MB[1] * P[9] + MB[2] * P[10];  // dL1/dL0
    state.jacobian[7] =
        MB[0] * P[16] + MB[1] * P[17] + MB[2] * P[18];  // dL1/dL1
    state.jacobian[8] =
        MB[0] * P[24] + MB[1] * P[25] + MB[2] * P[26];  // dL1/dPhi
    state.jacobian[9] =
        MB[0] * P[32] + MB[1] * P[33] + MB[2] * P[34];  // dL1/dThe
    state.jacobian[10] =
        MB[0] * P[40] + MB[1] * P[41] + MB[2] * P[42];  // dL1/dCM
    state.jacobian[11] = 0.;                            // dL1/dT

    state.jacobian[12] = P3 * P[13] - P4 * P[12];  // dPhi/dL0
    state.jacobian[13] = P3 * P[21] - P4 * P[20];  // dPhi/dL1
    state.jacobian[14] = P3 * P[29] - P4 * P[28];  // dPhi/dPhi
    state.jacobian[15] = P3 * P[37] - P4 * P[36];  // dPhi/dThe
    state.jacobian[16] = P3 * P[45] - P4 * P[44];  // dPhi/dCM
    state.jacobian[17] = 0.;                       // dPhi/dT

    state.jacobian[18] = C * P[14];  // dThe/dL0
    state.jacobian[19] = C * P[22];  // dThe/dL1
    state.jacobian[20] = C * P[30];  // dThe/dPhi
    state.jacobian[21] = C * P[38];  // dThe/dThe
    state.jacobian[22] = C * P[46];  // dThe/dCM
    state.jacobian[23] = 0.;         // dThe/dT

    state.jacobian[24] = 0.;     // dCM /dL0
    state.jacobian[25] = 0.;     // dCM /dL1
    state.jacobian[26] = 0.;     // dCM /dPhi
    state.jacobian[27] = 0.;     // dCM /dTheta
    state.jacobian[28] = P[47];  // dCM /dCM
    state.jacobian[29] = 0.;     // dCM/dT

    state.jacobian[30] = P[11];  // dT/dL0
    state.jacobian[31] = P[19];  // dT/dL1
    state.jacobian[32] = P[27];  // dT/dPhi
    state.jacobian[33] = P[35];  // dT/dThe
    state.jacobian[34] = P[43];  // dT/dCM
    state.jacobian[35] = P[51];  // dT/dT

    Eigen::Map<Eigen::Matrix<double, eBoundSize, eBoundSize, Eigen::RowMajor>>
        J(state.jacobian);
    state.cov = J * (*state.covariance) * J.transpose();
  }

  /// Perform the actual step on the state
  ///
  /// @param state is the provided stepper state (caller keeps thread locality)
  /// @param propDir is the direction of propagation
  /// @param material is the material properties
  /// @return the result of the step
  Result<double> step(State& state, Direction propDir,
                      const IVolumeMaterial* material) const {
    (void)material;

    // we use h for keeping the nominclature with the original atlas code
    auto h = state.stepSize.value() * propDir;
    bool Jac = state.useJacobian;

    double* R = &(state.pVector[0]);  // Coordinates
    double* A = &(state.pVector[4]);  // Directions
    double* sA = &(state.pVector[56]);
    // Invert mometum/2.
    double Pi = 0.5 * state.pVector[7];
    //    double dltm = 0.0002 * .03;
    Vector3 f0, f;

    // if new field is required get it
    if (state.newfield) {
      const Vector3 pos(R[0], R[1], R[2]);
      // This is sd.B_first in EigenStepper
      auto fRes = getField(state, pos);
      if (!fRes.ok()) {
        return fRes.error();
      }
      f0 = *fRes;
    } else {
      f0 = state.field;
    }

    bool Helix = false;
    // if (std::abs(S) < m_cfg.helixStep) Helix = true;

    std::size_t nStepTrials = 0;
    while (h != 0.) {
      nStepTrials++;

      // PS2 is h/(2*momentum) in EigenStepper
      double S3 = (1. / 3.) * h, S4 = .25 * h, PS2 = Pi * h;

      // First point
      //
      // H0 is (h/(2*momentum) * sd.B_first) in EigenStepper
      double H0[3] = {f0[0] * PS2, f0[1] * PS2, f0[2] * PS2};
      // { A0, B0, C0 } is (h/2 * sd.k1) in EigenStepper
      double A0 = A[1] * H0[2] - A[2] * H0[1];
      double B0 = A[2] * H0[0] - A[0] * H0[2];
      double C0 = A[0] * H0[1] - A[1] * H0[0];
      // { A2, B2, C2 } is (h/2 * sd.k1 + direction) in EigenStepper
      double A2 = A0 + A[0];
      double B2 = B0 + A[1];
      double C2 = C0 + A[2];
      // { A1, B1, C1 } is (h/2 * sd.k1 + 2*direction) in EigenStepper
      double A1 = A2 + A[0];
      double B1 = B2 + A[1];
      double C1 = C2 + A[2];

      // Second point
      //
      if (!Helix) {
        // This is pos1 in EigenStepper
        const Vector3 pos(R[0] + A1 * S4, R[1] + B1 * S4, R[2] + C1 * S4);
        // This is sd.B_middle in EigenStepper
        auto fRes = getField(state, pos);
        if (!fRes.ok()) {
          return fRes.error();
        }
        f = *fRes;
      } else {
        f = f0;
      }

      // H1 is (h/(2*momentum) * sd.B_middle) in EigenStepper
      double H1[3] = {f[0] * PS2, f[1] * PS2, f[2] * PS2};
      // { A3, B3, C3 } is (direction + h/2 * sd.k2) in EigenStepper
      double A3 = (A[0] + B2 * H1[2]) - C2 * H1[1];
      double B3 = (A[1] + C2 * H1[0]) - A2 * H1[2];
      double C3 = (A[2] + A2 * H1[1]) - B2 * H1[0];
      // { A4, B4, C4 } is (direction + h/2 * sd.k3) in EigenStepper
      double A4 = (A[0] + B3 * H1[2]) - C3 * H1[1];
      double B4 = (A[1] + C3 * H1[0]) - A3 * H1[2];
      double C4 = (A[2] + A3 * H1[1]) - B3 * H1[0];
      // { A5, B5, C5 } is (direction + h * sd.k3) in EigenStepper
      double A5 = 2. * A4 - A[0];
      double B5 = 2. * B4 - A[1];
      double C5 = 2. * C4 - A[2];

      // Last point
      //
      if (!Helix) {
        // This is pos2 in EigenStepper
        const Vector3 pos(R[0] + h * A4, R[1] + h * B4, R[2] + h * C4);
        // This is sd.B_last in Eigen stepper
        auto fRes = getField(state, pos);
        if (!fRes.ok()) {
          return fRes.error();
        }
        f = *fRes;
      } else {
        f = f0;
      }

      // H2 is (h/(2*momentum) * sd.B_last) in EigenStepper
      double H2[3] = {f[0] * PS2, f[1] * PS2, f[2] * PS2};
      // { A6, B6, C6 } is (h/2 * sd.k4) in EigenStepper
      double A6 = B5 * H2[2] - C5 * H2[1];
      double B6 = C5 * H2[0] - A5 * H2[2];
      double C6 = A5 * H2[1] - B5 * H2[0];

      // Test approximation quality on give step and possible step reduction
      //
      // This is (h2 * (sd.k1 - sd.k2 - sd.k3 + sd.k4).template lpNorm<1>())
      // in EigenStepper
      double EST =
          2. * h *
          (std::abs((A1 + A6) - (A3 + A4)) + std::abs((B1 + B6) - (B3 + B4)) +
           std::abs((C1 + C6) - (C3 + C4)));
      if (std::abs(EST) > std::abs(state.options.stepTolerance)) {
        h = h * .5;
        // neutralize the sign of h again
        state.stepSize.setAccuracy(h * propDir);
        //        dltm = 0.;
        continue;
      }

      //      if (EST < dltm) h *= 2.;

      // Parameters calculation
      //
      double A00 = A[0], A11 = A[1], A22 = A[2];

      A[0] = 2. * A3 + (A0 + A5 + A6);
      A[1] = 2. * B3 + (B0 + B5 + B6);
      A[2] = 2. * C3 + (C0 + C5 + C6);

      double D = (A[0] * A[0] + A[1] * A[1]) + (A[2] * A[2] - 9.);
      double Sl = 2. / h;
      D = (1. / 3.) - ((1. / 648.) * D) * (12. - D);

      R[0] += (A2 + A3 + A4) * S3;
      R[1] += (B2 + B3 + B4) * S3;
      R[2] += (C2 + C3 + C4) * S3;
      A[0] *= D;
      A[1] *= D;
      A[2] *= D;
      sA[0] = A6 * Sl;
      sA[1] = B6 * Sl;
      sA[2] = C6 * Sl;

      double mass = particleHypothesis(state).mass();
      double momentum = absoluteMomentum(state);

      // Evaluate the time propagation
      double dtds = std::sqrt(1 + mass * mass / (momentum * momentum));
      state.pVector[3] += h * dtds;
      state.pVector[59] = dtds;
      state.field = f;
      state.newfield = false;

      if (Jac) {
        double dtdl = h * mass * mass * qOverP(state) / dtds;
        state.pVector[43] += dtdl;

        // Jacobian calculation
        //
        double* d2A = &state.pVector[28];
        double* d3A = &state.pVector[36];
        double* d4A = &state.pVector[44];
        double d2A0 = H0[2] * d2A[1] - H0[1] * d2A[2];
        double d2B0 = H0[0] * d2A[2] - H0[2] * d2A[0];
        double d2C0 = H0[1] * d2A[0] - H0[0] * d2A[1];
        double d3A0 = H0[2] * d3A[1] - H0[1] * d3A[2];
        double d3B0 = H0[0] * d3A[2] - H0[2] * d3A[0];
        double d3C0 = H0[1] * d3A[0] - H0[0] * d3A[1];
        double d4A0 = (A0 + H0[2] * d4A[1]) - H0[1] * d4A[2];
        double d4B0 = (B0 + H0[0] * d4A[2]) - H0[2] * d4A[0];
        double d4C0 = (C0 + H0[1] * d4A[0]) - H0[0] * d4A[1];
        double d2A2 = d2A0 + d2A[0];
        double d2B2 = d2B0 + d2A[1];
        double d2C2 = d2C0 + d2A[2];
        double d3A2 = d3A0 + d3A[0];
        double d3B2 = d3B0 + d3A[1];
        double d3C2 = d3C0 + d3A[2];
        double d4A2 = d4A0 + d4A[0];
        double d4B2 = d4B0 + d4A[1];
        double d4C2 = d4C0 + d4A[2];
        double d0 = d4A[0] - A00;
        double d1 = d4A[1] - A11;
        double d2 = d4A[2] - A22;
        double d2A3 = (d2A[0] + d2B2 * H1[2]) - d2C2 * H1[1];
        double d2B3 = (d2A[1] + d2C2 * H1[0]) - d2A2 * H1[2];
        double d2C3 = (d2A[2] + d2A2 * H1[1]) - d2B2 * H1[0];
        double d3A3 = (d3A[0] + d3B2 * H1[2]) - d3C2 * H1[1];
        double d3B3 = (d3A[1] + d3C2 * H1[0]) - d3A2 * H1[2];
        double d3C3 = (d3A[2] + d3A2 * H1[1]) - d3B2 * H1[0];
        double d4A3 = ((A3 + d0) + d4B2 * H1[2]) - d4C2 * H1[1];
        double d4B3 = ((B3 + d1) + d4C2 * H1[0]) - d4A2 * H1[2];
        double d4C3 = ((C3 + d2) + d4A2 * H1[1]) - d4B2 * H1[0];
        double d2A4 = (d2A[0] + d2B3 * H1[2]) - d2C3 * H1[1];
        double d2B4 = (d2A[1] + d2C3 * H1[0]) - d2A3 * H1[2];
        double d2C4 = (d2A[2] + d2A3 * H1[1]) - d2B3 * H1[0];
        double d3A4 = (d3A[0] + d3B3 * H1[2]) - d3C3 * H1[1];
        double d3B4 = (d3A[1] + d3C3 * H1[0]) - d3A3 * H1[2];
        double d3C4 = (d3A[2] + d3A3 * H1[1]) - d3B3 * H1[0];
        double d4A4 = ((A4 + d0) + d4B3 * H1[2]) - d4C3 * H1[1];
        double d4B4 = ((B4 + d1) + d4C3 * H1[0]) - d4A3 * H1[2];
        double d4C4 = ((C4 + d2) + d4A3 * H1[1]) - d4B3 * H1[0];
        double d2A5 = 2. * d2A4 - d2A[0];
        double d2B5 = 2. * d2B4 - d2A[1];
        double d2C5 = 2. * d2C4 - d2A[2];
        double d3A5 = 2. * d3A4 - d3A[0];
        double d3B5 = 2. * d3B4 - d3A[1];
        double d3C5 = 2. * d3C4 - d3A[2];
        double d4A5 = 2. * d4A4 - d4A[0];
        double d4B5 = 2. * d4B4 - d4A[1];
        double d4C5 = 2. * d4C4 - d4A[2];
        double d2A6 = d2B5 * H2[2] - d2C5 * H2[1];
        double d2B6 = d2C5 * H2[0] - d2A5 * H2[2];
        double d2C6 = d2A5 * H2[1] - d2B5 * H2[0];
        double d3A6 = d3B5 * H2[2] - d3C5 * H2[1];
        double d3B6 = d3C5 * H2[0] - d3A5 * H2[2];
        double d3C6 = d3A5 * H2[1] - d3B5 * H2[0];
        double d4A6 = d4B5 * H2[2] - d4C5 * H2[1];
        double d4B6 = d4C5 * H2[0] - d4A5 * H2[2];
        double d4C6 = d4A5 * H2[1] - d4B5 * H2[0];

        double* dR = &state.pVector[24];
        dR[0] += (d2A2 + d2A3 + d2A4) * S3;
        dR[1] += (d2B2 + d2B3 + d2B4) * S3;
        dR[2] += (d2C2 + d2C3 + d2C4) * S3;
        d2A[0] = ((d2A0 + 2. * d2A3) + (d2A5 + d2A6)) * (1. / 3.);
        d2A[1] = ((d2B0 + 2. * d2B3) + (d2B5 + d2B6)) * (1. / 3.);
        d2A[2] = ((d2C0 + 2. * d2C3) + (d2C5 + d2C6)) * (1. / 3.);

        dR = &state.pVector[32];
        dR[0] += (d3A2 + d3A3 + d3A4) * S3;
        dR[1] += (d3B2 + d3B3 + d3B4) * S3;
        dR[2] += (d3C2 + d3C3 + d3C4) * S3;
        d3A[0] = ((d3A0 + 2. * d3A3) + (d3A5 + d3A6)) * (1. / 3.);
        d3A[1] = ((d3B0 + 2. * d3B3) + (d3B5 + d3B6)) * (1. / 3.);
        d3A[2] = ((d3C0 + 2. * d3C3) + (d3C5 + d3C6)) * (1. / 3.);

        dR = &state.pVector[40];
        dR[0] += (d4A2 + d4A3 + d4A4) * S3;
        dR[1] += (d4B2 + d4B3 + d4B4) * S3;
        dR[2] += (d4C2 + d4C3 + d4C4) * S3;
        d4A[0] = ((d4A0 + 2. * d4A3) + (d4A5 + d4A6 + A6)) * (1. / 3.);
        d4A[1] = ((d4B0 + 2. * d4B3) + (d4B5 + d4B6 + B6)) * (1. / 3.);
        d4A[2] = ((d4C0 + 2. * d4C3) + (d4C5 + d4C6 + C6)) * (1. / 3.);
      }

      break;
    }

    state.pathAccumulated += h;
    ++state.nSteps;
    state.nStepTrials += nStepTrials;

    return h;
  }

  /// Method that reset the Jacobian to the Identity for when no bound state are
  /// available
  ///
  /// @param [in,out] state State of the stepper
  void setIdentityJacobian(State& state) const {
    state.jacobian[0] = 1.;  // dL0/dL0
    state.jacobian[1] = 0.;  // dL0/dL1
    state.jacobian[2] = 0.;  // dL0/dPhi
    state.jacobian[3] = 0.;  // dL0/dThe
    state.jacobian[4] = 0.;  // dL0/dCM
    state.jacobian[5] = 0.;  // dL0/dT

    state.jacobian[6] = 0.;   // dL1/dL0
    state.jacobian[7] = 1.;   // dL1/dL1
    state.jacobian[8] = 0.;   // dL1/dPhi
    state.jacobian[9] = 0.;   // dL1/dThe
    state.jacobian[10] = 0.;  // dL1/dCM
    state.jacobian[11] = 0.;  // dL1/dT

    state.jacobian[12] = 0.;  // dPhi/dL0
    state.jacobian[13] = 0.;  // dPhi/dL1
    state.jacobian[14] = 1.;  // dPhi/dPhi
    state.jacobian[15] = 0.;  // dPhi/dThe
    state.jacobian[16] = 0.;  // dPhi/dCM
    state.jacobian[17] = 0.;  // dPhi/dT

    state.jacobian[18] = 0.;  // dThe/dL0
    state.jacobian[19] = 0.;  // dThe/dL1
    state.jacobian[20] = 0.;  // dThe/dPhi
    state.jacobian[21] = 1.;  // dThe/dThe
    state.jacobian[22] = 0.;  // dThe/dCM
    state.jacobian[23] = 0.;  // dThe/dT

    state.jacobian[24] = 0.;  // dCM /dL0
    state.jacobian[25] = 0.;  // dCM /dL1
    state.jacobian[26] = 0.;  // dCM /dPhi
    state.jacobian[27] = 0.;  // dCM /dTheta
    state.jacobian[28] = 1.;  // dCM /dCM
    state.jacobian[29] = 0.;  // dCM/dT

    state.jacobian[30] = 0.;  // dT/dL0
    state.jacobian[31] = 0.;  // dT/dL1
    state.jacobian[32] = 0.;  // dT/dPhi
    state.jacobian[33] = 0.;  // dT/dThe
    state.jacobian[34] = 0.;  // dT/dCM
    state.jacobian[35] = 1.;  // dT/dT
  }

 private:
  std::shared_ptr<const MagneticFieldProvider> m_bField;

  /// Overstep limit: could/should be dynamic
  double m_overstepLimit = 100 * UnitConstants::um;
};

}  // namespace Acts
