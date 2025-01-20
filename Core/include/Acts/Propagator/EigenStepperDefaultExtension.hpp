// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <array>

namespace Acts {

class IVolumeMaterial;

/// @brief Default evaluator of the k_i's and elements of the transport matrix
/// D of the RKN4 stepping. This is a pure implementation by textbook.
struct EigenStepperDefaultExtension {
  /// @brief Evaluator of the k_i's of the RKN4. For the case of i = 0 this
  /// step sets up qop, too.
  ///
  /// @tparam i Index of the k_i, i = [0, 3]
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param [in] state State of the stepper
  /// @param [in] stepper Stepper of the propagation
  /// @param [in] volumeMaterial Material of the volume
  /// @param [out] knew Next k_i that is evaluated
  /// @param [in] bField B-Field at the evaluation position
  /// @param [out] kQoP k_i elements of the momenta
  /// @param [in] h Step size (= 0. ^ 0.5 * StepSize ^ StepSize)
  /// @param [in] kprev Evaluated k_{i - 1}
  ///
  /// @return Boolean flag if the calculation is valid
  template <int i, typename stepper_t>
  bool k(const stepper_t::State& state, const stepper_t& stepper,
         const IVolumeMaterial* volumeMaterial, Vector3& knew,
         const Vector3& bField, std::array<double, 4>& kQoP,
         const double h = 0., const Vector3& kprev = Vector3::Zero())
    requires(i >= 0 && i <= 3)
  {
    (void)volumeMaterial;

    auto qop = stepper.qOverP(state);
    // First step does not rely on previous data
    if constexpr (i == 0) {
      knew = qop * stepper.direction(state).cross(bField);
      kQoP = {0., 0., 0., 0.};
    } else {
      knew = qop * (stepper.direction(state) + h * kprev).cross(bField);
    }
    return true;
  }

  /// @brief Veto function after a RKN4 step was accepted by judging on the
  /// error of the step. Since the textbook does not deliver further vetos,
  /// this is a dummy function.
  ///
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param [in] state State of the stepper
  /// @param [in] stepper Stepper of the propagation
  /// @param [in] volumeMaterial Material of the volume
  /// @param [in] h Step size
  ///
  /// @return Boolean flag if the calculation is valid
  template <typename stepper_t>
  bool finalize(stepper_t::State& state, const stepper_t& stepper,
                const IVolumeMaterial* volumeMaterial, const double h) const {
    (void)volumeMaterial;

    propagateTime(state, stepper, h);
    return true;
  }

  /// @brief Veto function after a RKN4 step was accepted by judging on the
  /// error of the step. Since the textbook does not deliver further vetos,
  /// this is just for the evaluation of the transport matrix.
  ///
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param [in] state State of the stepper
  /// @param [in] stepper Stepper of the propagation
  /// @param [in] volumeMaterial Material of the volume
  /// @param [in] h Step size
  /// @param [out] D Transport matrix
  ///
  /// @return Boolean flag if the calculation is valid
  template <typename stepper_t>
  bool finalize(stepper_t::State& state, const stepper_t& stepper,
                const IVolumeMaterial* volumeMaterial, const double h,
                FreeMatrix& D) const {
    (void)volumeMaterial;

    propagateTime(state, stepper, h);
    return transportMatrix(state, stepper, h, D);
  }

 private:
  /// @brief Propagation function for the time coordinate
  ///
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param [in, out] state State of the stepper
  /// @param [in] stepper Stepper of the propagation
  /// @param [in] h Step size
  template <typename stepper_t>
  void propagateTime(stepper_t::State& state, const stepper_t& stepper,
                     const double h) const {
    /// This evaluation is based on dt/ds = 1/v = 1/(beta * c) with the velocity
    /// v, the speed of light c and beta = v/c. This can be re-written as dt/ds
    /// = sqrt(m^2/p^2 + c^{-2}) with the mass m and the momentum p.
    auto m = stepper.particleHypothesis(state).mass();
    auto p = stepper.absoluteMomentum(state);
    auto dtds = std::sqrt(1 + m * m / (p * p));
    state.pars[eFreeTime] += h * dtds;
    if (state.covTransport) {
      state.derivative(3) = dtds;
    }
  }

  /// @brief Calculates the transport matrix D for the jacobian
  ///
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param [in] state State of the stepper
  /// @param [in] stepper Stepper of the propagation
  /// @param [in] h Step size
  /// @param [out] D Transport matrix
  ///
  /// @return Boolean flag if evaluation is valid
  template <typename stepper_t>
  bool transportMatrix(stepper_t::State& state, const stepper_t& stepper,
                       const double h, FreeMatrix& D) const {
    /// The calculations are based on ATL-SOFT-PUB-2009-002. The update of the
    /// Jacobian matrix is requires only the calculation of eq. 17 and 18.
    /// Since the terms of eq. 18 are currently 0, this matrix is not needed
    /// in the calculation. The matrix A from eq. 17 consists out of 3
    /// different parts. The first one is given by the upper left 3x3 matrix
    /// that are calculated by the derivatives dF/dT (called dFdT) and dG/dT
    /// (calls dGdT). The second is given by the top 3 lines of the rightmost
    /// column. This is calculated by dFdL and dGdL. The remaining non-zero term
    /// is calculated directly. The naming of the variables is explained in eq.
    /// 11 and are directly related to the initial problem in eq. 7.
    /// The evaluation is based by propagating the parameters T and lambda as
    /// given in eq. 16 and evaluating the derivations for matrix A.
    /// @note The translation for u_{n+1} in eq. 7 is in this case a
    /// 3-dimensional vector without a dependency of Lambda or lambda neither in
    /// u_n nor in u_n'. The second and fourth eq. in eq. 14 have the constant
    /// offset matrices h * Id and Id respectively. This involves that the
    /// constant offset does not exist for rectangular matrix dGdu' (due to the
    /// missing Lambda part) and only exists for dFdu' in dlambda/dlambda.

    auto m = state.particleHypothesis.mass();
    auto& sd = state.stepData;
    auto dir = stepper.direction(state);
    auto qop = stepper.qOverP(state);
    auto p = stepper.absoluteMomentum(state);
    auto dtds = std::sqrt(1 + m * m / (p * p));

    D = FreeMatrix::Identity();

    double half_h = h * 0.5;
    // This sets the reference to the sub matrices
    // dFdx is already initialised as (3x3) idendity
    auto dFdT = D.block<3, 3>(0, 4);
    auto dFdL = D.block<3, 1>(0, 7);
    // dGdx is already initialised as (3x3) zero
    auto dGdT = D.block<3, 3>(4, 4);
    auto dGdL = D.block<3, 1>(4, 7);

    SquareMatrix3 dk1dT = SquareMatrix3::Zero();
    SquareMatrix3 dk2dT = SquareMatrix3::Identity();
    SquareMatrix3 dk3dT = SquareMatrix3::Identity();
    SquareMatrix3 dk4dT = SquareMatrix3::Identity();

    Vector3 dk1dL = Vector3::Zero();
    Vector3 dk2dL = Vector3::Zero();
    Vector3 dk3dL = Vector3::Zero();
    Vector3 dk4dL = Vector3::Zero();

    // For the case without energy loss
    dk1dL = dir.cross(sd.B_first);
    dk2dL = (dir + half_h * sd.k1).cross(sd.B_middle) +
            qop * half_h * dk1dL.cross(sd.B_middle);
    dk3dL = (dir + half_h * sd.k2).cross(sd.B_middle) +
            qop * half_h * dk2dL.cross(sd.B_middle);
    dk4dL =
        (dir + h * sd.k3).cross(sd.B_last) + qop * h * dk3dL.cross(sd.B_last);

    dk1dT(0, 1) = sd.B_first.z();
    dk1dT(0, 2) = -sd.B_first.y();
    dk1dT(1, 0) = -sd.B_first.z();
    dk1dT(1, 2) = sd.B_first.x();
    dk1dT(2, 0) = sd.B_first.y();
    dk1dT(2, 1) = -sd.B_first.x();
    dk1dT *= qop;

    dk2dT += half_h * dk1dT;
    dk2dT = qop * VectorHelpers::cross(dk2dT, sd.B_middle);

    dk3dT += half_h * dk2dT;
    dk3dT = qop * VectorHelpers::cross(dk3dT, sd.B_middle);

    dk4dT += h * dk3dT;
    dk4dT = qop * VectorHelpers::cross(dk4dT, sd.B_last);

    dFdT.setIdentity();
    dFdT += h / 6. * (dk1dT + dk2dT + dk3dT);
    dFdT *= h;

    dFdL = (h * h) / 6. * (dk1dL + dk2dL + dk3dL);

    dGdT += h / 6. * (dk1dT + 2. * (dk2dT + dk3dT) + dk4dT);

    dGdL = h / 6. * (dk1dL + 2. * (dk2dL + dk3dL) + dk4dL);

    D(3, 7) = h * m * m * qop / dtds;
    return true;
  }
};

}  // namespace Acts
