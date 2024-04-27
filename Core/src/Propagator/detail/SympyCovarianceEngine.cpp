// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/SympyCovarianceEngine.hpp"

#include "Acts/Propagator/detail/JacobianEngine.hpp"

namespace {

template <typename T>
void transportCovarianceToBoundImpl(const T* C, const T* J_full, T* new_C) {
  const auto x0 = C[0] * J_full[0] + C[6] * J_full[6] + C[12] * J_full[12] +
                  C[18] * J_full[18] + C[24] * J_full[24];
  const auto x1 = C[6] * J_full[0] + C[7] * J_full[6] + C[13] * J_full[12] +
                  C[19] * J_full[18] + C[25] * J_full[24];
  const auto x2 = C[12] * J_full[0] + C[13] * J_full[6] + C[14] * J_full[12] +
                  C[20] * J_full[18] + C[26] * J_full[24];
  const auto x3 = C[18] * J_full[0] + C[19] * J_full[6] + C[20] * J_full[12] +
                  C[21] * J_full[18] + C[27] * J_full[24];
  const auto x4 = C[24] * J_full[0] + C[25] * J_full[6] + C[26] * J_full[12] +
                  C[27] * J_full[18] + C[28] * J_full[24];
  const auto x5 = C[0] * J_full[1] + C[6] * J_full[7] + C[12] * J_full[13] +
                  C[18] * J_full[19] + C[24] * J_full[25];
  const auto x6 = C[6] * J_full[1] + C[7] * J_full[7] + C[13] * J_full[13] +
                  C[19] * J_full[19] + C[25] * J_full[25];
  const auto x7 = C[12] * J_full[1] + C[13] * J_full[7] + C[14] * J_full[13] +
                  C[20] * J_full[19] + C[26] * J_full[25];
  const auto x8 = C[18] * J_full[1] + C[19] * J_full[7] + C[20] * J_full[13] +
                  C[21] * J_full[19] + C[27] * J_full[25];
  const auto x9 = C[24] * J_full[1] + C[25] * J_full[7] + C[26] * J_full[13] +
                  C[27] * J_full[19] + C[28] * J_full[25];
  const auto x10 = C[0] * J_full[2] + C[6] * J_full[8] + C[12] * J_full[14] +
                   C[18] * J_full[20] + C[24] * J_full[26];
  const auto x11 = C[6] * J_full[2] + C[7] * J_full[8] + C[13] * J_full[14] +
                   C[19] * J_full[20] + C[25] * J_full[26];
  const auto x12 = C[12] * J_full[2] + C[13] * J_full[8] + C[14] * J_full[14] +
                   C[20] * J_full[20] + C[26] * J_full[26];
  const auto x13 = C[18] * J_full[2] + C[19] * J_full[8] + C[20] * J_full[14] +
                   C[21] * J_full[20] + C[27] * J_full[26];
  const auto x14 = C[24] * J_full[2] + C[25] * J_full[8] + C[26] * J_full[14] +
                   C[27] * J_full[20] + C[28] * J_full[26];
  const auto x15 = C[0] * J_full[3] + C[6] * J_full[9] + C[12] * J_full[15] +
                   C[18] * J_full[21] + C[24] * J_full[27];
  const auto x16 = C[6] * J_full[3] + C[7] * J_full[9] + C[13] * J_full[15] +
                   C[19] * J_full[21] + C[25] * J_full[27];
  const auto x17 = C[12] * J_full[3] + C[13] * J_full[9] + C[14] * J_full[15] +
                   C[20] * J_full[21] + C[26] * J_full[27];
  const auto x18 = C[18] * J_full[3] + C[19] * J_full[9] + C[20] * J_full[15] +
                   C[21] * J_full[21] + C[27] * J_full[27];
  const auto x19 = C[24] * J_full[3] + C[25] * J_full[9] + C[26] * J_full[15] +
                   C[27] * J_full[21] + C[28] * J_full[27];
  const auto x20 = C[24] * J_full[5] + C[25] * J_full[11] + C[26] * J_full[17] +
                   C[27] * J_full[23] + C[28] * J_full[29] + C[34];
  const auto x21 = C[0] * J_full[5] + C[6] * J_full[11] + C[12] * J_full[17] +
                   C[18] * J_full[23] + C[24] * J_full[29] + C[30];
  const auto x22 = C[6] * J_full[5] + C[7] * J_full[11] + C[13] * J_full[17] +
                   C[19] * J_full[23] + C[25] * J_full[29] + C[31];
  const auto x23 = C[12] * J_full[5] + C[13] * J_full[11] + C[14] * J_full[17] +
                   C[20] * J_full[23] + C[26] * J_full[29] + C[32];
  const auto x24 = C[18] * J_full[5] + C[19] * J_full[11] + C[20] * J_full[17] +
                   C[21] * J_full[23] + C[27] * J_full[29] + C[33];
  new_C[0] = x0 * J_full[0] + x1 * J_full[6] + x2 * J_full[12] +
             x3 * J_full[18] + x4 * J_full[24];
  new_C[1] = x5 * J_full[0] + x6 * J_full[6] + x7 * J_full[12] +
             x8 * J_full[18] + x9 * J_full[24];
  new_C[2] = x10 * J_full[0] + x11 * J_full[6] + x12 * J_full[12] +
             x13 * J_full[18] + x14 * J_full[24];
  new_C[3] = x15 * J_full[0] + x16 * J_full[6] + x17 * J_full[12] +
             x18 * J_full[18] + x19 * J_full[24];
  new_C[4] = x4;
  new_C[5] = x20 * J_full[24] + x21 * J_full[0] + x22 * J_full[6] +
             x23 * J_full[12] + x24 * J_full[18];
  new_C[6] = x0 * J_full[1] + x1 * J_full[7] + x2 * J_full[13] +
             x3 * J_full[19] + x4 * J_full[25];
  new_C[7] = x5 * J_full[1] + x6 * J_full[7] + x7 * J_full[13] +
             x8 * J_full[19] + x9 * J_full[25];
  new_C[8] = x10 * J_full[1] + x11 * J_full[7] + x12 * J_full[13] +
             x13 * J_full[19] + x14 * J_full[25];
  new_C[9] = x15 * J_full[1] + x16 * J_full[7] + x17 * J_full[13] +
             x18 * J_full[19] + x19 * J_full[25];
  new_C[10] = x9;
  new_C[11] = x20 * J_full[25] + x21 * J_full[1] + x22 * J_full[7] +
              x23 * J_full[13] + x24 * J_full[19];
  new_C[12] = x0 * J_full[2] + x1 * J_full[8] + x2 * J_full[14] +
              x3 * J_full[20] + x4 * J_full[26];
  new_C[13] = x5 * J_full[2] + x6 * J_full[8] + x7 * J_full[14] +
              x8 * J_full[20] + x9 * J_full[26];
  new_C[14] = x10 * J_full[2] + x11 * J_full[8] + x12 * J_full[14] +
              x13 * J_full[20] + x14 * J_full[26];
  new_C[15] = x15 * J_full[2] + x16 * J_full[8] + x17 * J_full[14] +
              x18 * J_full[20] + x19 * J_full[26];
  new_C[16] = x14;
  new_C[17] = x20 * J_full[26] + x21 * J_full[2] + x22 * J_full[8] +
              x23 * J_full[14] + x24 * J_full[20];
  new_C[18] = x0 * J_full[3] + x1 * J_full[9] + x2 * J_full[15] +
              x3 * J_full[21] + x4 * J_full[27];
  new_C[19] = x5 * J_full[3] + x6 * J_full[9] + x7 * J_full[15] +
              x8 * J_full[21] + x9 * J_full[27];
  new_C[20] = x10 * J_full[3] + x11 * J_full[9] + x12 * J_full[15] +
              x13 * J_full[21] + x14 * J_full[27];
  new_C[21] = x15 * J_full[3] + x16 * J_full[9] + x17 * J_full[15] +
              x18 * J_full[21] + x19 * J_full[27];
  new_C[22] = x19;
  new_C[23] = x20 * J_full[27] + x21 * J_full[3] + x22 * J_full[9] +
              x23 * J_full[15] + x24 * J_full[21];
  new_C[24] = x4;
  new_C[25] = x9;
  new_C[26] = x14;
  new_C[27] = x19;
  new_C[28] = C[28];
  new_C[29] = x20;
  new_C[30] = x0 * J_full[5] + x1 * J_full[11] + x2 * J_full[17] +
              x3 * J_full[23] + x4 * J_full[29] + C[30] * J_full[0] +
              C[31] * J_full[6] + C[32] * J_full[12] + C[33] * J_full[18] +
              C[34] * J_full[24];
  new_C[31] = x5 * J_full[5] + x6 * J_full[11] + x7 * J_full[17] +
              x8 * J_full[23] + x9 * J_full[29] + C[30] * J_full[1] +
              C[31] * J_full[7] + C[32] * J_full[13] + C[33] * J_full[19] +
              C[34] * J_full[25];
  new_C[32] = x10 * J_full[5] + x11 * J_full[11] + x12 * J_full[17] +
              x13 * J_full[23] + x14 * J_full[29] + C[30] * J_full[2] +
              C[31] * J_full[8] + C[32] * J_full[14] + C[33] * J_full[20] +
              C[34] * J_full[26];
  new_C[33] = x15 * J_full[5] + x16 * J_full[11] + x17 * J_full[17] +
              x18 * J_full[23] + x19 * J_full[29] + C[30] * J_full[3] +
              C[31] * J_full[9] + C[32] * J_full[15] + C[33] * J_full[21] +
              C[34] * J_full[27];
  new_C[34] = x20;
  new_C[35] = x20 * J_full[29] + x21 * J_full[5] + x22 * J_full[11] +
              x23 * J_full[17] + x24 * J_full[23] + C[30] * J_full[5] +
              C[31] * J_full[11] + C[32] * J_full[17] + C[33] * J_full[23] +
              C[34] * J_full[29] + C[35];
}

}  // namespace

namespace Acts::detail {

/// Some type defs
using Jacobian = BoundMatrix;
using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;
using CurvilinearState =
    std::tuple<CurvilinearTrackParameters, Jacobian, double>;

Result<BoundState> sympy::boundState(
    const GeometryContext& geoContext, const Surface& surface,
    BoundSquareMatrix& boundCovariance, BoundMatrix& fullTransportJacobian,
    FreeMatrix& freeTransportJacobian, FreeVector& freeToPathDerivatives,
    BoundToFreeMatrix& boundToFreeJacobian, FreeVector& freeParameters,
    const ParticleHypothesis& particleHypothesis, bool covTransport,
    double accumulatedPath,
    const FreeToBoundCorrection& freeToBoundCorrection) {
  // Create the bound parameters
  Result<BoundVector> bv =
      transformFreeToBoundParameters(freeParameters, surface, geoContext);
  if (!bv.ok()) {
    return bv.error();
  }

  // Covariance transport
  std::optional<BoundSquareMatrix> cov = std::nullopt;
  if (covTransport) {
    // Calculate the jacobian and transport the covarianceMatrix to final local.
    // Then reinitialize the transportJacobian, derivatives and the
    // boundToFreeJacobian
    transportCovarianceToBound(geoContext, surface, boundCovariance,
                               fullTransportJacobian, freeTransportJacobian,
                               freeToPathDerivatives, boundToFreeJacobian,
                               freeParameters, freeToBoundCorrection);
    cov = boundCovariance;
  }

  // Create the bound state
  return std::make_tuple(
      BoundTrackParameters(surface.getSharedPtr(), *bv, std::move(cov),
                           particleHypothesis),
      fullTransportJacobian, accumulatedPath);
}

CurvilinearState sympy::curvilinearState(
    BoundSquareMatrix& boundCovariance, BoundMatrix& fullTransportJacobian,
    FreeMatrix& freeTransportJacobian, FreeVector& freeToPathDerivatives,
    BoundToFreeMatrix& boundToFreeJacobian, const FreeVector& freeParameters,
    const ParticleHypothesis& particleHypothesis, bool covTransport,
    double accumulatedPath) {
  const Vector3& direction = freeParameters.segment<3>(eFreeDir0);

  // Covariance transport
  std::optional<BoundSquareMatrix> cov = std::nullopt;
  if (covTransport) {
    // Calculate the jacobian and transport the covarianceMatrix to final local.
    // Then reinitialize the transportJacobian, derivatives and the
    // boundToFreeJacobian
    transportCovarianceToCurvilinear(
        boundCovariance, fullTransportJacobian, freeTransportJacobian,
        freeToPathDerivatives, boundToFreeJacobian, direction);
    cov = boundCovariance;
  }

  // Create the curvilinear parameters
  Vector4 pos4 = Vector4::Zero();
  pos4[ePos0] = freeParameters[eFreePos0];
  pos4[ePos1] = freeParameters[eFreePos1];
  pos4[ePos2] = freeParameters[eFreePos2];
  pos4[eTime] = freeParameters[eFreeTime];
  CurvilinearTrackParameters curvilinearParams(
      pos4, direction, freeParameters[eFreeQOverP], std::move(cov),
      particleHypothesis);
  // Create the curvilinear state
  return std::make_tuple(std::move(curvilinearParams), fullTransportJacobian,
                         accumulatedPath);
}

void sympy::transportCovarianceToBound(
    const GeometryContext& geoContext, const Surface& surface,
    BoundSquareMatrix& boundCovariance, BoundMatrix& fullTransportJacobian,
    FreeMatrix& freeTransportJacobian, FreeVector& freeToPathDerivatives,
    BoundToFreeMatrix& boundToFreeJacobian, FreeVector& freeParameters,
    const FreeToBoundCorrection& freeToBoundCorrection) {
  // Calculate the full jacobian from local parameters at the start surface to
  // current bound parameters
  boundToBoundTransportJacobian(geoContext, surface, freeParameters,
                                boundToFreeJacobian, freeTransportJacobian,
                                freeToPathDerivatives, fullTransportJacobian);

  bool correction = false;
  if (freeToBoundCorrection) {
    BoundToFreeMatrix startBoundToFinalFreeJacobian =
        freeTransportJacobian * boundToFreeJacobian;
    FreeSquareMatrix freeCovariance = startBoundToFinalFreeJacobian *
                                      boundCovariance *
                                      startBoundToFinalFreeJacobian.transpose();

    auto transformer =
        detail::CorrectedFreeToBoundTransformer(freeToBoundCorrection);
    auto correctedRes =
        transformer(freeParameters, freeCovariance, surface, geoContext);

    if (correctedRes.has_value()) {
      auto correctedValue = correctedRes.value();
      BoundVector boundParams = std::get<BoundVector>(correctedValue);
      // 1. Update the free parameters with the corrected bound parameters
      freeParameters =
          transformBoundToFreeParameters(surface, geoContext, boundParams);

      // 2. Update the bound covariance
      boundCovariance = std::get<BoundSquareMatrix>(correctedValue);

      correction = true;
    }
  }

  if (!correction) {
    // Apply the actual covariance transport to get covariance of the current
    // bound parameters
    BoundMatrix newBoundCovariance;
    transportCovarianceToBoundImpl(boundCovariance.data(),
                                   fullTransportJacobian.data(),
                                   newBoundCovariance.data());
    boundCovariance = newBoundCovariance;
  }

  // Reinitialize jacobian components:
  // ->The transportJacobian is reinitialized to Identity
  // ->The derivatives is reinitialized to Zero
  // ->The boundToFreeJacobian is initialized to that at the current surface
  reinitializeJacobians(geoContext, surface, freeTransportJacobian,
                        freeToPathDerivatives, boundToFreeJacobian,
                        freeParameters);
}

void sympy::transportCovarianceToCurvilinear(
    BoundSquareMatrix& boundCovariance, BoundMatrix& fullTransportJacobian,
    FreeMatrix& freeTransportJacobian, FreeVector& freeToPathDerivatives,
    BoundToFreeMatrix& boundToFreeJacobian, const Vector3& direction) {
  // Calculate the full jacobian from local parameters at the start surface to
  // current curvilinear parameters
  boundToCurvilinearTransportJacobian(
      direction, boundToFreeJacobian, freeTransportJacobian,
      freeToPathDerivatives, fullTransportJacobian);

  // Apply the actual covariance transport to get covariance of the current
  // curvilinear parameters
  BoundMatrix newBoundCovariance;
  transportCovarianceToBoundImpl(boundCovariance.data(),
                                 fullTransportJacobian.data(),
                                 newBoundCovariance.data());
  boundCovariance = newBoundCovariance;

  // Reinitialize jacobian components:
  // ->The free transportJacobian is reinitialized to Identity
  // ->The path derivatives is reinitialized to Zero
  // ->The boundToFreeJacobian is reinitialized to that at the current
  // curvilinear surface
  reinitializeJacobians(freeTransportJacobian, freeToPathDerivatives,
                        boundToFreeJacobian, direction);
}

}  // namespace Acts::detail
