// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"

#include <type_traits>

namespace Acts::detail {

/// @brief check and correct covariance matrix
///
/// @tparam CovMatrix_t The type of covariance matrix
/// @tparam NumIter The number of iterations to run the correction
///
/// Invocation:
///   - covariance_helper<CovMatrix_t, numIter>::validate(covariance)
///    The 'covariance' is checked against semi-positivedefiniteness
///    and limited number of iterations to replace it with the
///    closest semi-positivedefinite are made if it's not
///
/// @return The (corrected) covariance is semi-positivedefinite or not
template <typename CovMatrix_t, signed int NumIter = 1>
struct CovarianceHelper {
  /// check if the covariance is semi-positive and correction is attempted
  static bool validate(CovMatrix_t& covariance) {
    if (covariance.hasNaN()) {
      return false;
    }
    std::size_t nIteration = 0;
    while (nIteration < NumIter) {
      if (isSemiPositive(covariance)) {
        return true;
      } else {
        Eigen::JacobiSVD<CovMatrix_t> svdCov(
            covariance, Eigen::ComputeFullU | Eigen::ComputeFullV);
        CovMatrix_t S = svdCov.singularValues().asDiagonal();
        CovMatrix_t V = svdCov.matrixV();
        CovMatrix_t H = V * S * V.transpose();
        covariance = (covariance + H) / 2;
        nIteration++;
      }
    }
    /// check again after the iterations
    return isSemiPositive(covariance);
  }

  /// check if the covariance is semi-positive
  static bool isSemiPositive(const CovMatrix_t& covariance) {
    if (covariance.hasNaN()) {
      return false;
    }
    Eigen::LDLT<CovMatrix_t> ldltCov(covariance);
    return ldltCov.isPositive();
  }

  /// check if the covariance is positive
  static bool isPositive(const CovMatrix_t& covariance) {
    if (covariance.hasNaN()) {
      return false;
    }
    Eigen::LLT<CovMatrix_t> lltCov(covariance);
    return lltCov.info() == Eigen::Success;
  }
};

}  // namespace Acts::detail
