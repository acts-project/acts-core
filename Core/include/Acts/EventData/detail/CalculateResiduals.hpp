// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/detail/ParameterTraits.hpp"

#include <cassert>

#include <Eigen/Core>

namespace Acts {
namespace detail {

/// Residuals between bound reference parameters and a measured subspace.
///
/// @tparam index_container_t SequenceContainer for measured inidices
/// @tparam measured_t Measured parameters vector type
/// @tparam residuals_t Residuals vector type
/// @param[in] size Size of the measured parameters subspace
/// @param[in] indices Indices of measured subspace, must have `size` entries
/// @param[in] reference Reference bound parameters
/// @param[in] measured Measured parameters subspace
/// @param[out] residuals Resulting residuals in the measured subspace
///
/// @note The separate `size` parameter is also used to allow the selection of
///   the correct residuals methods depending on the parameters type.
template <typename index_container_t, typename measured_t, typename residuals_t>
inline void calculateResiduals(BoundIndices size,
                               const index_container_t& indices,
                               const BoundVector& reference,
                               const Eigen::MatrixBase<measured_t>& measured,
                               Eigen::MatrixBase<residuals_t>& residuals) {
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(measured_t);
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(residuals_t);
  assert((size <= eBoundSize) and "Measured subspace is too large");
  assert((size <= measured.size()) and "Inconsistent measured size");
  assert((size <= residuals.size()) and "Inconsistent residuals size");

  for (size_t i = 0; i < size; ++i) {
    // this is neither elegant nor smart but it is the simplest solution.
    //
    // only phi must be handled specially here. the theta limits can only be
    // correctly handled if phi is updated, too. since we can not ensure that
    // both are available, it is probably less error-prone to treat theta as a
    // regular, unrestricted parameter here.
    if (i == eBoundPhi) {
      residuals[i] = ParameterTraits<BoundIndices, eBoundPhi>::getDifference(
          measured[i], reference[indices[i]]);
    } else {
      residuals[i] = measured[i] - reference[indices[i]];
    }
  }
}

/// Residuals between free reference parameters and a measured subspace.
///
/// @tparam index_container_t SequenceContainer for measured inidices
/// @tparam measured_t Measured parameters vector type
/// @tparam residuals_t Residuals vector type
/// @param[in] size Size of the measured parameters subspace
/// @param[in] indices Indices of measured subspace, must have `size` entries
/// @param[in] reference Reference free parameters
/// @param[in] measured Measured parameters subspace
/// @param[out] residuals Resulting residuals in the measured subspace
///
/// @note The separate `size` parameter is also used to allow the selection of
///   the correct residuals methods depending on the parameters type.
template <typename index_container_t, typename measured_t, typename residuals_t>
inline void calculateResiduals(FreeIndices size,
                               const index_container_t& indices,
                               const FreeVector& reference,
                               const Eigen::MatrixBase<measured_t>& measured,
                               Eigen::MatrixBase<residuals_t>& residuals) {
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(measured_t);
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(residuals_t);
  assert((size <= eFreeSize) and "Measured subspace is too large");
  assert((size <= measured.size()) and "Inconsistent measured size");
  assert((size <= residuals.size()) and "Inconsistent residuals size");

  for (size_t i = 0; i < size; ++i) {
    // all free parameters are unrestricted. no need to call parameter traits
    residuals[i] = measured[i] - reference[indices[i]];
  }
}

}  // namespace detail
}  // namespace Acts
