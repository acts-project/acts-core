// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GainMatrixUpdater.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"

#include <algorithm>
#include <cstddef>
#include <utility>

#include <Eigen/src/Core/MatrixBase.h>

namespace Acts {

std::tuple<double, std::error_code> GainMatrixUpdater::visitMeasurement(
    InternalTrackState trackState, Direction direction,
    const Logger& logger) const {
  // default-constructed error represents success, i.e. an invalid error code
  std::error_code error;
  double chi2 = 0;

  visit_measurement(trackState.calibratedSize, [&](auto N) -> bool {
    constexpr std::size_t kMeasurementSize = decltype(N)::value;
    using ParametersVector = ActsVector<kMeasurementSize>;
    using CovarianceMatrix = ActsSquareMatrix<kMeasurementSize>;

    typename TrackStateTraits<kMeasurementSize, true>::Calibrated calibrated{
        trackState.calibrated};
    typename TrackStateTraits<kMeasurementSize, true>::CalibratedCovariance
        calibratedCovariance{trackState.calibratedCovariance};

    ACTS_VERBOSE("Measurement dimension: " << kMeasurementSize);
    ACTS_VERBOSE("Calibrated measurement: " << calibrated.transpose());
    ACTS_VERBOSE("Calibrated measurement covariance:\n"
                 << calibratedCovariance);

    const auto H = trackState.projector
                       .template topLeftCorner<kMeasurementSize, eBoundSize>()
                       .eval();

    ACTS_VERBOSE("Measurement projector H:\n" << H);

    const auto K = (trackState.predictedCovariance * H.transpose() *
                    (H * trackState.predictedCovariance * H.transpose() +
                     calibratedCovariance)
                        .inverse())
                       .eval();

    ACTS_VERBOSE("Gain Matrix K:\n" << K);

    if (K.hasNaN()) {
      // set to error abort execution
      error = (direction == Direction::Forward)
                  ? KalmanFitterError::ForwardUpdateFailed
                  : KalmanFitterError::BackwardUpdateFailed;
      return false;
    }

    trackState.filtered =
        trackState.predicted + K * (calibrated - H * trackState.predicted);
    trackState.filteredCovariance = (BoundSquareMatrix::Identity() - K * H) *
                                    trackState.predictedCovariance;
    ACTS_VERBOSE("Filtered parameters: " << trackState.filtered.transpose());
    ACTS_VERBOSE("Filtered covariance:\n" << trackState.filteredCovariance);

    ParametersVector residual;
    residual = calibrated - H * trackState.filtered;
    ACTS_VERBOSE("Residual: " << residual.transpose());

    CovarianceMatrix m =
        ((CovarianceMatrix::Identity() - H * K) * calibratedCovariance).eval();

    chi2 = (residual.transpose() * m.inverse() * residual).value();

    ACTS_VERBOSE("Chi2: " << chi2);
    return true;
  });

  return {chi2, error};
}

}  // namespace Acts
