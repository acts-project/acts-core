// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/MbfSmoother.hpp"

namespace Acts {

void MbfSmoother::calculateSmoothed(InternalTrackState& ts,
                                    const BoundMatrix& big_lambda_hat,
                                    const BoundVector& small_lambda_hat) const {
  ts.smoothedCovariance = ts.filteredCovariance - ts.filteredCovariance *
                                                      big_lambda_hat *
                                                      ts.filteredCovariance;
  ts.smoothed = ts.filtered - ts.filteredCovariance * small_lambda_hat;
}

void MbfSmoother::visitNonMeasurement(const InternalTrackState& ts,
                                      BoundMatrix& big_lambda_hat,
                                      BoundVector& small_lambda_hat) const {
  const auto F = ts.jacobian;

  big_lambda_hat = F.transpose() * big_lambda_hat * F;
  small_lambda_hat = F.transpose() * small_lambda_hat;
}

void MbfSmoother::visitMeasurement(const InternalTrackState& ts,
                                   BoundMatrix& big_lambda_hat,
                                   BoundVector& small_lambda_hat) const {
  visit_measurement(ts.calibratedSize, [&](auto N) -> void {
    constexpr std::size_t kMeasurementSize = decltype(N)::value;

    typename TrackStateTraits<kMeasurementSize, true>::Calibrated calibrated{
        ts.calibrated};
    typename TrackStateTraits<kMeasurementSize, true>::CalibratedCovariance
        calibratedCovariance{ts.calibratedCovariance};

    const auto H =
        ts.projector.template topLeftCorner<kMeasurementSize, eBoundSize>()
            .eval();

    const auto S =
        (H * ts.predictedCovariance * H.transpose() + calibratedCovariance)
            .eval();
    // TODO Sinv could be cached by the filter step
    const auto Sinv = S.inverse().eval();

    // TODO K could be cached by the filter step
    const auto K = (ts.predictedCovariance * H.transpose() * Sinv).eval();

    const auto C = (BoundMatrix::Identity() - K * H).eval();
    const auto y = (calibrated - H * ts.predicted).eval();

    const auto big_lambda_tilde =
        (H.transpose() * Sinv * H + C.transpose() * big_lambda_hat * C).eval();
    const auto small_lambda_tilde =
        (-H.transpose() * Sinv * y + C.transpose() * small_lambda_hat).eval();

    const auto F = ts.jacobian;

    big_lambda_hat = F.transpose() * big_lambda_tilde * F;
    small_lambda_hat = F.transpose() * small_lambda_tilde;
  });
}

}  // namespace Acts