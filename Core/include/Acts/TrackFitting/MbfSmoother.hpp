// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cassert>
#include <cstddef>
#include <system_error>

namespace Acts {

/// Kalman trajectory smoother based on the Modified Bryson–Frazier smoother.
///
/// This implements not a single smoothing step, but the full backwards
/// smoothing procedure for a filtered, forward trajectory using the stored
/// linearization.
class MbfSmoother {
 public:
  struct InternalTrackState {
    using Projector =
        typename TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                                  false>::Projector;
    using Jacobian =
        typename TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                                  false>::Covariance;
    using Parameters =
        typename TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                                  false>::Parameters;
    using Covariance =
        typename TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                                  false>::Covariance;

    // This is used to build a covariance matrix view in the .cpp file
    unsigned int calibratedSize{0};
    const double* calibrated{nullptr};
    const double* calibratedCovariance{nullptr};
    Projector projector;

    Jacobian jacobian{nullptr};

    Parameters predicted{nullptr};
    Covariance predictedCovariance{nullptr};
    Parameters filtered{nullptr};
    Covariance filteredCovariance{nullptr};
    Parameters smoothed{nullptr};
    Covariance smoothedCovariance{nullptr};

    InternalTrackState() = default;

    template <typename TrackStateProxy>
    explicit InternalTrackState(TrackStateProxy ts)
        : jacobian(ts.jacobian()),
          predicted(ts.predicted()),
          predictedCovariance(ts.predictedCovariance()),
          filtered(ts.filtered()),
          filteredCovariance(ts.filteredCovariance()),
          smoothed(ts.smoothed()),
          smoothedCovariance(ts.smoothedCovariance()) {
      if (ts.typeFlags().test(TrackStateFlag::MeasurementFlag)) {
        calibratedSize = ts.calibratedSize();
        // Note that we pass raw pointers here which are used in the correct
        // shape later
        calibrated = {ts.effectiveCalibrated().data()};
        calibratedCovariance = {ts.effectiveCalibratedCovariance().data()};
        projector = {ts.projector()};
      }
    }
  };

  /// Run the Kalman smoothing for one trajectory.
  ///
  /// @param[in] gctx The geometry context to be used
  /// @param[in,out] trajectory The trajectory to be smoothed
  /// @param[in] entryIndex The index of state to start the smoothing
  /// @param[in] logger Where to write logging information to
  template <typename traj_t>
  Result<void> operator()(const GeometryContext& gctx, traj_t& trajectory,
                          std::size_t entryIndex,
                          const Logger& logger = getDummyLogger()) const {
    (void)gctx;
    (void)logger;

    using TrackStateProxy = typename traj_t::TrackStateProxy;

    TrackStateProxy start_ts = trajectory.getTrackState(entryIndex);

    BoundMatrix big_lambda_hat = BoundMatrix::Zero();
    BoundVector small_lambda_hat = BoundVector::Zero();

    trajectory.applyBackwards(start_ts.index(), [&](TrackStateProxy ts) {
      // ensure the track state has a smoothed component
      ts.addComponents(TrackStatePropMask::Smoothed);

      InternalTrackState internalTrackState(ts);

      calculateSmoothed(internalTrackState, big_lambda_hat, small_lambda_hat);

      if (!ts.hasPrevious()) {
        return;
      }

      if (ts.typeFlags().test(TrackStateFlag::MeasurementFlag)) {
        visitMeasurement(internalTrackState, big_lambda_hat, small_lambda_hat);
      } else {
        visitNonMeasurement(internalTrackState, big_lambda_hat,
                            small_lambda_hat);
      }
    });

    return Result<void>::success();
  }

 private:
  void calculateSmoothed(InternalTrackState& ts,
                         const BoundMatrix& big_lambda_hat,
                         const BoundVector& small_lambda_hat) const;
  void visitNonMeasurement(const InternalTrackState& ts,
                           BoundMatrix& big_lambda_hat,
                           BoundVector& small_lambda_hat) const;
  void visitMeasurement(const InternalTrackState& ts,
                        BoundMatrix& big_lambda_hat,
                        BoundVector& small_lambda_hat) const;
};

}  // namespace Acts