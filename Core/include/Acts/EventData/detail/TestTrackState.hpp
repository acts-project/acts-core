// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/detail/GenerateParameters.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

#include <random>

namespace Acts::detail::Test {

struct TestTrackState {
  std::shared_ptr<Surface> surface;
  TestSourceLink sourceLink;
  BoundTrackParameters predicted;
  BoundTrackParameters filtered;
  BoundTrackParameters smoothed;
  BoundMatrix jacobian;
  double chi2;
  double pathLength;

  // Generate a random TestTrackState.
  //
  // @param rng Random number generator
  // @param std::size_t nMeasurement either 1 or 2
  template <typename rng_t>
  TestTrackState(rng_t& rng, std::size_t measdim)
      : surface(
            CurvilinearSurface(Vector3::Zero(), Vector3::UnitZ()).surface()),
        // set bogus parameters first since they are not default-constructible
        predicted(surface, someBoundParametersA(), std::nullopt,
                  ParticleHypothesis::pion()),
        filtered(surface, someBoundParametersA(), std::nullopt,
                 ParticleHypothesis::pion()),
        smoothed(surface, someBoundParametersA(), std::nullopt,
                 ParticleHypothesis::pion()),
        jacobian(BoundMatrix::Identity()),
        chi2(std::chi_squared_distribution<double>(measdim)(rng)),
        pathLength(std::uniform_real_distribution<double>(
            1 * Acts::UnitConstants::mm, 10 * Acts::UnitConstants::mm)(rng)) {
    // set a random geometry identifier to uniquely identify each surface
    Acts::GeometryIdentifier geoId{
        std::uniform_int_distribution<GeometryIdentifier::Value>()(rng)};
    surface->assignGeometryId(geoId);

    // create source link w/ inline 1d or 2d measurement data
    if (measdim == 1u) {
      auto [par, cov] = generateParametersCovariance<double, 1u>(rng);
      sourceLink = TestSourceLink(eBoundLoc0, par[0], cov(0, 0), geoId);
    } else if (measdim == 2u) {
      auto [par, cov] = generateParametersCovariance<double, 2u>(rng);
      sourceLink = TestSourceLink(eBoundLoc1, eBoundQOverP, par, cov, geoId);
    } else {
      throw std::runtime_error("invalid number of measurement dimensions");
    }

    // create track parameters
    auto [trkPar, trkCov] = generateBoundParametersCovariance(rng, {});
    // trkPar[eBoundPhi] = 45_degree;
    // trkPar[eBoundTheta] = 90_degree;
    // trkPar[eBoundQOverP] = 5.;
    // predicted
    predicted = BoundTrackParameters(surface, trkPar, trkCov,
                                     ParticleHypothesis::pion());
    // filtered, modified q/p, reduced covariance
    // trkPar[eBoundQOverP] = 10.;
    filtered = BoundTrackParameters(surface, trkPar, 0.75 * trkCov,
                                    ParticleHypothesis::pion());
    // smoothed, modified q/p, further reduced covariance
    // trkPar[eBoundQOverP] = 15.;
    smoothed = BoundTrackParameters(surface, trkPar, 0.5 * trkCov,
                                    ParticleHypothesis::pion());

    // propagation jacobian is identity + corrections
    for (Eigen::Index c = 0; c < jacobian.cols(); ++c) {
      for (Eigen::Index r = 0; r < jacobian.rows(); ++r) {
        jacobian(c, r) +=
            std::uniform_real_distribution<double>(-0.125, 0.125)(rng);
      }
    }
  }
};

// Fill a TrackStateProxy with values from a TestTrackState.
//
// @param[in] pc TestTrackState with the input values
// @param[in] mask Specifies which components are used/filled
// @param[out] ts TrackStateProxy which is filled
// @param [in] measdim Dimension of the measurement
template <typename trajectory_t, typename track_state_t>
void fillTrackState(const TestTrackState& pc, TrackStatePropMask mask,
                    track_state_t& ts) {
  // always set the reference surface
  ts.setReferenceSurface(pc.predicted.referenceSurface().getSharedPtr());

  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Predicted)) {
    ts.predicted() = pc.predicted.parameters();
    assert(pc.predicted.covariance().has_value());
    ts.predictedCovariance() = *(pc.predicted.covariance());
  }
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Filtered)) {
    ts.filtered() = pc.filtered.parameters();
    assert(pc.filtered.covariance().has_value());
    ts.filteredCovariance() = *(pc.filtered.covariance());
  }
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Smoothed)) {
    ts.smoothed() = pc.smoothed.parameters();
    assert(pc.smoothed.covariance().has_value());
    ts.smoothedCovariance() = *(pc.smoothed.covariance());
  }
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Jacobian)) {
    ts.jacobian() = pc.jacobian;
  }
  ts.chi2() = pc.chi2;
  ts.pathLength() = pc.pathLength;
  // source link defines the uncalibrated measurement
  // create calibrated measurements from source link
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Calibrated)) {
    testSourceLinkCalibrator<trajectory_t>(Acts::GeometryContext{},
                                           Acts::CalibrationContext{},
                                           SourceLink{pc.sourceLink}, ts);
    assert(ts.hasUncalibratedSourceLink());
  }
}

}  // namespace Acts::detail::Test
