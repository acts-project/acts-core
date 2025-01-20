// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParametersConcept.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"

namespace Acts {

/// Curvilinear track parameters for a single track.
///
/// @tparam particle_hypothesis_t Helper type to interpret the particle charge/momentum
///
/// This is intended as a user-facing data class that adds additional accessors
/// and charge/momentum interpretation on top of the pure parameters vector. All
/// parameters and their corresponding covariance matrix are stored in
/// curvilinear parametrization.
///
/// @see GenericBoundTrackParameters
template <typename particle_hypothesis_t>
class GenericCurvilinearTrackParameters
    : public GenericBoundTrackParameters<particle_hypothesis_t> {
  using Base = GenericBoundTrackParameters<particle_hypothesis_t>;

 public:
  using ParametersVector = BoundVector;
  using CovarianceMatrix = BoundSquareMatrix;
  using ParticleHypothesis = particle_hypothesis_t;

  /// Construct from four-position, direction, and qOverP.
  ///
  /// @param pos4 Track position/time four-vector
  /// @param dir Track direction three-vector; normalization is ignored.
  /// @param qOverP Charge over momentum
  /// @param cov Curvilinear bound parameters covariance matrix
  /// @param particleHypothesis Particle hypothesis
  GenericCurvilinearTrackParameters(const Vector4& pos4, const Vector3& dir,
                                    double qOverP,
                                    std::optional<CovarianceMatrix> cov,
                                    ParticleHypothesis particleHypothesis)
      : Base(CurvilinearSurface(pos4.segment<3>(ePos0), dir).surface(),
             transformFreeToCurvilinearParameters(pos4[eTime], dir, qOverP),
             std::move(cov), std::move(particleHypothesis)) {}

  /// Construct from four-position, angles, and qOverP.
  ///
  /// @param pos4 Track position/time four-vector
  /// @param phi Transverse track direction angle
  /// @param theta Longitudinal track direction angle
  /// @param qOverP Charge over momentum
  /// @param cov Curvilinear bound parameters covariance matrix
  /// @param particleHypothesis Particle hypothesis
  GenericCurvilinearTrackParameters(const Vector4& pos4, double phi,
                                    double theta, double qOverP,
                                    std::optional<CovarianceMatrix> cov,
                                    ParticleHypothesis particleHypothesis)
      : Base(CurvilinearSurface(pos4.segment<3>(ePos0),
                                makeDirectionFromPhiTheta(phi, theta))
                 .surface(),
             transformFreeToCurvilinearParameters(pos4[eTime], phi, theta,
                                                  qOverP),
             std::move(cov), std::move(particleHypothesis)) {}

  /// Converts a bound track parameter with a different hypothesis.
  template <typename other_particle_hypothesis_t>
  explicit GenericCurvilinearTrackParameters(
      const GenericCurvilinearTrackParameters<other_particle_hypothesis_t>&
          other)
      : GenericCurvilinearTrackParameters(other.fourPosition(),
                                          other.particleHypothesis(),
                                          other.covariance()) {}

  /// Converts an unknown bound track parameter.
  template <typename other_track_parameter_t>
  static GenericCurvilinearTrackParameters create(
      const other_track_parameter_t& other) {
    static_assert(BoundTrackParametersConcept<other_track_parameter_t>);

    return GenericCurvilinearTrackParameters(
        other.fourPosition(), other.particleHypothesis(), other.covariance());
  }

  /// Parameters are not default constructible due to the charge type.
  GenericCurvilinearTrackParameters() = delete;
  GenericCurvilinearTrackParameters(const GenericCurvilinearTrackParameters&) =
      default;
  GenericCurvilinearTrackParameters(GenericCurvilinearTrackParameters&&) =
      default;
  ~GenericCurvilinearTrackParameters() = default;
  GenericCurvilinearTrackParameters& operator=(
      const GenericCurvilinearTrackParameters&) = default;
  GenericCurvilinearTrackParameters& operator=(
      GenericCurvilinearTrackParameters&&) = default;

  using GenericBoundTrackParameters<ParticleHypothesis>::fourPosition;
  using GenericBoundTrackParameters<ParticleHypothesis>::position;

  /// Space-time position four-vector.
  Vector4 fourPosition() const {
    return GenericBoundTrackParameters<ParticleHypothesis>::fourPosition({});
  }
  /// Spatial position three-vector.
  Vector3 position() const {
    return GenericBoundTrackParameters<ParticleHypothesis>::position({});
  }

  /// Reflect the parameters.
  /// @return Reflected parameters.
  GenericCurvilinearTrackParameters<ParticleHypothesis> reflect() const {
    GenericCurvilinearTrackParameters<ParticleHypothesis> reflected = *this;
    reflected.reflectInPlace();
    return reflected;
  }
};

}  // namespace Acts
