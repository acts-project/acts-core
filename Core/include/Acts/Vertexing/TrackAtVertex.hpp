// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"

#include <functional>

namespace Acts {

/// @class TrackAtVertex
///
/// @brief Defines a track at vertex object
///
/// @tparam input_track_t Track object type

struct FittedMomentum {
  FittedMomentum(Vector3 mom, std::optional<ActsSymMatrix<3>> cov)
      : momentum(mom), covariance(cov) {}

  FittedMomentum() : momentum(Vector3::Zero()), covariance(std::nullopt) {}

  Vector3 momentum;
  // TODO: do we even need the optional?
  std::optional<ActsSymMatrix<3>> covariance;
};

template <typename input_track_t>

struct TrackAtVertex {
  /// Deleted default constructor
  TrackAtVertex() = delete;

  /// @brief Parameterized constructor
  ///
  /// @param chi2perTrack Chi2 of track
  /// @param paramsAtVertex Fitted perigee parameter
  /// @param originalTrack Original perigee parameter
  TrackAtVertex(double chi2perTrack, const BoundTrackParameters& paramsAtVertex,
                const input_track_t* originalTrack)
      : fittedParams(paramsAtVertex),
        originalParams(originalTrack),
        chi2Track(chi2perTrack) {}

  TrackAtVertex(double chi2perTrack, const BoundTrackParameters& paramsAtVertex,
                const FittedMomentum& fittedMom, const input_track_t* originalTrack)
      : fittedParams(paramsAtVertex),
        fittedMomentum(fittedMom),
        originalParams(originalTrack),
        chi2Track(chi2perTrack) {}

  /// @brief Constructor with default chi2
  ///
  /// @param paramsAtVertex Fitted perigee parameter
  /// @param originalTrack Original perigee parameter
  TrackAtVertex(const BoundTrackParameters& paramsAtVertex,
                const input_track_t* originalTrack)
      : fittedParams(paramsAtVertex), originalParams(originalTrack) {}

  TrackAtVertex(const BoundTrackParameters& paramsAtVertex, const FittedMomentum& fittedMom,
                const input_track_t* originalTrack)
      : fittedParams(paramsAtVertex), fittedMomentum(fittedMom), originalParams(originalTrack) {}

  /// Fitted perigee
  BoundTrackParameters fittedParams;

  /// Momentum after vertex fit
  FittedMomentum fittedMomentum;

  /// Original input parameters
  const input_track_t* originalParams;

  /// Chi2 of track
  double chi2Track = 0;

  /// Number degrees of freedom
  /// Note: Can be different from integer value
  /// since annealing can result in effective
  /// non-interger values
  double ndf = 0;

  /// Value of the compatibility of the track to the actual vertex, based
  /// on the estimation of the 3d distance between the track and the vertex
  double vertexCompatibility = 0;

  /// Weight of track in fit
  double trackWeight = 1;

  /// The linearized state of the track at vertex
  LinearizedTrack linearizedState;

  /// Is already linearized
  bool isLinearized = false;
};

}  // namespace Acts
