// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/KalmanGlobalCovariance.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsAlignment/Kernel/AlignmentError.hpp"
#include "ActsAlignment/Kernel/detail/AlignmentEngine.hpp"

#include <limits>
#include <map>
#include <queue>
#include <vector>

namespace ActsAlignment {
using AlignedTransformUpdater =
    std::function<bool(Acts::DetectorElementBase*, const Acts::GeometryContext&,
                       const Acts::Transform3&)>;
///
/// @brief Options for align() call
///
/// @tparam fit_options_t The fit options type
///
template <typename fit_options_t>
struct AlignmentOptions {
  /// Deleted default constructor
  AlignmentOptions() = delete;

  /// AlignmentOptions
  ///
  /// @param fOptions The fit options
  /// @param logger_ The logger wrapper
  /// @param aTransformUpdater The updater to update aligned transform
  /// @param aDetElements The alignable detector elements
  /// @param chi2CufOff The alignment chi2 tolerance
  /// @param deltaChi2CutOff The change of chi2 within a few iterations
  /// @param maxIters The alignment maximum iterations

  AlignmentOptions(
      const fit_options_t& fOptions,
      const AlignedTransformUpdater& aTransformUpdater,
      Acts::LoggerWrapper logger_,
      const std::vector<Acts::DetectorElementBase*>& aDetElements = {},
      double chi2CutOff = 0.5,
      const std::pair<size_t, double>& deltaChi2CutOff = {5, 0.01},
      size_t maxIters = 5,
      const std::map<unsigned int, AlignmentMask>& iterState = {})
      : fitOptions(fOptions),
        alignedTransformUpdater(aTransformUpdater),
        alignedDetElements(aDetElements),
        averageChi2ONdfCutOff(chi2CutOff),
        deltaAverageChi2ONdfCutOff(deltaChi2CutOff),
        maxIterations(maxIters),
        iterationState(iterState),
        logger(logger_) {}

  // The fit options
  fit_options_t fitOptions;

  /// The updater to the aligned transform
  AlignedTransformUpdater alignedTransformUpdater = nullptr;

  // The detector elements to be aligned
  std::vector<Acts::DetectorElementBase*> alignedDetElements;

  // The alignment tolerance to determine if the alignment is coveraged
  double averageChi2ONdfCutOff = 0.5;

  // The delta of average chi2/ndf within a couple of iterations to determine if
  // alignment is converged
  std::pair<size_t, double> deltaAverageChi2ONdfCutOff = {5, 0.01};

  // The maximum number of iterations to run alignment
  size_t maxIterations = 5;

  // The alignment mask for different iterations
  std::map<unsigned int, AlignmentMask> iterationState;

  /// Logger
  Acts::LoggerWrapper logger;
};

/// @brief Alignment result struct
///
struct AlignmentResult {
  // The change of alignment parameters
  Acts::ActsDynamicVector deltaAlignmentParameters;

  // The aligned parameters for detector elements
  std::unordered_map<Acts::DetectorElementBase*, Acts::Transform3>
      alignedParameters;

  // The covariance of alignment parameters
  Acts::ActsDynamicMatrix alignmentCovariance;

  // The avarage chi2/ndf (ndf is the measurement dim)
  double averageChi2ONdf = std::numeric_limits<double>::max();

  // The delta chi2
  double deltaChi2 = std::numeric_limits<double>::max();

  // The chi2
  double chi2 = 0;

  // The measurement dimension from all tracks
  size_t measurementDim = 0;

  // The alignment degree of freedom
  size_t alignmentDof = 0;

  // The number of tracks used for alignment
  size_t numTracks = 0;

  // The indexed alignable surfaces
  std::unordered_map<const Acts::Surface*, size_t> idxedAlignSurfaces;

  Acts::Result<void> result{Acts::Result<void>::success()};
};

/// @brief KalmanFitter-based alignment implementation
///
/// @tparam fitter_t Type of the fitter class
template <typename fitter_t>
struct Alignment {
  /// Default constructor is deleted
  Alignment() = delete;

  /// Constructor from arguments
  Alignment(fitter_t fitter) : m_fitter(std::move(fitter)) {}

  /// @brief evaluate alignment state for a single track
  ///
  /// @tparam source_link_t Source link type identifying uncalibrated input
  /// measurements.
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam fit_options_t The fit options type
  ///
  /// @param gctx The current geometry context object
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param fitOptions The fit Options steering the fit
  /// @param idxedAlignSurfaces The idxed surfaces to be aligned
  /// @param alignMask The alignment mask (same for all detector element for the
  /// moment)
  ///
  /// @result The alignment state for a single track
  template <typename source_link_t, typename start_parameters_t,
            typename fit_options_t>
  Acts::Result<detail::TrackAlignmentState> evaluateTrackAlignmentState(
      const Acts::GeometryContext& gctx,
      const std::vector<source_link_t>& sourcelinks,
      const start_parameters_t& sParameters, const fit_options_t& fitOptions,
      const std::unordered_map<const Acts::Surface*, size_t>&
          idxedAlignSurfaces,
      const AlignmentMask& alignMask,
      Acts::LoggerWrapper logger = Acts::getDummyLogger()) const;

  /// @brief calculate the alignment parameters delta
  ///
  /// @tparam trajectory_container_t The trajectories container type
  /// @tparam start_parameters_t The initial parameters container type
  /// @tparam fit_options_t The fit options type
  ///
  /// @param trajectoryCollection The collection of trajectories as input of
  /// fitting
  /// @param startParametersCollection The collection of starting parameters as
  /// input of fitting
  /// @param fitOptions The fit Options steering the fit
  /// @param alignResult [in, out] The aligned result
  /// @param alignMask The alignment mask (same for all measurements now)
  template <typename trajectory_container_t,
            typename start_parameters_container_t, typename fit_options_t>
  void calculateAlignmentParameters(
      const trajectory_container_t& trajectoryCollection,
      const start_parameters_container_t& startParametersCollection,
      const fit_options_t& fitOptions, AlignmentResult& alignResult,
      const AlignmentMask& alignMask = AlignmentMask::All,
      Acts::LoggerWrapper logger = Acts::getDummyLogger()) const;

  /// @brief update the detector element alignment parameters
  ///
  /// @param gctx The geometry context
  /// @param alignedDetElements The detector elements to be aligned
  /// @param alignedTransformUpdater The updater for updating the aligned
  /// @param alignResult [in, out] The aligned result
  Acts::Result<void> updateAlignmentParameters(
      const Acts::GeometryContext& gctx,
      const std::vector<Acts::DetectorElementBase*>& alignedDetElements,
      const AlignedTransformUpdater& alignedTransformUpdater,
      AlignmentResult& alignResult,
      Acts::LoggerWrapper logger = Acts::getDummyLogger()) const;

  /// @brief Alignment implementation
  ///
  /// @tparam trajectory_container_t The trajectories container type
  /// @tparam start_parameters_t The initial parameters container type
  /// @tparam fit_options_t The fit options type
  ///
  /// @param trajectoryCollection The collection of trajectories as input of
  /// fitting
  /// @param startParametersCollection The collection of starting parameters as
  /// input of fitting
  /// @param alignOptions The alignment options
  ///
  /// @result The alignment result
  template <typename trajectory_container_t,
            typename start_parameters_container_t, typename fit_options_t>
  Acts::Result<AlignmentResult> align(
      const trajectory_container_t& trajectoryCollection,
      const start_parameters_container_t& startParametersCollection,
      const AlignmentOptions<fit_options_t>& alignOptions) const;

 private:
  // The fitter
  fitter_t m_fitter;
};
}  // namespace ActsAlignment

#include "ActsAlignment/Kernel/Alignment.ipp"
