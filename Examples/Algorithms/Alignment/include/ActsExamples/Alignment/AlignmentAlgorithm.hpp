// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/MagneticField/MagneticField.hpp"

#include <functional>
#include <memory>
#include <vector>

namespace ActsExamples {

class AlignmentAlgorithm final : public BareAlgorithm {
 public:
  using AlignmentResult = Acts::Result<ActsAlignment::AlignmentResult>;
  /// Alignment function that takes sets of input measurements, initial
  /// trackstate and alignment options and returns some alignment-specific
  /// result.
  using TrackFitterOptions =
      Acts::KalmanFitterOptions<MeasurementCalibrator, Acts::VoidOutlierFinder>;

  /// Alignment function that takes the above parameters and runs alignment
  /// @note This is separated into a virtual interface to keep compilation units
  /// small
  class AlignmentFunction {
   public:
    virtual ~AlignmentFunction() = default;
    virtual AlignmentResult operator()(
        const std::vector<std::vector<IndexSourceLink>>&,
        const TrackParametersContainer&,
        const ActsAlignment::AlignmentOptions<TrackFitterOptions>&) const = 0;
  };

  /// Create the alignment function implementation.
  ///
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyways.
  static std::shared_ptr<AlignmentFunction> makeAlignmentFunction(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField);

  struct Config {
    /// Input measurements collection.
    std::string inputMeasurements;
    /// Input source links collection.
    std::string inputSourceLinks;
    /// Input proto tracks collection, i.e. groups of hit indices.
    std::string inputProtoTracks;
    /// Input initial track parameter estimates for for each proto track.
    std::string inputInitialTrackParameters;
    /// Output aligned parameters collection.
    std::string outputAlignmentParameters;
    /// Type erased fitter function.
    std::shared_ptr<AlignmentFunction> align;
    /// The alignd transform updater
    ActsAlignment::AlignedTransformUpdater alignedTransformUpdater;
    /// The surfaces (with detector elements) to be aligned
    std::vector<Acts::DetectorElementBase*> alignedDetElements;
    /// The alignment mask at each iteration
    std::map<unsigned int, std::bitset<6>> iterationState;
    /// Cutoff value for average chi2/ndf
    double chi2ONdfCutOff = 0.10;
    /// Cutoff value for delta of average chi2/ndf within a couple of iterations
    std::pair<size_t, double> deltaChi2ONdfCutOff = {10, 0.00001};
    /// Maximum number of iterations
    size_t maxNumIterations = 100;
    /// Number of tracks to be used for alignment
    int maxNumTracks = -1;
  };

  /// Constructor of the alignment algorithm
  ///
  /// @param cfg is the config struct to configure the algorihtm
  /// @param level is the logging level
  AlignmentAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Framework execute method of the alignment algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algporithm flow
  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& ctx) const final override;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
