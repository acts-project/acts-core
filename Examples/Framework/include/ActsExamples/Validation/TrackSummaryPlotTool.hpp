// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Utilities/Helpers.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <map>
#include <memory>
#include <string>

namespace ActsExamples {

// Tools to make track info plots to show tracking track info.
class TrackSummaryPlotTool {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, PlotHelpers::Binning> varBinning = {
        {"Eta", PlotHelpers::Binning("#eta", 40, -4, 4)},
        {"Phi", PlotHelpers::Binning("#phi", 100, -3.15, 3.15)},
        {"Pt", PlotHelpers::Binning("pT [GeV/c]", 40, 0, 100)},
        {"Num", PlotHelpers::Binning("N", 30, -0.5, 29.5)}};
  };

  /// @brief Nested Cache struct
  struct TrackSummaryPlotCache {
    TProfile* nStates_vs_eta;  ///< Number of total states vs eta
    TProfile*
        nMeasurements_vs_eta;    ///< Number of non-outlier measurements vs eta
    TProfile* nHoles_vs_eta;     ///< Number of holes vs eta
    TProfile* nOutliers_vs_eta;  ///< Number of outliers vs eta
    TProfile* nSharedHits_vs_eta;  ///< Number of Shared Hits vs eta
    TProfile* nStates_vs_pt;       ///< Number of total states vs pt
    TProfile*
        nMeasurements_vs_pt;      ///< Number of non-outlier measurements vs pt
    TProfile* nHoles_vs_pt;       ///< Number of holes vs pt
    TProfile* nOutliers_vs_pt;    ///< Number of outliers vs pt
    TProfile* nSharedHits_vs_pt;  ///< Number of Shared Hits vs pt
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  TrackSummaryPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief book the track info plots
  ///
  /// @param trackSummaryPlotCache the cache for track info plots
  void book(TrackSummaryPlotCache& trackSummaryPlotCache) const;

  /// @brief fill reco track info w.r.t. fitted track parameters
  ///
  /// @param trackSummaryPlotCache cache object for track info plots
  /// @param fittedParameters fitted track parameters of this track
  /// @param nStates number of track states
  /// @param nMeasurements number of measurements
  /// @param nOutliers number of outliers
  /// @param nHoles number of holes
  void fill(TrackSummaryPlotCache& trackSummaryPlotCache,
            const Acts::BoundTrackParameters& fittedParameters, size_t nStates,
            size_t nMeasurments, size_t Outliers, size_t nHoles,
            size_t nSharedHits) const;

  /// @brief write the track info plots to file
  ///
  /// @param trackSummaryPlotCache cache object for track info plots
  void write(const TrackSummaryPlotCache& trackSummaryPlotCache) const;

  /// @brief delete the track info plots
  ///
  /// @param trackSummaryPlotCache cache object for track info plots
  void clear(TrackSummaryPlotCache& trackSummaryPlotCache) const;

 private:
  Config m_cfg;                                  ///< The Config class
  std::unique_ptr<const Acts::Logger> m_logger;  ///< The logging instance

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
