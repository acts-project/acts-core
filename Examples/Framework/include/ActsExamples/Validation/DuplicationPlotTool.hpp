// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Helpers.hpp"

#include <cstddef>
#include <map>
#include <memory>
#include <string>

class TEfficiency;
class TProfile;
namespace ActsFatras {
class Particle;
}  // namespace ActsFatras

namespace ActsExamples {

// Tools to make duplication rate and duplication number plots to show tracking
// duplication.
//
// The duplication is investigated for those truth-matched reco tracks. If there
// are a few reco tracks matched to the same truth particle, the reco track with
// the highest matching probability is tagges as 'real' and the others are
// 'duplicated'.
class DuplicationPlotTool {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, PlotHelpers::Binning> varBinning = {
        {"Eta", PlotHelpers::Binning("#eta", 40, -4, 4)},
        {"Phi", PlotHelpers::Binning("#phi", 100, -3.15, 3.15)},
        {"Pt", PlotHelpers::Binning("p_{T} [GeV/c]", 40, 0, 100)},
        {"Num", PlotHelpers::Binning("N", 30, -0.5, 29.5)}};
  };

  /// @brief Nested Cache struct
  struct DuplicationPlotCache {
    TProfile* nDuplicated_vs_pT;   ///< Number of duplicated tracks vs pT
    TProfile* nDuplicated_vs_eta;  ///< Number of duplicated tracks vs eta
    TProfile* nDuplicated_vs_phi;  ///< Number of duplicated tracks vs phi
    TProfile2D*
        nDuplicated_vs_pT_eta;  ///< Number of duplicated tracks vs pT and eta
    TEfficiency* duplicationRate_vs_pT;   ///< Tracking duplication rate vs pT
    TEfficiency* duplicationRate_vs_eta;  ///< Tracking duplication rate vs eta
    TEfficiency* duplicationRate_vs_phi;  ///< Tracking duplication rate vs phi
    TEfficiency*
        duplicationRate_vs_pT_eta;  ///< 2d tracking duplication rate vs pT
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  DuplicationPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief book the duplication plots
  ///
  /// @param duplicationPlotCache the cache for duplication plots
  void book(DuplicationPlotCache& duplicationPlotCache) const;

  /// @brief fill duplication rate w.r.t. fitted track parameters
  ///
  /// @param duplicationPlotCache cache object for duplication plots
  /// @param fittedParameters fitted track parameters of this track
  /// @param status the (truth-matched) reconstructed track is duplicated or not
  void fill(DuplicationPlotCache& duplicationPlotCache,
            const Acts::BoundTrackParameters& fittedParameters,
            bool status) const;

  /// @brief fill number of duplicated tracks for a truth particle seed
  ///
  /// @param duplicationPlotCache cache object for duplication plots
  /// @param truthParticle the truth Particle
  /// @param nDuplicatedTracks the number of duplicated tracks
  void fill(DuplicationPlotCache& duplicationPlotCache,
            const SimParticleState& truthParticle,
            std::size_t nDuplicatedTracks) const;

  /// @brief write the duplication plots to file
  ///
  /// @param duplicationPlotCache cache object for duplication plots
  void write(const DuplicationPlotCache& duplicationPlotCache) const;

  /// @brief delete the duplication plots
  ///
  /// @param duplicationPlotCache cache object for duplication plots
  void clear(DuplicationPlotCache& duplicationPlotCache) const;

 private:
  Config m_cfg;                                  ///< The Config class
  std::unique_ptr<const Acts::Logger> m_logger;  ///< The logging instance

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
