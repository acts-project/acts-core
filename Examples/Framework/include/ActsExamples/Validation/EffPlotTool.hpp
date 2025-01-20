// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Helpers.hpp"

#include <map>
#include <memory>
#include <string>

class TEfficiency;
namespace ActsFatras {
class Particle;
}  // namespace ActsFatras

namespace ActsExamples {

// Tools to make efficiency plots to show tracking efficiency.
// For the moment, the efficiency is taken as the fraction of successfully
// smoothed track over all tracks
class EffPlotTool {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, PlotHelpers::Binning> varBinning = {
        {"Eta", PlotHelpers::Binning("#eta", 40, -4, 4)},
        {"Phi", PlotHelpers::Binning("#phi", 100, -3.15, 3.15)},
        {"Pt", PlotHelpers::Binning("p_{T} [GeV/c]", 40, 0, 100)},
        {"Z0", PlotHelpers::Binning("z_{0} [mm]", 50, -200, 200)},
        {"DeltaR", PlotHelpers::Binning("#Delta R", 100, 0, 0.3)}};
  };

  /// @brief Nested Cache struct
  struct EffPlotCache {
    TEfficiency* trackEff_vs_pT{nullptr};   ///< Tracking efficiency vs pT
    TEfficiency* trackEff_vs_eta{nullptr};  ///< Tracking efficiency vs eta
    TEfficiency* trackEff_vs_phi{nullptr};  ///< Tracking efficiency vs phi
    TEfficiency* trackEff_vs_z0{nullptr};   ///< Tracking efficiency vs z0
    TEfficiency* trackEff_vs_DeltaR{
        nullptr};  ///< Tracking efficiency vs distance to the closest truth
                   ///< particle
    TEfficiency* track2dEff_vs_pT_eta{
        nullptr};  ///< 2d tracking efficiency vs pT and eta
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  EffPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief book the efficiency plots
  ///
  /// @param effPlotCache the cache for efficiency plots
  void book(EffPlotCache& effPlotCache) const;

  /// @brief fill efficiency plots
  ///
  /// @param effPlotCache cache object for efficiency plots
  /// @param truthParticle the truth Particle
  /// @param deltaR the distance to the closest truth particle
  /// @param status the reconstruction status
  void fill(EffPlotCache& effPlotCache, const SimParticleState& truthParticle,
            double deltaR, bool status) const;

  /// @brief write the efficiency plots to file
  ///
  /// @param effPlotCache cache object for efficiency plots
  void write(const EffPlotCache& effPlotCache) const;

  /// @brief delete the efficiency plots
  ///
  /// @param effPlotCache cache object for efficiency plots
  void clear(EffPlotCache& effPlotCache) const;

 private:
  Config m_cfg;                                  ///< The Config class
  std::unique_ptr<const Acts::Logger> m_logger;  ///< The logging instance

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
