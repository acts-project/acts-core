// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <TProfile.h>

ActsExamples::TrackSummaryPlotTool::TrackSummaryPlotTool(
    const ActsExamples::TrackSummaryPlotTool::Config& cfg,
    Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("TrackSummaryPlotTool", lvl)) {}

void ActsExamples::TrackSummaryPlotTool::book(
    TrackSummaryPlotTool::TrackSummaryPlotCache& trackSummaryPlotCache) const {
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
  PlotHelpers::Binning bNum = m_cfg.varBinning.at("Num");
  ACTS_DEBUG("Initialize the histograms for track info plots");
  // number of track states versus eta
  trackSummaryPlotCache.nStates_vs_eta = PlotHelpers::bookProf(
      "nStates_vs_eta", "Number of total states vs. #eta", bEta, bNum);
  // number of measurements versus eta
  trackSummaryPlotCache.nMeasurements_vs_eta = PlotHelpers::bookProf(
      "nMeasurements_vs_eta", "Number of measurements vs. #eta", bEta, bNum);
  // number of holes versus eta
  trackSummaryPlotCache.nHoles_vs_eta = PlotHelpers::bookProf(
      "nHoles_vs_eta", "Number of holes vs. #eta", bEta, bNum);
  // number of outliers versus eta
  trackSummaryPlotCache.nOutliers_vs_eta = PlotHelpers::bookProf(
      "nOutliers_vs_eta", "Number of outliers vs. #eta", bEta, bNum);
  // number of Shared Hits versus eta
  trackSummaryPlotCache.nSharedHits_vs_eta = PlotHelpers::bookProf(
      "nSharedHits_vs_eta", "Number of Shared Hits vs. #eta", bEta, bNum);
  // number of track states versus pt
  trackSummaryPlotCache.nStates_vs_pt = PlotHelpers::bookProf(
      "nStates_vs_pT", "Number of total states vs. pT", bPt, bNum);
  // number of measurements versus pt
  trackSummaryPlotCache.nMeasurements_vs_pt = PlotHelpers::bookProf(
      "nMeasurements_vs_pT", "Number of measurements vs. pT", bPt, bNum);
  // number of holes versus pt
  trackSummaryPlotCache.nHoles_vs_pt = PlotHelpers::bookProf(
      "nHoles_vs_pT", "Number of holes vs. pT", bPt, bNum);
  // number of outliers versus pt
  trackSummaryPlotCache.nOutliers_vs_pt = PlotHelpers::bookProf(
      "nOutliers_vs_pT", "Number of outliers vs. pT", bPt, bNum);
  // number of Shared Hits versus pt
  trackSummaryPlotCache.nSharedHits_vs_pt = PlotHelpers::bookProf(
      "nSharedHits_vs_pT", "Number of Shared Hits vs. pT", bPt, bNum);
  // number of track states versus eta and pt
  trackSummaryPlotCache.nStates_vs_eta_pt = PlotHelpers::bookProf2D(
      "nStates_vs_eta_pt", "Number of total states vs. #eta and pT", bEta, bPt,
      bNum);
  // number of measurements versus eta and pt
  trackSummaryPlotCache.nMeasurements_vs_eta_pt = PlotHelpers::bookProf2D(
      "nMeasurements_vs_eta_pt", "Number of measurements vs. #eta and pT", bEta,
      bPt, bNum);
  // number of holes versus eta and pt
  trackSummaryPlotCache.nHoles_vs_eta_pt = PlotHelpers::bookProf2D(
      "nHoles_vs_eta_pt", "Number of holes vs. #eta and pT", bEta, bPt, bNum);
  // number of outliers versus eta and pt
  trackSummaryPlotCache.nOutliers_vs_eta_pt = PlotHelpers::bookProf2D(
      "nOutliers_vs_eta_pt", "Number of outliers vs. #eta and pT", bEta, bPt,
      bNum);
  // number of Shared Hits versus eta and pt
  trackSummaryPlotCache.nSharedHits_vs_eta_pt = PlotHelpers::bookProf2D(
      "nSharedHits_vs_eta_pt", "Number of Shared Hits vs. #eta and pT", bEta,
      bPt, bNum);
}

void ActsExamples::TrackSummaryPlotTool::clear(
    TrackSummaryPlotCache& trackSummaryPlotCache) const {
  delete trackSummaryPlotCache.nStates_vs_eta;
  delete trackSummaryPlotCache.nMeasurements_vs_eta;
  delete trackSummaryPlotCache.nOutliers_vs_eta;
  delete trackSummaryPlotCache.nHoles_vs_eta;
  delete trackSummaryPlotCache.nSharedHits_vs_eta;
  delete trackSummaryPlotCache.nStates_vs_pt;
  delete trackSummaryPlotCache.nMeasurements_vs_pt;
  delete trackSummaryPlotCache.nOutliers_vs_pt;
  delete trackSummaryPlotCache.nHoles_vs_pt;
  delete trackSummaryPlotCache.nSharedHits_vs_pt;
  delete trackSummaryPlotCache.nStates_vs_eta_pt;
  delete trackSummaryPlotCache.nMeasurements_vs_eta_pt;
  delete trackSummaryPlotCache.nOutliers_vs_eta_pt;
  delete trackSummaryPlotCache.nHoles_vs_eta_pt;
  delete trackSummaryPlotCache.nSharedHits_vs_eta_pt;
}

void ActsExamples::TrackSummaryPlotTool::write(
    const TrackSummaryPlotTool::TrackSummaryPlotCache& trackSummaryPlotCache)
    const {
  ACTS_DEBUG("Write the plots to output file.");
  trackSummaryPlotCache.nStates_vs_eta->Write();
  trackSummaryPlotCache.nMeasurements_vs_eta->Write();
  trackSummaryPlotCache.nOutliers_vs_eta->Write();
  trackSummaryPlotCache.nHoles_vs_eta->Write();
  trackSummaryPlotCache.nSharedHits_vs_eta->Write();
  trackSummaryPlotCache.nStates_vs_pt->Write();
  trackSummaryPlotCache.nMeasurements_vs_pt->Write();
  trackSummaryPlotCache.nOutliers_vs_pt->Write();
  trackSummaryPlotCache.nHoles_vs_pt->Write();
  trackSummaryPlotCache.nSharedHits_vs_pt->Write();
  trackSummaryPlotCache.nStates_vs_eta_pt->Write();
  trackSummaryPlotCache.nMeasurements_vs_eta_pt->Write();
  trackSummaryPlotCache.nOutliers_vs_eta_pt->Write();
  trackSummaryPlotCache.nHoles_vs_eta_pt->Write();
  trackSummaryPlotCache.nSharedHits_vs_eta_pt->Write();
}

void ActsExamples::TrackSummaryPlotTool::fill(
    TrackSummaryPlotTool::TrackSummaryPlotCache& trackSummaryPlotCache,
    const Acts::BoundTrackParameters& fittedParameters, std::size_t nStates,
    std::size_t nMeasurements, std::size_t nOutliers, std::size_t nHoles,
    std::size_t nSharedHits) const {
  using Acts::VectorHelpers::eta;
  using Acts::VectorHelpers::perp;
  const auto momentum = fittedParameters.momentum();
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  PlotHelpers::fillProf(trackSummaryPlotCache.nStates_vs_eta, fit_eta, nStates);
  PlotHelpers::fillProf(trackSummaryPlotCache.nMeasurements_vs_eta, fit_eta,
                        nMeasurements);
  PlotHelpers::fillProf(trackSummaryPlotCache.nOutliers_vs_eta, fit_eta,
                        nOutliers);
  PlotHelpers::fillProf(trackSummaryPlotCache.nHoles_vs_eta, fit_eta, nHoles);
  PlotHelpers::fillProf(trackSummaryPlotCache.nSharedHits_vs_eta, fit_eta,
                        nSharedHits);

  PlotHelpers::fillProf(trackSummaryPlotCache.nStates_vs_pt, fit_pT, nStates);
  PlotHelpers::fillProf(trackSummaryPlotCache.nMeasurements_vs_pt, fit_pT,
                        nMeasurements);
  PlotHelpers::fillProf(trackSummaryPlotCache.nOutliers_vs_pt, fit_pT,
                        nOutliers);
  PlotHelpers::fillProf(trackSummaryPlotCache.nHoles_vs_pt, fit_pT, nHoles);
  PlotHelpers::fillProf(trackSummaryPlotCache.nSharedHits_vs_pt, fit_pT,
                        nSharedHits);

  PlotHelpers::fillProf2D(trackSummaryPlotCache.nStates_vs_eta, fit_eta, fit_pT,
                          nStates);
  PlotHelpers::fillProf2D(trackSummaryPlotCache.nMeasurements_vs_eta, fit_eta,
                          fit_pT, nMeasurements);
  PlotHelpers::fillProf2D(trackSummaryPlotCache.nOutliers_vs_eta, fit_eta,
                          fit_pT, nOutliers);
  PlotHelpers::fillProf2D(trackSummaryPlotCache.nHoles_vs_eta, fit_eta, fit_pT,
                          nHoles);
  PlotHelpers::fillProf2D(trackSummaryPlotCache.nSharedHits_vs_eta, fit_eta,
                          fit_pT, nSharedHits);
}
