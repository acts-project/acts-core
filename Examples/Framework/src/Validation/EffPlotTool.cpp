// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/EffPlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <TEfficiency.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

ActsExamples::EffPlotTool::EffPlotTool(
    const ActsExamples::EffPlotTool::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("EffPlotTool", lvl)) {}

void ActsExamples::EffPlotTool::book(
    EffPlotTool::EffPlotCache& effPlotCache) const {
  PlotHelpers::Binning bPhi = m_cfg.varBinning.at("Phi");
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
  PlotHelpers::Binning bDeltaR = m_cfg.varBinning.at("DeltaR");
  PlotHelpers::Binning bZ0 = m_cfg.varBinning.at("Z0");
  ACTS_DEBUG("Initialize the histograms for efficiency plots");
  // efficiency vs pT
  effPlotCache.trackEff_vs_pT = PlotHelpers::bookEff(
      "trackeff_vs_pT", "Tracking efficiency;Truth p_{T} [GeV/c];Efficiency",
      bPt);
  // efficiency vs eta
  effPlotCache.trackEff_vs_eta = PlotHelpers::bookEff(
      "trackeff_vs_eta", "Tracking efficiency;Truth #eta;Efficiency", bEta);
  // efficiency vs phi
  effPlotCache.trackEff_vs_phi = PlotHelpers::bookEff(
      "trackeff_vs_phi", "Tracking efficiency;Truth #phi;Efficiency", bPhi);
  // efficiency vs z0
  effPlotCache.trackEff_vs_z0 = PlotHelpers::bookEff(
      "trackeff_vs_z0", "Tracking efficiency;Truth z_{0} [mm];Efficiency", bZ0);
  // efficiancy vs distance to the closest truth particle
  effPlotCache.trackEff_vs_DeltaR = PlotHelpers::bookEff(
      "trackeff_vs_DeltaR",
      "Tracking efficiency;Closest track #Delta R;Efficiency", bDeltaR);
  // efficiency vs pT and eta
  effPlotCache.track2dEff_vs_pT_eta = PlotHelpers::bookEff(
      "trackeff_vs_pT_eta",
      "Tracking efficiency;Truth p_{T} [GeV/c];Truth #eta;Efficiency", bPt,
      bEta);
}

void ActsExamples::EffPlotTool::clear(EffPlotCache& effPlotCache) const {
  delete effPlotCache.trackEff_vs_pT;
  delete effPlotCache.trackEff_vs_eta;
  delete effPlotCache.trackEff_vs_phi;
  delete effPlotCache.trackEff_vs_z0;
  delete effPlotCache.trackEff_vs_DeltaR;
  delete effPlotCache.track2dEff_vs_pT_eta;
}

void ActsExamples::EffPlotTool::write(
    const EffPlotTool::EffPlotCache& effPlotCache) const {
  ACTS_DEBUG("Write the plots to output file.");
  effPlotCache.trackEff_vs_pT->Write();
  effPlotCache.trackEff_vs_eta->Write();
  effPlotCache.trackEff_vs_phi->Write();
  effPlotCache.trackEff_vs_z0->Write();
  effPlotCache.trackEff_vs_DeltaR->Write();
  effPlotCache.track2dEff_vs_pT_eta->Write();
}

void ActsExamples::EffPlotTool::fill(EffPlotTool::EffPlotCache& effPlotCache,
                                     const SimParticleState& truthParticle,
                                     double deltaR, bool status) const {
  const auto t_phi = phi(truthParticle.direction());
  const auto t_eta = eta(truthParticle.direction());
  const auto t_pT = truthParticle.transverseMomentum();
  const auto t_z0 = truthParticle.position().z();
  const auto t_deltaR = deltaR;

  PlotHelpers::fillEff(effPlotCache.trackEff_vs_pT, t_pT, status);
  PlotHelpers::fillEff(effPlotCache.trackEff_vs_eta, t_eta, status);
  PlotHelpers::fillEff(effPlotCache.trackEff_vs_phi, t_phi, status);
  PlotHelpers::fillEff(effPlotCache.trackEff_vs_z0, t_z0, status);
  PlotHelpers::fillEff(effPlotCache.trackEff_vs_DeltaR, t_deltaR, status);
  PlotHelpers::fillEff(effPlotCache.track2dEff_vs_pT_eta, t_pT, t_eta, status);
}
