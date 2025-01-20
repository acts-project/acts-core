// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/FakeRatePlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <TEfficiency.h>
#include <TH2.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

ActsExamples::FakeRatePlotTool::FakeRatePlotTool(
    const ActsExamples::FakeRatePlotTool::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("FakeRatePlotTool", lvl)) {}

void ActsExamples::FakeRatePlotTool::book(
    FakeRatePlotTool::FakeRatePlotCache& fakeRatePlotCache) const {
  PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bPhi = m_cfg.varBinning.at("Phi");
  PlotHelpers::Binning bNum = m_cfg.varBinning.at("Num");
  ACTS_DEBUG("Initialize the histograms for fake rate plots");

  // number of reco tracks vs pT scatter plots
  fakeRatePlotCache.nReco_vs_pT = PlotHelpers::bookHisto(
      "nRecoTracks_vs_pT", "Number of reconstructed track candidates", bPt,
      bNum);
  // number of truth-matched tracks vs pT scatter plots
  fakeRatePlotCache.nTruthMatched_vs_pT = PlotHelpers::bookHisto(
      "nTruthMatchedTracks_vs_pT", "Number of truth-matched track candidates",
      bPt, bNum);
  // number of fake tracks vs pT scatter plots
  fakeRatePlotCache.nFake_vs_pT = PlotHelpers::bookHisto(
      "nFakeTracks_vs_pT", "Number of fake track candidates", bPt, bNum);

  // number of reco tracks vs eta scatter plots
  fakeRatePlotCache.nReco_vs_eta = PlotHelpers::bookHisto(
      "nRecoTracks_vs_eta", "Number of reconstructed track candidates", bEta,
      bNum);
  // number of truth-matched tracks vs eta scatter plots
  fakeRatePlotCache.nTruthMatched_vs_eta = PlotHelpers::bookHisto(
      "nTruthMatchedTracks_vs_eta", "Number of truth-matched track candidates",
      bEta, bNum);
  // number of fake tracks vs eta scatter plots
  fakeRatePlotCache.nFake_vs_eta = PlotHelpers::bookHisto(
      "nFakeTracks_vs_eta", "Number of fake track candidates", bEta, bNum);

  // fake rate vs pT
  fakeRatePlotCache.fakeRate_vs_pT = PlotHelpers::bookEff(
      "fakerate_vs_pT", "Tracking fake rate;p_{T} [GeV/c];Fake rate", bPt);
  // fake rate vs eta
  fakeRatePlotCache.fakeRate_vs_eta = PlotHelpers::bookEff(
      "fakerate_vs_eta", "Tracking fake rate;#eta;Fake rate", bEta);
  // fake rate vs phi
  fakeRatePlotCache.fakeRate_vs_phi = PlotHelpers::bookEff(
      "fakerate_vs_phi", "Tracking fake rate;#phi;Fake rate", bPhi);
  // fake rate vs pT and eta
  fakeRatePlotCache.fakeRate_vs_pT_eta = PlotHelpers::bookEff(
      "fakerate_vs_pT_eta", "Tracking fake rate;p_{T] [GeV/c];#eta;Fake rate",
      bPt);
}

void ActsExamples::FakeRatePlotTool::clear(
    FakeRatePlotCache& fakeRatePlotCache) const {
  delete fakeRatePlotCache.nReco_vs_pT;
  delete fakeRatePlotCache.nTruthMatched_vs_pT;
  delete fakeRatePlotCache.nFake_vs_pT;
  delete fakeRatePlotCache.nReco_vs_eta;
  delete fakeRatePlotCache.nTruthMatched_vs_eta;
  delete fakeRatePlotCache.nFake_vs_eta;
  delete fakeRatePlotCache.fakeRate_vs_pT;
  delete fakeRatePlotCache.fakeRate_vs_eta;
  delete fakeRatePlotCache.fakeRate_vs_phi;
  delete fakeRatePlotCache.fakeRate_vs_pT_eta;
}

void ActsExamples::FakeRatePlotTool::write(
    const FakeRatePlotTool::FakeRatePlotCache& fakeRatePlotCache) const {
  ACTS_DEBUG("Write the plots to output file.");
  fakeRatePlotCache.nReco_vs_pT->Write();
  fakeRatePlotCache.nTruthMatched_vs_pT->Write();
  fakeRatePlotCache.nFake_vs_pT->Write();
  fakeRatePlotCache.nReco_vs_eta->Write();
  fakeRatePlotCache.nTruthMatched_vs_eta->Write();
  fakeRatePlotCache.nFake_vs_eta->Write();
  fakeRatePlotCache.fakeRate_vs_pT->Write();
  fakeRatePlotCache.fakeRate_vs_eta->Write();
  fakeRatePlotCache.fakeRate_vs_phi->Write();
  fakeRatePlotCache.fakeRate_vs_pT_eta->Write();
}

void ActsExamples::FakeRatePlotTool::fill(
    FakeRatePlotTool::FakeRatePlotCache& fakeRatePlotCache,
    const Acts::BoundTrackParameters& fittedParameters, bool status) const {
  const auto momentum = fittedParameters.momentum();
  const double fit_phi = phi(momentum);
  const double fit_eta = eta(momentum);
  const double fit_pT = perp(momentum);

  PlotHelpers::fillEff(fakeRatePlotCache.fakeRate_vs_pT, fit_pT, status);
  PlotHelpers::fillEff(fakeRatePlotCache.fakeRate_vs_eta, fit_eta, status);
  PlotHelpers::fillEff(fakeRatePlotCache.fakeRate_vs_phi, fit_phi, status);
  PlotHelpers::fillEff(fakeRatePlotCache.fakeRate_vs_pT_eta, fit_pT, fit_eta,
                       status);
}

void ActsExamples::FakeRatePlotTool::fill(
    FakeRatePlotTool::FakeRatePlotCache& fakeRatePlotCache,
    const SimParticleState& truthParticle, std::size_t nTruthMatchedTracks,
    std::size_t nFakeTracks) const {
  const auto t_eta = eta(truthParticle.direction());
  const auto t_pT = truthParticle.transverseMomentum();

  PlotHelpers::fillHisto(fakeRatePlotCache.nReco_vs_pT, t_pT,
                         nTruthMatchedTracks + nFakeTracks);
  PlotHelpers::fillHisto(fakeRatePlotCache.nTruthMatched_vs_pT, t_pT,
                         nTruthMatchedTracks);
  PlotHelpers::fillHisto(fakeRatePlotCache.nFake_vs_pT, t_pT, nFakeTracks);

  PlotHelpers::fillHisto(fakeRatePlotCache.nReco_vs_eta, t_eta,
                         nTruthMatchedTracks + nFakeTracks);
  PlotHelpers::fillHisto(fakeRatePlotCache.nTruthMatched_vs_eta, t_eta,
                         nTruthMatchedTracks);
  PlotHelpers::fillHisto(fakeRatePlotCache.nFake_vs_eta, t_eta, nFakeTracks);
}
