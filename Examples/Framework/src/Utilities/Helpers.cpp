// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/Helpers.hpp"

#include <cassert>

#include <TAxis.h>
#include <TEfficiency.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

namespace ActsExamples::PlotHelpers {
TH1F* bookHisto(const char* histName, const char* histTitle,
                const Binning& varBinning) {
  TH1F* hist =
      new TH1F(histName, histTitle, varBinning.nBins(), varBinning.data());
  hist->GetXaxis()->SetTitle(varBinning.title().c_str());
  hist->GetYaxis()->SetTitle("Entries");
  hist->Sumw2();
  return hist;
}

TH2F* bookHisto(const char* histName, const char* histTitle,
                const Binning& varXBinning, const Binning& varYBinning) {
  TH2F* hist =
      new TH2F(histName, histTitle, varXBinning.nBins(), varXBinning.data(),
               varYBinning.nBins(), varYBinning.data());
  hist->GetXaxis()->SetTitle(varXBinning.title().c_str());
  hist->GetYaxis()->SetTitle(varYBinning.title().c_str());
  hist->Sumw2();
  return hist;
}

void fillHisto(TH1F* hist, float value, float weight) {
  assert(hist != nullptr);
  hist->Fill(value, weight);
}

void fillHisto(TH2F* hist, float xValue, float yValue, float weight) {
  assert(hist != nullptr);
  hist->Fill(xValue, yValue, weight);
}

void anaHisto(TH1D* inputHist, int j, TH1F* meanHist, TH1F* widthHist) {
  // evaluate mean and width via the Gauss fit
  assert(inputHist != nullptr);
  if (inputHist->GetEntries() > 0) {
    TFitResultPtr r = inputHist->Fit("gaus", "QS0");
    if ((r.Get() != nullptr) && ((r->Status() % 1000) == 0)) {
      // fill the mean and width into 'j'th bin of the meanHist and widthHist,
      // respectively
      meanHist->SetBinContent(j, r->Parameter(1));
      meanHist->SetBinError(j, r->ParError(1));
      widthHist->SetBinContent(j, r->Parameter(2));
      widthHist->SetBinError(j, r->ParError(2));
    }
  }
}

TEfficiency* bookEff(const char* effName, const char* effTitle,
                     const Binning& varBinning) {
  TEfficiency* efficiency =
      new TEfficiency(effName, effTitle, varBinning.nBins(), varBinning.data());
  return efficiency;
}

TEfficiency* bookEff(const char* effName, const char* effTitle,
                     const Binning& varXBinning, const Binning& varYBinning) {
  TEfficiency* efficiency = new TEfficiency(
      effName, effTitle, varXBinning.nBins(), varXBinning.data(),
      varYBinning.nBins(), varYBinning.data());
  return efficiency;
}

void fillEff(TEfficiency* efficiency, float value, bool status) {
  assert(efficiency != nullptr);
  efficiency->Fill(status, value);
}

void fillEff(TEfficiency* efficiency, float xValue, float yValue, bool status) {
  assert(efficiency != nullptr);
  efficiency->Fill(status, xValue, yValue);
}

TProfile* bookProf(const char* profName, const char* profTitle,
                   const Binning& varXBinning, const Binning& varYBinning) {
  TProfile* prof =
      new TProfile(profName, profTitle, varXBinning.nBins(), varXBinning.data(),
                   varYBinning.low(), varYBinning.high());
  prof->GetXaxis()->SetTitle(varXBinning.title().c_str());
  prof->GetYaxis()->SetTitle(varYBinning.title().c_str());
  return prof;
}

TProfile2D* bookProf2D(const char* profName, const char* profTitle,
                       const Binning& varXBinning, const Binning& varYBinning,
                       const Binning& varZBinning) {
  TProfile2D* prof =
      new TProfile2D(profName, profTitle, varXBinning.nBins, varXBinning.min,
                     varXBinning.max, varYBinning.nBins, varYBinning.min,
                     varYBinning.max, varZBinning.min, varZBinning.max);
  prof->GetXaxis()->SetTitle(varXBinning.title.c_str());
  prof->GetYaxis()->SetTitle(varYBinning.title.c_str());
  prof->GetZaxis()->SetTitle(varZBinning.title.c_str());
  return prof;
}

void fillProf(TProfile* profile, float xValue, float yValue, float weight) {
  assert(profile != nullptr);
  profile->Fill(xValue, yValue, weight);
}

void fillProf2D(TProfile2D* profile, float xValue, float yValue, float zValue,
                float weight) {
  assert(profile != nullptr);
  profile->Fill(xValue, yValue, zValue, weight);
}

}  // namespace ActsExamples::PlotHelpers
