// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
#include <limits>

#include <TColor.h>
#include <TDirectory.h>
#include <TH1F.h>
#include <TString.h>

// Helper function:
// function to set up the histogram style
template <typename hist_t>
void setHistStyle(hist_t* hist, short color = 1) {
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.);
  hist->GetYaxis()->SetTitleOffset(1.8);
  hist->GetXaxis()->SetNdivisions(505);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.8);
  hist->SetLineWidth(2);
  hist->SetTitle("");
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
}

// Helper function:
// function to set up the efficiency histogram style
template <typename eff_t>
void setEffStyle(eff_t* eff, short color = 1) {
  eff->SetMarkerStyle(20);
  eff->SetMarkerSize(0.8);
  eff->SetLineWidth(2);
  eff->SetLineColor(color);
  eff->SetMarkerColor(color);
}

// Helper function:
// set color pallette
template <typename hist_t>
void adaptColorPalette(hist_t* h, float rmin, float rmax, float rgood,
                       float rwindow, int n) {
  // min - max is the range of the axis
  float rel_good = (rgood - rmin) / (rmax - rmin);
  float rel_window = rwindow / (rmax - rmin);

  // Stops are
  const int number = 5;
  double red[number] = {0., 0., 0., 1., 1.};
  double green[number] = {0., 1., 1., 1., 0.};
  double blue[number] = {1., 1., 0., 0., 0.};
  double stops[number] = {0., rel_good - rel_window, rel_good,
                          rel_good + rel_window, 1.};
  h->SetContour(n);

  TColor::CreateGradientColorTable(number, stops, red, green, blue, n);
}

// Helper function:
// increase eff range by a scale factor. Note that it assumes the eff has
// already been drawn
template <typename eff_t>
void adaptEffRange(eff_t* eff, float minScale = 1, float maxScale = 1.1) {
  gPad->Update();
  auto ymin = gPad->GetUymin();
  auto ymax = gPad->GetUymax();
  auto graph = eff->GetPaintedGraph();
  graph->SetMinimum(ymin * minScale);
  graph->SetMaximum(ymax * maxScale);
  gPad->Modified();
  gPad->Update();
}

struct ParameterHandle {
  /// A tag name
  std::string tag = "";

  /// Title and names: residual
  std::string residualStr = "";
  std::string residualUnit = "";

  /// Title and names: error
  std::string errorStr = "";

  /// The rangeDrawStr draw string
  std::string rangeDrawStr = "";
  std::string rangeMaxStr = "";
  std::string rangeCutStr = "";

  /// The range array
  std::array<float, 2> range = {0., 0.};

  /// Value function that allows to create
  /// combined parameters
  std::function<float(ULong64_t)> value;

  /// The associated error accessor
  std::function<float(ULong64_t)> error;

  /// The acceptance
  std::function<bool(ULong64_t)> accept;

  TH1F* rangeHist = nullptr;

  TH1F* residualHist = nullptr;

  TH1F* pullHist = nullptr;

  ULong64_t accepted = 0;

  // Fill the entry
  void fill(unsigned int entry) {
    if (accept(entry)) {
      // Access the value, error
      float v = value(entry);
      residualHist->Fill(v);
      pullHist->Fill(v / error(entry));
      // Count the accessor
      ++accepted;
    }
  };
};

// This is a combined acceptor
struct AcceptCombination {
  std::function<bool(ULong64_t)> one;

  std::function<bool(ULong64_t)> two;

  /// returns true if value is within range
  /// @param entry the entry in the tree
  bool operator()(ULong64_t entry) { return (one(entry) and two(entry)); }
};

// This Struct is to accept all values
struct AcceptAll {
  // Call operator always returns true
  bool operator()(ULong64_t /*event*/) { return true; }
};

// This Struct is to accept a certain range
struct AcceptRange {
  std::vector<float>* value = nullptr;

  std::array<float, 2> range = {0., 0.};

  /// returns true if value is within range
  /// @param entry the entry in the tree
  bool operator()(ULong64_t entry) {
    if (value) {
      float v = value->at(entry);
      return (range[0] <= v and range[1] > v);
    }
    return false;
  }
};

// This is a direct type accessor
struct DirectAccessor {
  std::vector<float>* value = nullptr;

  /// returns the calculated Residual
  /// @param entry the entry in the tree
  float operator()(ULong64_t entry) {
    if (value) {
      float v = value->at(entry);
      return v;
    }
    return std::numeric_limits<float>::infinity();
  }
};

// This is a residual type accessor
struct ResidualAccessor {
  std::vector<float>* value = nullptr;

  std::vector<float>* reference = nullptr;

  /// returns the calculated Residual
  /// @param entry the entry in the tree
  float operator()(ULong64_t entry) {
    if (value and reference) {
      float v = value->at(entry);
      float r = reference->at(entry);
      return (v - r);
    }
    return std::numeric_limits<float>::infinity();
  }
};

// This is a  qop residual accessor
struct QopResidualAccessor {
  std::vector<float>* qop_value = nullptr;

  std::vector<int>* reference_charge = nullptr;

  std::vector<float>* reference_p = nullptr;

  /// returns the calculated Residual
  /// @param entry the entry in the tree
  float operator()(ULong64_t entry) {
    if (qop_value and reference_charge and reference_p) {
      float v = qop_value->at(entry);
      float q_true = reference_charge->at(entry);
      float p_true = reference_p->at(entry);
      return (v - q_true / p_true);
    }
    return std::numeric_limits<float>::infinity();
  }
};

/// This the pT residual accessor
struct PtResidualAccessor {
  std::vector<float>* qop_value = nullptr;

  std::vector<float>* theta_value = nullptr;

  std::vector<float>* reference_pt = nullptr;

  /// returns the calculated Residual
  /// @param entry the entry in the tree
  float operator()(ULong64_t entry) {
    if (qop_value and theta_value and reference_pt) {
      float p = 1. / std::abs(qop_value->at(entry));
      float theta = theta_value->at(entry);
      float pt_true = reference_pt->at(entry);
      return (p * sin(theta) - pt_true);
    }
    return std::numeric_limits<float>::infinity();
  }
};

// This is a  qop residual accessor
struct PtErrorAccessor {
  std::vector<float>* qop_value = nullptr;
  std::vector<float>* qop_error = nullptr;

  std::vector<float>* theta_value = nullptr;
  std::vector<float>* theta_error = nullptr;

  /// returns the calculated Residual
  /// @param entry the entry in the tree
  float operator()(ULong64_t entry) {
    if (qop_value and qop_error and theta_value and theta_error) {
      float qop_v = qop_value->at(entry);
      float qop_e = qop_error->at(entry);
      float theta_v = theta_value->at(entry);
      float theta_e = theta_error->at(entry);
      return std::cos(theta_v) / qop_v * theta_e -
             std::sin(theta_v) / (qop_v * qop_v) * qop_e;
    }
    return std::numeric_limits<float>::infinity();
  }
};

/// Range estimation
template <typename dir_t, typename tree_t>
void estimateRange(ParameterHandle& handle, dir_t& directory, tree_t& tree,
                   unsigned long peakEntries, unsigned int hBarcode) {
  // Change into the Directory
  directory.cd();
  TString rangeHist = handle.rangeDrawStr;
  rangeHist += ">>";
  // Hist name snipped
  TString rangeHN = "hrg_";
  rangeHN += hBarcode;
  // Full histogram
  rangeHist += rangeHN;
  rangeHist += handle.rangeMaxStr;

  // Do the drawing
  tree.Draw(rangeHist.Data(), handle.rangeCutStr.c_str(), "", peakEntries);
  handle.rangeHist = dynamic_cast<TH1F*>(gDirectory->Get(rangeHN.Data()));
  if (handle.rangeHist != nullptr) {
    float rms = handle.rangeHist->GetRMS();
    handle.range = {-rms, rms};
  }
}

void bookHistograms(ParameterHandle& handle, float pullRange,
                    unsigned int hBins, unsigned int hBarcode) {
  // Residual histogram
  TString rName = std::string("res_") + handle.tag;
  rName += hBarcode;
  handle.residualHist =
      new TH1F(rName.Data(), handle.tag.c_str(), hBins,
               pullRange * handle.range[0], pullRange * handle.range[1]);
  std::string xAxisTitle =
      handle.residualStr + std::string(" ") + handle.residualUnit;
  handle.residualHist->GetXaxis()->SetTitle(xAxisTitle.c_str());
  handle.residualHist->GetYaxis()->SetTitle("Entries");

  // Pull histogram
  TString pName = std::string("pull_") + handle.tag;
  pName += hBarcode;
  handle.pullHist =
      new TH1F(pName.Data(), (std::string("pull ") + handle.tag).c_str(), hBins,
               -pullRange, pullRange);
  xAxisTitle = std::string("(") + handle.residualStr + std::string(")/") +
               handle.errorStr;
  handle.pullHist->GetXaxis()->SetTitle(xAxisTitle.c_str());
  handle.pullHist->GetYaxis()->SetTitle("Entries");
}

/// Helper method to get and opentially overwrite the entries to be processed
///
/// @param tree is the TTree/TChain in question
/// @param configuredEntries is a configuraiton parameter
///
/// @return the number of entries
template <typename tree_t>
unsigned long estimateEntries(const tree_t& tree,
                              unsigned long configuredEntries) {
  unsigned long entries = static_cast<unsigned long>(tree.GetEntries());
  if (configuredEntries > 0 and configuredEntries < entries) {
    entries = configuredEntries;
  }
  return entries;
}
