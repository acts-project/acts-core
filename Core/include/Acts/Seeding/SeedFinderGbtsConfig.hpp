// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// TODO: update to C++17 style

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/GbtsBase.hpp"  //definition of Trigsispacepoint base and trigtriplets
#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"

#include <memory>

// core algorithm so in acts namespace
namespace Acts {

template <typename SpacePoint>
struct SeedFinderGbtsConfig {
  // // how many sigmas of scattering angle should be considered?
  float sigmaScattering = 5;

  // Seed cut
  float minPt = 400. * Acts::UnitConstants::MeV;

  //   //detector ROI
  //   // derived values, set on SeedFinder construction
  float highland = 0;
  float maxScatteringAngle2 = 0;
  //   bool isInInternalUnits = false;
  /// for load space points

  // Parameter which can loosen the tolerance of the track seed to form a
  // helix. This is useful for e.g. misaligned seeding.
  float helixCutTolerance = 1.;

  float m_phiSliceWidth{};    // initialised in loadSpacePoints function
  float m_nMaxPhiSlice = 53;  // used to calculate phi slices
  bool m_useClusterWidth =
      false;  // bool for use of cluster width in loadSpacePoints function
  std::string ConnectorInputFile;  // Path to the connector configuration file
                                   // that defines the layer connections
  std::vector<TrigInDetSiLayer> m_layerGeometry;

  // for runGbts_TrackFinder
  bool m_LRTmode = true;  // eventually want to set from full chain
  bool m_useEtaBinning =
      true;  // bool to use eta binning from geometry structure
  bool m_doubletFilterRZ = true;  // bool applies new Z cuts on doublets
  float m_minDeltaRadius = 2.0;   // min dr for doublet
  float m_tripletD0Max = 4.0;     // D0 cut for triplets
  unsigned int m_maxTripletBufferLength =
      3;                        // maximum number of space points per triplet
  int MaxEdges = 2000000;       // max number of Gbts edges/doublets
  float cut_dphi_max = 0.012;   // phi cut for triplets
  float cut_dcurv_max = 0.001;  // curv cut for triplets
  float cut_tau_ratio_max = 0.007;  // tau cut for doublets and triplets
  float maxOuterRadius = 550.0;     // used to calculate Z cut on doublets
  float m_PtMin = 1000.0;
  float m_tripletPtMinFrac = 0.3;
  float m_tripletPtMin = m_PtMin * m_tripletPtMinFrac;  // Limit on triplet pt
  double ptCoeff =
      0.29997 * 1.9972 / 2.0;  // ~0.3*B/2 - assumes nominal field of 2*T

  // ROI:
  bool containsPhi() {
    return false;
    // need to implement this function
  }

  ////
  // 2 member functions
  SeedFinderGbtsConfig calculateDerivedQuantities() const {
    // thorw statement if the isInternalUnits member is false, ie if dont call
    // this function
    SeedFinderGbtsConfig config = *this;
    // use a formula to calculate scattering

    return config;
  }

  SeedFinderGbtsConfig toInternalUnits() const {
    // throw statement if the isInternalUnits member is false, ie if dont call
    // this function
    SeedFinderGbtsConfig config = *this;
    // divides inputs by 1mm, all ones input
    // changes member inInInternalUnits to true
    return config;
  }

};  // end of config struct

}  // namespace Acts
