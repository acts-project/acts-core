// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Common.hpp>

template <typename TrackContainer>
Acts::FastJet::TrackJetSequence<TrackContainer> Acts::FastJet::makeTrackJets(
    TrackContainer& tracks, fastjet::JetDefinition jetDef) {
  std::vector<fastjet::PseudoJet> inputs;

  for (std::size_t i = 0; i < tracks.size(); i++) {
    typename TrackContainer::ConstTrackProxy track = tracks.getTrack(i);
    Acts::Vector3 p = track.momentum();
    float m = track.particleHypothesis().mass();

    float px = p[Acts::eMom0];
    float py = p[Acts::eMom1];
    float pz = p[Acts::eMom2];
    float e = std::sqrt(m * m + px * px + py * py + pz * pz);

    inputs.emplace_back(px, py, pz, e);
    inputs.back().set_user_index(i);
  }

  fastjet::ClusterSequence cs(inputs, jetDef);

  return TrackJetSequence(std::move(cs), tracks);
}

template <typename TrackContainer>
std::vector<fastjet::PseudoJet>
Acts::FastJet::TrackJetSequence<TrackContainer>::jets(float ptMin,
                                                      float etaMax) {
  fastjet::Selector sel_eta = fastjet::SelectorAbsEtaMax(etaMax);
  return sel_eta(m_clusterSeq.inclusive_jets(ptMin));
}

template <typename TrackContainer>
std::vector<typename TrackContainer::TrackProxy>
Acts::FastJet::TrackJetSequence<TrackContainer>::tracksInJet(
    const fastjet::PseudoJet& jet, std::optional<float> coreR) {
  fastjet::Selector sel = fastjet::SelectorIdentity();
  if (coreR.has_value()) {
    sel = fastjet::SelectorCircle(coreR.value());
    sel.set_reference(jet);
  }

  std::vector<typename TrackContainer::TrackProxy> tracks;
  for (fastjet::PseudoJet& cst : sel(jet.constituents())) {
    tracks.push_back(m_inputTracks.getTrack(cst.user_index()));
  }

  return tracks;
}
