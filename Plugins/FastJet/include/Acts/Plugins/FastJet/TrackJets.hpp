// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Definitions/Units.hpp>

#include <optional>
#include <vector>

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

namespace Acts::FastJet {

template <typename TrackContainer>
class TrackJetSequence {
 public:
  /// Get all the track jets passing the pT & eta cuts
  ///
  /// @param ptMin the minimum jet pT in GeV
  /// @param etaMax the maximum jet absolute eta
  ///
  /// @return a vector of fastjet::PseudoJet objects
  std::vector<fastjet::PseudoJet> jets(float ptMin = 20 *
                                                     Acts::UnitConstants::GeV,
                                       float etaMax = 2.5);

  /// Get the tracks making up a track-jet
  ///
  /// @param jet the jet from which to get the constituent tracks
  /// @param coreR optional radius inside which to get the tracks
  ///
  /// @return a vector of TrackProxy
  std::vector<typename TrackContainer::TrackProxy> tracksInJet(
      const fastjet::PseudoJet& jet, std::optional<float> coreR = {});

  /// Main constructor, but using the "makeTrackJets" function is recommended
  ///
  /// @param clusterSeq the fastjet::ClusterSequence object
  /// @param inputTracks the input tracks that make up the sequence
  TrackJetSequence(fastjet::ClusterSequence clusterSeq,
                   TrackContainer& inputTracks)
      : m_clusterSeq{std::move(clusterSeq)}, m_inputTracks{inputTracks} {}

 private:
  fastjet::ClusterSequence m_clusterSeq;
  TrackContainer& m_inputTracks;
};

/// Default jet definition: Anti-kt with a radius of 0.4
const fastjet::JetDefinition DefaultJetDefinition =
    fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4);

/// Create a sequence of track jets
///
/// @param tracks the input tracks
/// @jetDef the jet definition to use, defaults to "DefaultJetDefinition"
template <typename TrackContainer>
TrackJetSequence<TrackContainer> makeTrackJets(
    TrackContainer& tracks,
    fastjet::JetDefinition jetDef = DefaultJetDefinition);

}  // namespace Acts::FastJet

#include "TrackJets.ipp"
