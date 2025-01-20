// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <Acts/Plugins/FastJet/TrackJets.hpp>

class ParticleHypothesis {
 public:
  float mass() { return 139.57061 * Acts::UnitConstants::MeV; }
};

class Track {
 public:
  Track(float pt, float eta, float phi) {
    m_momentum[0] = pt * std::cos(phi);
    m_momentum[1] = pt * std::sin(phi);
    m_momentum[1] = pt * std::sinh(eta);
  }
  Acts::Vector3 momentum() const { return m_momentum; }
  ParticleHypothesis particleHypothesis() const { return ParticleHypothesis(); }

 private:
  Acts::Vector3 m_momentum;
};

bool operator==(Track const& lhs, Track const& rhs) {
  return lhs.momentum() == rhs.momentum();
}

class TrackContainer {
 public:
  using TrackProxy = Track;

  TrackContainer() {}
  void insert(Track track) { m_vec.push_back(std::move(track)); }
  std::size_t size() { return m_vec.size(); }

  using ConstTrackProxy = const Track&;
  ConstTrackProxy getTrack(std::size_t i) {
    if (i < size()) {
      return m_vec[i];
    }
    throw std::runtime_error("Too few tracks");
  }

 private:
  std::vector<Track> m_vec;
};

BOOST_AUTO_TEST_CASE(SingleTrack) {
  TrackContainer tracks;
  tracks.insert(Track(100, 0, 0));

  Acts::FastJet::TrackJetSequence jetSeq = Acts::FastJet::makeTrackJets(tracks);
  std::vector<fastjet::PseudoJet> jets = jetSeq.jets();

  BOOST_CHECK_EQUAL(jets.size(), 1);
  BOOST_CHECK_EQUAL(jets[0].constituents().size(), 1);
  BOOST_CHECK_EQUAL(jets[0].constituents()[0].user_index(), 0);
  BOOST_CHECK_CLOSE(jets[0].pt(), 100, 1e-3);
  BOOST_CHECK_CLOSE(jets[0].eta(), 0, 1e-3);
  BOOST_CHECK_CLOSE(jets[0].phi(), 0, 1e-3);
  BOOST_CHECK_CLOSE(jets[0].m(), ParticleHypothesis().mass(), 1);
}

BOOST_AUTO_TEST_CASE(TwoTracksTwoJets) {
  TrackContainer tracks;
  tracks.insert(Track(100, 0, 0.0));
  tracks.insert(Track(100, 0, M_PI));

  Acts::FastJet::TrackJetSequence jetSeq = Acts::FastJet::makeTrackJets(tracks);
  std::vector<fastjet::PseudoJet> jets = jetSeq.jets();

  BOOST_CHECK_EQUAL(jets.size(), 2);

  std::vector<Track> trks_0 = jetSeq.tracksInJet(jets[0]);
  BOOST_CHECK_EQUAL(trks_0.size(), 1);
  BOOST_CHECK(trks_0[0] == tracks.getTrack(0) ||
              trks_0[0] == tracks.getTrack(1));

  std::vector<Track> trks_1 = jetSeq.tracksInJet(jets[1]);
  BOOST_CHECK_EQUAL(trks_1.size(), 1);
  BOOST_CHECK(trks_1[0] == tracks.getTrack(0) ||
              trks_1[0] == tracks.getTrack(1));
  BOOST_CHECK(trks_0[0] != trks_1[0]);
}

BOOST_AUTO_TEST_CASE(TwoTracksOneJet) {
  TrackContainer tracks;
  tracks.insert(Track(100, 0, 0.0));
  tracks.insert(Track(100, 0, 0.2));

  Acts::FastJet::TrackJetSequence jetSeq = Acts::FastJet::makeTrackJets(tracks);
  std::vector<fastjet::PseudoJet> jets = jetSeq.jets();

  BOOST_CHECK_EQUAL(jets.size(), 1);

  std::vector<Track> trks_0 = jetSeq.tracksInJet(jets[0]);
  BOOST_CHECK_EQUAL(trks_0.size(), 2);
  BOOST_CHECK(trks_0[0] == tracks.getTrack(0) ||
              trks_0[0] == tracks.getTrack(1));
  BOOST_CHECK(trks_0[1] == tracks.getTrack(0) ||
              trks_0[1] == tracks.getTrack(1));
  BOOST_CHECK(trks_0[0] != trks_0[1]);
}

BOOST_AUTO_TEST_CASE(TracksInJetCore) {
  TrackContainer tracks;
  tracks.insert(Track(100, 0, 0));
  tracks.insert(Track(10, 0.05, 0));
  tracks.insert(Track(10, -0.05, 0));
  tracks.insert(Track(10, 0.2, 0));
  tracks.insert(Track(10, -0.2, 0));

  Acts::FastJet::TrackJetSequence jetSeq = Acts::FastJet::makeTrackJets(tracks);
  std::vector<fastjet::PseudoJet> jets = jetSeq.jets();

  BOOST_REQUIRE_EQUAL(jets.size(), 1);

  std::vector<Track> trks = jetSeq.tracksInJet(jets[0], 0.1);
  BOOST_CHECK_EQUAL(trks.size(), 3);

  BOOST_CHECK(std::find(trks.begin(), trks.end(), tracks.getTrack(0)) !=
              trks.end());
  BOOST_CHECK(std::find(trks.begin(), trks.end(), tracks.getTrack(1)) !=
              trks.end());
  BOOST_CHECK(std::find(trks.begin(), trks.end(), tracks.getTrack(2)) !=
              trks.end());
  BOOST_CHECK(std::find(trks.begin(), trks.end(), tracks.getTrack(3)) ==
              trks.end());
  BOOST_CHECK(std::find(trks.begin(), trks.end(), tracks.getTrack(4)) ==
              trks.end());
}
