// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "SeedingPerformanceWriter.hpp"

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <stdexcept>
#include <unordered_map>
#include <iostream>
#include <fstream>

#include <TFile.h>

namespace {
using SimParticleContainer = ActsExamples::SimParticleContainer;
using HitParticlesMap = ActsExamples::IndexMultimap<ActsFatras::Barcode>;
using ProtoTrackContainer = ActsExamples::ProtoTrackContainer;
}  // namespace

ActsExamples::SeedingPerformanceWriter::SeedingPerformanceWriter(
    ActsExamples::SeedingPerformanceWriter::Config config,
    Acts::Logging::Level level)
    : WriterT(config.inputProtoTracks, "SeedingPerformanceWriter", level),
      m_cfg(std::move(config)),
      m_effPlotTool(m_cfg.effPlotToolConfig, level),
      m_duplicationPlotTool(m_cfg.duplicationPlotToolConfig, level) {
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing hit-particles map input collection");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles collection");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  // the output file can not be given externally since TFile accesses to the
  // same file from multiple threads are unsafe.
  // must always be opened internally
  auto path = m_cfg.filePath;
  m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
  if (not m_outputFile) {
    throw std::invalid_argument("Could not open '" + path + "'");
  }
  // initialize the plot tools
  m_effPlotTool.book(m_effPlotCache);
  m_duplicationPlotTool.book(m_duplicationPlotCache);

  std::ofstream myfile;
  myfile.open("eff.csv");
  myfile << "event,phi,eta,pT,matched\n";
  myfile.close();
}

ActsExamples::SeedingPerformanceWriter::~SeedingPerformanceWriter() {
  m_effPlotTool.clear(m_effPlotCache);
  m_duplicationPlotTool.clear(m_duplicationPlotCache);
  if (m_outputFile) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::SeedingPerformanceWriter::endRun() {
  float eff = float(m_nTotalMatchedParticles) / m_nTotalParticles;
  float fakeRate = float(m_nTotalSeeds - m_nTotalMatchedSeeds) / m_nTotalSeeds;
  float duplicationRate =
      float(m_nTotalDuplicatedParticles) / m_nTotalMatchedParticles;
  float aveNDuplicatedSeeds =
      float(m_nTotalMatchedSeeds - m_nTotalMatchedParticles) /
      m_nTotalMatchedParticles;
  ACTS_DEBUG("nTotalSeeds               = " << m_nTotalSeeds);
  ACTS_DEBUG("nTotalMatchedSeeds        = " << m_nTotalMatchedSeeds);
  ACTS_DEBUG("nTotalParticles           = " << m_nTotalParticles);
  ACTS_DEBUG("nTotalMatchedParticles    = " << m_nTotalMatchedParticles);
  ACTS_DEBUG("nTotalDuplicatedParticles = " << m_nTotalDuplicatedParticles);

  ACTS_INFO("Efficiency (nMatchedParticles / nAllParticles) = " << eff << " (" << m_nTotalMatchedParticles << " / " << m_nTotalParticles << ")");
  ACTS_INFO("Fake rate (nUnMatchedSeeds / nAllSeeds) = " << fakeRate << " (" << (m_nTotalSeeds - m_nTotalMatchedSeeds) << " / " << m_nTotalSeeds << ")");
  ACTS_INFO(
      "Duplication rate (nDuplicatedMatchedParticles / nMatchedParticles) = "
      << duplicationRate);
  ACTS_INFO(
      "Average number of duplicated seeds ((nMatchedSeeds - nMatchedParticles) "
      "/ nMatchedParticles) = "
      << aveNDuplicatedSeeds);

  if (m_outputFile) {
    m_outputFile->cd();
    m_effPlotTool.write(m_effPlotCache);
    m_duplicationPlotTool.write(m_duplicationPlotCache);
    ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");
  }
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::SeedingPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const ProtoTrackContainer& tracks) {
  // Read truth information collections
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);

  size_t nSeeds = tracks.size();
  size_t nMatchedSeeds = 0;
  // Map from particles to how many times they were successfully found by a seed
  std::unordered_map<ActsFatras::Barcode, std::size_t> truthCount;

  for (size_t itrack = 0; itrack < tracks.size(); ++itrack) {
    const auto& track = tracks[itrack];
    std::vector<ParticleHitCount> particleHitCounts;
    identifyContributingParticles(hitParticlesMap, track, particleHitCounts);
    // All hits matched to the same particle
    if (particleHitCounts.size() == 1) {
      auto prt = particleHitCounts.at(0);
      auto it = truthCount.try_emplace(prt.particleId, 0u).first;
      it->second += 1;
      nMatchedSeeds++;
    }
  }

  int nMatchedParticles = 0;
  int nDuplicatedParticles = 0;
  // Fill the effeciency and fake rate plots

  std::ofstream myfile;
  myfile.open("eff.csv", std::ios_base::app);

  for (const auto& particle : particles) {
    const auto it1 = truthCount.find(particle.particleId());
    bool isMatched = false;
    int nMatchedSeedsForParticle = 0;
    if (it1 != truthCount.end()) {
      isMatched = true;
      nMatchedParticles++;
      nMatchedSeedsForParticle = it1->second;
      if (nMatchedSeedsForParticle > 1) {
        nDuplicatedParticles++;
      }
    }
    m_effPlotTool.fill(m_effPlotCache, particle, isMatched);
    m_duplicationPlotTool.fill(m_duplicationPlotCache, particle,
                               nMatchedSeedsForParticle - 1);

    const auto t_phi = Acts::VectorHelpers::phi(particle.unitDirection());
    const auto t_eta = Acts::VectorHelpers::eta(particle.unitDirection());
    const auto t_pT = particle.transverseMomentum();

    if (!isMatched) {
      std::cout << "Particle with ID " << particle.particleId() << " is unmatched!" << std::endl;
    }

    myfile << ctx.eventNumber << "," << t_phi << "," << t_eta << "," << t_pT << "," << isMatched << std::endl;
  }

  myfile.close();
  ACTS_DEBUG("Number of seeds: " << nSeeds);
  m_nTotalSeeds += nSeeds;
  m_nTotalMatchedSeeds += nMatchedSeeds;
  m_nTotalParticles += particles.size();
  m_nTotalMatchedParticles += nMatchedParticles;
  m_nTotalDuplicatedParticles += nDuplicatedParticles;

  return ProcessCode::SUCCESS;
}
