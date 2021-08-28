// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Performance/TrackFitterPerformanceWriter.hpp"

#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

using Acts::VectorHelpers::eta;

ActsExamples::TrackFitterPerformanceWriter::TrackFitterPerformanceWriter(
    ActsExamples::TrackFitterPerformanceWriter::Config config,
    Acts::Logging::Level level)
    : WriterT(config.inputTrajectories, "TrackFitterPerformanceWriter", level),
      m_cfg(std::move(config)),
      m_resPlotTool(m_cfg.resPlotToolConfig, level),
      m_effPlotTool(m_cfg.effPlotToolConfig, level),
      m_trackSummaryPlotTool(m_cfg.trackSummaryPlotToolConfig, level)

{
  // trajectories collection name is already checked by base ctor
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing hit-particles map input collection");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  // the output file can not be given externally since TFile accesses to the
  // same file from multiple threads are unsafe.
  // must always be opened internally
  auto path = m_cfg.filePath;
  m_outputFile = TFile::Open(path.c_str(), "RECREATE");
  if (not m_outputFile) {
    throw std::invalid_argument("Could not open '" + path + "'");
  }

  // initialize the residual and efficiency plots tool
  m_resPlotTool.book(m_resPlotCache);
  m_effPlotTool.book(m_effPlotCache);
  m_trackSummaryPlotTool.book(m_trackSummaryPlotCache);
}

ActsExamples::TrackFitterPerformanceWriter::~TrackFitterPerformanceWriter() {
  m_resPlotTool.clear(m_resPlotCache);
  m_effPlotTool.clear(m_effPlotCache);
  m_trackSummaryPlotTool.clear(m_trackSummaryPlotCache);

  if (m_outputFile) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::TrackFitterPerformanceWriter::endRun() {
  // fill residual and pull details into additional hists
  m_resPlotTool.refinement(m_resPlotCache);

  if (m_outputFile) {
    m_outputFile->cd();
    m_resPlotTool.write(m_resPlotCache);
    m_effPlotTool.write(m_effPlotCache);
    m_trackSummaryPlotTool.write(m_trackSummaryPlotCache);

    ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");
  }
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::TrackFitterPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const TrajectoriesContainer& trajectories) {
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;

  // Read truth input collections
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);

  // Truth particles with corresponding reconstructed tracks
  std::vector<ActsFatras::Barcode> reconParticleIds;
  reconParticleIds.reserve(particles.size());
  // For each particle within a track, how many hits did it contribute
  std::vector<ParticleHitCount> particleHitCounts;

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Loop over all trajectories
  for (size_t itraj = 0; itraj < trajectories.size(); ++itraj) {
    const auto& traj = trajectories[itraj];

    if (traj.empty()) {
      ACTS_WARNING("Empty trajectories object " << itraj);
      continue;
    }

    // The trajectory entry indices and the multiTrajectory
    const auto& trackTips = traj.tips();
    const auto& mj = traj.multiTrajectory();

    // Check the size of the trajectory entry indices. For track fitting, there
    // should be at most one trajectory
    if (trackTips.size() > 1) {
      ACTS_ERROR("Track fitting should not result in multiple trajectories.");
      return ProcessCode::ABORT;
    }
    // Get the entry index for the single trajectory
    auto trackTip = trackTips.front();

    // Select reco track with fitted parameters
    if (not traj.hasTrackParameters(trackTip)) {
      ACTS_WARNING("No fitted track parameters.");
      continue;
    }
    const auto& fittedParameters = traj.trackParameters(trackTip);

    // Get the majority truth particle for this trajectory
    identifyContributingParticles(hitParticlesMap, traj, trackTip,
                                  particleHitCounts);
    if (particleHitCounts.empty()) {
      ACTS_WARNING("No truth particle associated with this trajectory.");
      continue;
    }
    // Find the truth particle for the majority barcode
    const auto ip = particles.find(particleHitCounts.front().particleId);
    if (ip == particles.end()) {
      ACTS_WARNING("Majority particle not found in the particles collection.");
      continue;
    }

    // Record this majority particle ID of this trajectory
    reconParticleIds.push_back(ip->particleId());
    // Fill the residual plots
    m_resPlotTool.fill(m_resPlotCache, ctx.geoContext, *ip,
                       traj.trackParameters(trackTip));
    // Collect the trajectory summary info
    auto trajState =
        Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
    // Fill the trajectory summary info
    m_trackSummaryPlotTool.fill(m_trackSummaryPlotCache, fittedParameters,
                                trajState.nStates, trajState.nMeasurements,
                                trajState.nOutliers, trajState.nHoles,
                                trajState.nSharedHits);
  }

  // Fill the efficiency, defined as the ratio between number of tracks with
  // fitted parameter and total truth tracks (assumes one truth partilce has
  // one truth track)
  for (const auto& particle : particles) {
    bool isReconstructed = false;
    // Find if the particle has been reconstructed
    auto it = std::find(reconParticleIds.begin(), reconParticleIds.end(),
                        particle.particleId());
    if (it != reconParticleIds.end()) {
      isReconstructed = true;
    }
    m_effPlotTool.fill(m_effPlotCache, particle, isReconstructed);
  }

  return ProcessCode::SUCCESS;
}
