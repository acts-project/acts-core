// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootTrajectorySummaryWriter.hpp"

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <ios>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

ActsExamples::RootTrajectorySummaryWriter::RootTrajectorySummaryWriter(
    const ActsExamples::RootTrajectorySummaryWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT(config.inputTrajectories, "RootTrajectorySummaryWriter", level),
      m_cfg(config) {
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
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Setup ROOT I/O
  auto path = m_cfg.filePath;
  m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + path);
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr)
    throw std::bad_alloc();
  else {
    // I/O parameters
    m_outputTree->Branch("event_nr", &m_eventNr);
    m_outputTree->Branch("multiTraj_nr", &m_multiTrajNr);
    m_outputTree->Branch("subTraj_nr", &m_subTrajNr);

    m_outputTree->Branch("nStates", &m_nStates);
    m_outputTree->Branch("nMeasurements", &m_nMeasurements);
    m_outputTree->Branch("nOutliers", &m_nOutliers);
    m_outputTree->Branch("nHoles", &m_nHoles);
    m_outputTree->Branch("nSharedHits", &m_nSharedHits);
    m_outputTree->Branch("chi2Sum", &m_chi2Sum);
    m_outputTree->Branch("NDF", &m_NDF);
    m_outputTree->Branch("measurementChi2", &m_measurementChi2);
    m_outputTree->Branch("outlierChi2", &m_outlierChi2);
    m_outputTree->Branch("measurementVolume", &m_measurementVolume);
    m_outputTree->Branch("measurementLayer", &m_measurementLayer);
    m_outputTree->Branch("outlierVolume", &m_outlierVolume);
    m_outputTree->Branch("outlierLayer", &m_outlierLayer);

    m_outputTree->Branch("nMajorityHits", &m_nMajorityHits);
    m_outputTree->Branch("majorityParticleId", &m_majorityParticleId);
    m_outputTree->Branch("t_charge", &m_t_charge);
    m_outputTree->Branch("t_time", &m_t_time);
    m_outputTree->Branch("t_vx", &m_t_vx);
    m_outputTree->Branch("t_vy", &m_t_vy);
    m_outputTree->Branch("t_vz", &m_t_vz);
    m_outputTree->Branch("t_px", &m_t_px);
    m_outputTree->Branch("t_py", &m_t_py);
    m_outputTree->Branch("t_pz", &m_t_pz);
    m_outputTree->Branch("t_theta", &m_t_theta);
    m_outputTree->Branch("t_phi", &m_t_phi);
    m_outputTree->Branch("t_eta", &m_t_eta);
    m_outputTree->Branch("t_p", &m_t_p);
    m_outputTree->Branch("t_pT", &m_t_pT);
    m_outputTree->Branch("t_d0", &m_t_d0);
    m_outputTree->Branch("t_z0", &m_t_z0);

    m_outputTree->Branch("hasFittedParams", &m_hasFittedParams);
    m_outputTree->Branch("eLOC0_fit", &m_eLOC0_fit);
    m_outputTree->Branch("eLOC1_fit", &m_eLOC1_fit);
    m_outputTree->Branch("ePHI_fit", &m_ePHI_fit);
    m_outputTree->Branch("eTHETA_fit", &m_eTHETA_fit);
    m_outputTree->Branch("eQOP_fit", &m_eQOP_fit);
    m_outputTree->Branch("eT_fit", &m_eT_fit);
    m_outputTree->Branch("err_eLOC0_fit", &m_err_eLOC0_fit);
    m_outputTree->Branch("err_eLOC1_fit", &m_err_eLOC1_fit);
    m_outputTree->Branch("err_ePHI_fit", &m_err_ePHI_fit);
    m_outputTree->Branch("err_eTHETA_fit", &m_err_eTHETA_fit);
    m_outputTree->Branch("err_eQOP_fit", &m_err_eQOP_fit);
    m_outputTree->Branch("err_eT_fit", &m_err_eT_fit);
    m_outputTree->Branch("res_eLOC0_fit", &m_res_eLOC0_fit);
    m_outputTree->Branch("res_eLOC1_fit", &m_res_eLOC1_fit);
    m_outputTree->Branch("res_ePHI_fit", &m_res_ePHI_fit);
    m_outputTree->Branch("res_eTHETA_fit", &m_res_eTHETA_fit);
    m_outputTree->Branch("res_eQOP_fit", &m_res_eQOP_fit);
    m_outputTree->Branch("res_eT_fit", &m_res_eT_fit);
    m_outputTree->Branch("pull_eLOC0_fit", &m_pull_eLOC0_fit);
    m_outputTree->Branch("pull_eLOC1_fit", &m_pull_eLOC1_fit);
    m_outputTree->Branch("pull_ePHI_fit", &m_pull_ePHI_fit);
    m_outputTree->Branch("pull_eTHETA_fit", &m_pull_eTHETA_fit);
    m_outputTree->Branch("pull_eQOP_fit", &m_pull_eQOP_fit);
    m_outputTree->Branch("pull_eT_fit", &m_pull_eT_fit);
  }
}

ActsExamples::RootTrajectorySummaryWriter::~RootTrajectorySummaryWriter() {}

ActsExamples::ProcessCode ActsExamples::RootTrajectorySummaryWriter::endRun() {
  if (m_outputFile) {
    m_outputFile->cd();
    m_outputTree->Write();
    ACTS_INFO("Write parameters of trajectories to tree '"
              << m_cfg.treeName << "' in '" << m_cfg.filePath << "'");
    m_outputFile->Close();
  }
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootTrajectorySummaryWriter::writeT(
    const AlgorithmContext& ctx, const TrajectoriesContainer& trajectories) {
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;

  if (m_outputFile == nullptr)
    return ProcessCode::SUCCESS;

  // Read additional input collections
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);

  // For each particle within a track, how many hits did it contribute
  std::vector<ParticleHitCount> particleHitCounts;

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Get the event number
  m_eventNr = ctx.eventNumber;

  // Loop over the trajectories
  for (size_t itraj = 0; itraj < trajectories.size(); ++itraj) {
    const auto& traj = trajectories[itraj];

    if (traj.empty()) {
      ACTS_WARNING("Empty trajectories object " << itraj);
      continue;
    }

    // The trajectory index
    m_multiTrajNr.push_back(itraj);

    // The trajectory entry indices and the multiTrajectory
    const auto& mj = traj.multiTrajectory();
    const auto& trackTips = traj.tips();

    // Loop over the entry indices for the subtrajectories
    for (unsigned int isubtraj = 0; isubtraj < trackTips.size(); ++isubtraj) {
      // The subtrajectory index
      m_subTrajNr.push_back(isubtraj);
      // The entry index for this subtrajectory
      const auto& trackTip = trackTips[isubtraj];

      // Collect the trajectory summary info
      auto trajState =
          Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
      m_nStates.push_back(trajState.nStates);
      m_nMeasurements.push_back(trajState.nMeasurements);
      m_nOutliers.push_back(trajState.nOutliers);
      m_nHoles.push_back(trajState.nHoles);
      m_nSharedHits.push_back(trajState.nSharedHits);
      m_chi2Sum.push_back(trajState.chi2Sum);
      m_NDF.push_back(trajState.NDF);
      m_measurementChi2.push_back(trajState.measurementChi2);
      m_outlierChi2.push_back(trajState.measurementChi2);
      // They are stored as double (as the vector of vector of int is not known
      // to ROOT)
      m_measurementVolume.emplace_back(trajState.measurementVolume.begin(),
                                       trajState.measurementVolume.end());
      m_measurementLayer.emplace_back(trajState.measurementLayer.begin(),
                                      trajState.measurementLayer.end());
      m_outlierVolume.emplace_back(trajState.outlierVolume.begin(),
                                   trajState.outlierVolume.end());
      m_outlierLayer.emplace_back(trajState.outlierLayer.begin(),
                                  trajState.outlierLayer.end());

      // Initialize the truth particle info
      uint64_t majorityParticleId = NaNint;
      unsigned int nMajorityHits = NaNint;
      float t_charge = NaNint;
      float t_time = NaNfloat;
      float t_vx = NaNfloat;
      float t_vy = NaNfloat;
      float t_vz = NaNfloat;
      float t_px = NaNfloat;
      float t_py = NaNfloat;
      float t_pz = NaNfloat;
      float t_theta = NaNfloat;
      float t_phi = NaNfloat;
      float t_eta = NaNfloat;
      float t_p = NaNfloat;
      float t_pT = NaNfloat;
      float t_d0 = NaNfloat;
      float t_z0 = NaNfloat;

      // Get the perigee surface
      Acts::Surface* pSurface = nullptr;
      if (traj.hasTrackParameters(trackTip)) {
        const auto& boundParam = traj.trackParameters(trackTip);
        pSurface = const_cast<Acts::Surface*>(&boundParam.referenceSurface());
      }

      // Get the majority truth particle to this track
      identifyContributingParticles(hitParticlesMap, traj, trackTip,
                                    particleHitCounts);
      bool foundMajorityParticle = false;
      // Get the truth particle info
      if (not particleHitCounts.empty()) {
        // Get the barcode of the majority truth particle
        majorityParticleId = particleHitCounts.front().particleId.value();
        nMajorityHits = particleHitCounts.front().hitCount;

        // Find the truth particle via the barcode
        auto ip = particles.find(majorityParticleId);
        if (ip != particles.end()) {
          foundMajorityParticle = true;

          const auto& particle = *ip;
          ACTS_DEBUG(
              "Find the truth particle with barcode = " << majorityParticleId);
          // Get the truth particle info at vertex
          t_p = particle.absoluteMomentum();
          t_charge = particle.charge();
          t_time = particle.time();
          t_vx = particle.position().x();
          t_vy = particle.position().y();
          t_vz = particle.position().z();
          t_px = t_p * particle.unitDirection().x();
          t_py = t_p * particle.unitDirection().y();
          t_pz = t_p * particle.unitDirection().z();
          t_theta = theta(particle.unitDirection());
          t_phi = phi(particle.unitDirection());
          t_eta = eta(particle.unitDirection());
          t_pT = t_p * perp(particle.unitDirection());

          if (pSurface) {
            // get the truth perigee parameter
            auto lpResult = pSurface->globalToLocal(
                ctx.geoContext, particle.position(), particle.unitDirection());
            if (lpResult.ok()) {
              t_d0 = lpResult.value()[Acts::BoundIndices::eBoundLoc0];
              t_z0 = lpResult.value()[Acts::BoundIndices::eBoundLoc1];
            } else {
              ACTS_ERROR("Global to local transformation did not succeed.");
            }
          }
        } else {
          ACTS_WARNING("Truth particle with barcode = "
                       << majorityParticleId
                       << " not found in the input collection!");
        }
      }
      if (not foundMajorityParticle) {
        ACTS_WARNING("Truth particle for mj " << itraj << " subtraj "
                                              << isubtraj << " not found!");
      }

      // Push the corresponding truth particle info for the track.
      // Always push back even if majority particle not found
      m_majorityParticleId.push_back(majorityParticleId);
      m_nMajorityHits.push_back(nMajorityHits);
      m_t_charge.push_back(t_charge);
      m_t_time.push_back(t_time);
      m_t_vx.push_back(t_vx);
      m_t_vy.push_back(t_vy);
      m_t_vz.push_back(t_vz);
      m_t_px.push_back(t_px);
      m_t_py.push_back(t_py);
      m_t_pz.push_back(t_pz);
      m_t_theta.push_back(t_theta);
      m_t_phi.push_back(t_phi);
      m_t_eta.push_back(t_eta);
      m_t_p.push_back(t_p);
      m_t_pT.push_back(t_pT);
      m_t_d0.push_back(t_d0);
      m_t_z0.push_back(t_z0);

      // Initialize the fitted track parameters info
      std::array<float, Acts::eBoundSize> param = {
          NaNfloat, NaNfloat, NaNfloat, NaNfloat, NaNfloat, NaNfloat};
      std::array<float, Acts::eBoundSize> res = {NaNfloat, NaNfloat, NaNfloat,
                                                 NaNfloat, NaNfloat, NaNfloat};
      std::array<float, Acts::eBoundSize> error = {
          NaNfloat, NaNfloat, NaNfloat, NaNfloat, NaNfloat, NaNfloat};
      std::array<float, Acts::eBoundSize> pull = {NaNfloat, NaNfloat, NaNfloat,
                                                  NaNfloat, NaNfloat, NaNfloat};
      bool hasFittedParams = false;
      if (traj.hasTrackParameters(trackTip)) {
        hasFittedParams = true;
        const auto& boundParam = traj.trackParameters(trackTip);
        const auto& parameter = boundParam.parameters();
        for (unsigned int i = 0; i < Acts::eBoundSize; ++i) {
          param[i] = parameter[i];
        }

        res = {param[Acts::eBoundLoc0] - t_d0,
               param[Acts::eBoundLoc1] - t_z0,
               param[Acts::eBoundPhi] - t_phi,
               param[Acts::eBoundTheta] - t_theta,
               param[Acts::eBoundQOverP] - t_charge / t_p,
               param[Acts::eBoundTime] - t_time};

        if (boundParam.covariance().has_value()) {
          const auto& covariance = *boundParam.covariance();
          for (unsigned int i = 0; i < Acts::eBoundSize; ++i) {
            error[i] = std::sqrt(covariance(i, i));
            pull[i] = res[i] / error[i];
          }
        }
      }

      // Push the fitted track parameters.
      // Always push back even if no fitted track parameters
      m_eLOC0_fit.push_back(param[Acts::eBoundLoc0]);
      m_eLOC1_fit.push_back(param[Acts::eBoundLoc1]);
      m_ePHI_fit.push_back(param[Acts::eBoundPhi]);
      m_eTHETA_fit.push_back(param[Acts::eBoundTheta]);
      m_eQOP_fit.push_back(param[Acts::eBoundQOverP]);
      m_eT_fit.push_back(param[Acts::eBoundTime]);

      m_res_eLOC0_fit.push_back(res[Acts::eBoundLoc0]);
      m_res_eLOC1_fit.push_back(res[Acts::eBoundLoc1]);
      m_res_ePHI_fit.push_back(res[Acts::eBoundPhi]);
      m_res_eTHETA_fit.push_back(res[Acts::eBoundTheta]);
      m_res_eQOP_fit.push_back(res[Acts::eBoundQOverP]);
      m_res_eT_fit.push_back(res[Acts::eBoundTime]);

      m_err_eLOC0_fit.push_back(error[Acts::eBoundLoc0]);
      m_err_eLOC1_fit.push_back(error[Acts::eBoundLoc1]);
      m_err_ePHI_fit.push_back(error[Acts::eBoundPhi]);
      m_err_eTHETA_fit.push_back(error[Acts::eBoundTheta]);
      m_err_eQOP_fit.push_back(error[Acts::eBoundQOverP]);
      m_err_eT_fit.push_back(error[Acts::eBoundTime]);

      m_pull_eLOC0_fit.push_back(pull[Acts::eBoundLoc0]);
      m_pull_eLOC1_fit.push_back(pull[Acts::eBoundLoc1]);
      m_pull_ePHI_fit.push_back(pull[Acts::eBoundPhi]);
      m_pull_eTHETA_fit.push_back(pull[Acts::eBoundTheta]);
      m_pull_eQOP_fit.push_back(pull[Acts::eBoundQOverP]);
      m_pull_eT_fit.push_back(pull[Acts::eBoundTime]);

      m_hasFittedParams.push_back(hasFittedParams);
    }  // all subtrajectories
  }    // all trajectories

  // fill the variables
  m_outputTree->Fill();

  m_multiTrajNr.clear();
  m_subTrajNr.clear();
  m_nStates.clear();
  m_nMeasurements.clear();
  m_nOutliers.clear();
  m_nHoles.clear();
  m_nSharedHits.clear();
  m_chi2Sum.clear();
  m_NDF.clear();
  m_measurementChi2.clear();
  m_outlierChi2.clear();
  m_measurementVolume.clear();
  m_measurementLayer.clear();
  m_outlierVolume.clear();
  m_outlierLayer.clear();

  m_nMajorityHits.clear();
  m_majorityParticleId.clear();
  m_t_charge.clear();
  m_t_time.clear();
  m_t_vx.clear();
  m_t_vy.clear();
  m_t_vz.clear();
  m_t_px.clear();
  m_t_py.clear();
  m_t_pz.clear();
  m_t_theta.clear();
  m_t_phi.clear();
  m_t_p.clear();
  m_t_pT.clear();
  m_t_eta.clear();
  m_t_d0.clear();
  m_t_z0.clear();

  m_hasFittedParams.clear();
  m_eLOC0_fit.clear();
  m_eLOC1_fit.clear();
  m_ePHI_fit.clear();
  m_eTHETA_fit.clear();
  m_eQOP_fit.clear();
  m_eT_fit.clear();
  m_err_eLOC0_fit.clear();
  m_err_eLOC1_fit.clear();
  m_err_ePHI_fit.clear();
  m_err_eTHETA_fit.clear();
  m_err_eQOP_fit.clear();
  m_err_eT_fit.clear();
  m_res_eLOC0_fit.clear();
  m_res_eLOC1_fit.clear();
  m_res_ePHI_fit.clear();
  m_res_eTHETA_fit.clear();
  m_res_eQOP_fit.clear();
  m_res_eT_fit.clear();
  m_pull_eLOC0_fit.clear();
  m_pull_eLOC1_fit.clear();
  m_pull_ePHI_fit.clear();
  m_pull_eTHETA_fit.clear();
  m_pull_eQOP_fit.clear();
  m_pull_eT_fit.clear();

  return ProcessCode::SUCCESS;
}
