// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <mutex>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {

/// @class RootTrajectorySummaryWriter
///
/// Write out the information (including number of measurements, outliers, holes
/// etc. , fitted track parameters and corresponding majority truth particle
/// info) of the reconstructed trajectories into a TTree
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
///
/// Each entry in the TTree corresponds to all reconstructed trajectories in one
/// single event. The event number is part of the written data.
///
/// A common file can be provided for to the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing
/// file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootTrajectorySummaryWriter final
    : public WriterT<TrajectoriesContainer> {
 public:
  struct Config {
    /// Input (fitted) trajectories collection
    std::string inputTrajectories;
    /// Input particles collection.
    std::string inputParticles;
    /// Input hit-particles map collection.
    std::string inputMeasurementParticlesMap;
    /// Output filename.
    std::string filePath = "tracksummary.root";
    /// Name of the output tree.
    std::string treeName = "tracksummary";
    /// File access mode.
    std::string fileMode = "RECREATE";
  };

  /// Constructor
  ///
  /// @param config Configuration struct
  /// @param level Message level declaration
  RootTrajectorySummaryWriter(const Config& config, Acts::Logging::Level level);
  ~RootTrajectorySummaryWriter() final override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  /// @param [in] trajectories are what to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const TrajectoriesContainer& trajectories) final override;

 private:
  Config m_cfg;             ///< The config class
  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  TFile* m_outputFile{nullptr};  ///< The output file
  TTree* m_outputTree{nullptr};  ///< The output tree
  uint32_t m_eventNr{0};         ///< The event number
  std::vector<uint32_t>
      m_multiTrajNr;  ///< The multi-trajectory numbers in event
  std::vector<unsigned int>
      m_subTrajNr;  ///< The multi-trajectory sub-trajectory number in event

  std::vector<unsigned int> m_nStates;        ///< The number of states
  std::vector<unsigned int> m_nMeasurements;  ///< The number of measurements
  std::vector<unsigned int> m_nOutliers;      ///< The number of outliers
  std::vector<unsigned int> m_nHoles;         ///< The number of holes
  std::vector<unsigned int> m_nSharedHits;    ///< The number of shared hits
  std::vector<float> m_chi2Sum;               ///< The total chi2
  std::vector<unsigned int>
      m_NDF;  ///< The number of ndf of the measurements+outliers
  std::vector<std::vector<double>>
      m_measurementChi2;  ///< The chi2 on all measurement states
  std::vector<std::vector<double>>
      m_outlierChi2;  ///< The chi2 on all outlier states
  std::vector<std::vector<double>>
      m_measurementVolume;  ///< The volume id of the measurements
  std::vector<std::vector<double>>
      m_measurementLayer;  ///< The layer id of the measurements
  std::vector<std::vector<double>>
      m_outlierVolume;  ///< The volume id of the outliers
  std::vector<std::vector<double>>
      m_outlierLayer;  ///< The layer id of the outliers

  // The majority truth particle info
  std::vector<unsigned int>
      m_nMajorityHits;  ///< The number of hits from majority particle
  std::vector<uint64_t>
      m_majorityParticleId;      ///< The particle Id of the majority particle
  std::vector<int> m_t_charge;   ///< Charge of majority particle
  std::vector<float> m_t_time;   ///< Time of majority particle
  std::vector<float> m_t_vx;     ///< Vertex x positions of majority particle
  std::vector<float> m_t_vy;     ///< Vertex y positions of majority particle
  std::vector<float> m_t_vz;     ///< Vertex z positions of majority particle
  std::vector<float> m_t_px;     ///< Initial momenta px of majority particle
  std::vector<float> m_t_py;     ///< Initial momenta py of majority particle
  std::vector<float> m_t_pz;     ///< Initial momenta pz of majority particle
  std::vector<float> m_t_theta;  ///< Initial momenta theta of majority particle
  std::vector<float> m_t_phi;    ///< Initial momenta phi of majority particle
  std::vector<float> m_t_p;      ///< Initial abs momenta of majority particle
  std::vector<float> m_t_pT;     ///< Initial momenta pT of majority particle
  std::vector<float> m_t_eta;    ///< Initial momenta eta of majority particle
  std::vector<float>
      m_t_d0;  ///< The extrapolated truth transverse impact parameter
  std::vector<float>
      m_t_z0;  ///< The extrapolated truth longitudinal impact parameter

  std::vector<bool> m_hasFittedParams;  ///< If the track has fitted parameter
  // The fitted parameters
  std::vector<float> m_eLOC0_fit;   ///< Fitted parameters eBoundLoc0 of track
  std::vector<float> m_eLOC1_fit;   ///< Fitted parameters eBoundLoc1 of track
  std::vector<float> m_ePHI_fit;    ///< Fitted parameters ePHI of track
  std::vector<float> m_eTHETA_fit;  ///< Fitted parameters eTHETA of track
  std::vector<float> m_eQOP_fit;    ///< Fitted parameters eQOP of track
  std::vector<float> m_eT_fit;      ///< Fitted parameters eT of track
  // The error of fitted parameters
  std::vector<float> m_err_eLOC0_fit;  ///< Fitted parameters eLOC err of track
  std::vector<float>
      m_err_eLOC1_fit;  ///< Fitted parameters eBoundLoc1 err of track
  std::vector<float> m_err_ePHI_fit;  ///< Fitted parameters ePHI err of track
  std::vector<float>
      m_err_eTHETA_fit;               ///< Fitted parameters eTHETA err of track
  std::vector<float> m_err_eQOP_fit;  ///< Fitted parameters eQOP err of track
  std::vector<float> m_err_eT_fit;    ///< Fitted parameters eT err of track
  // The residual of fitted parameters
  std::vector<float> m_res_eLOC0_fit;  ///< Fitted parameters eLOC res of track
  std::vector<float>
      m_res_eLOC1_fit;  ///< Fitted parameters eBoundLoc1 res of track
  std::vector<float> m_res_ePHI_fit;  ///< Fitted parameters ePHI res of track
  std::vector<float>
      m_res_eTHETA_fit;               ///< Fitted parameters eTHETA res of track
  std::vector<float> m_res_eQOP_fit;  ///< Fitted parameters eQOP res of track
  std::vector<float> m_res_eT_fit;    ///< Fitted parameters eT res of track
  // The pull of fitted parameters
  std::vector<float>
      m_pull_eLOC0_fit;  ///< Fitted parameters eLOC pull of track
  std::vector<float>
      m_pull_eLOC1_fit;  ///< Fitted parameters eBoundLoc1 pull of track
  std::vector<float> m_pull_ePHI_fit;  ///< Fitted parameters ePHI pull of track
  std::vector<float>
      m_pull_eTHETA_fit;  ///< Fitted parameters eTHETA pull of track
  std::vector<float> m_pull_eQOP_fit;  ///< Fitted parameters eQOP pull of track
  std::vector<float> m_pull_eT_fit;    ///< Fitted parameters eT pull of track
};

}  // namespace ActsExamples
