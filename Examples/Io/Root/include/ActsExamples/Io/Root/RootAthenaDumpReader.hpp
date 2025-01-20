// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include <ActsExamples/EventData/Cluster.hpp>
#include <ActsExamples/EventData/SimParticle.hpp>
#include <ActsExamples/EventData/Track.hpp>

#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "TBranch.h"

class TChain;

namespace ActsExamples {

/// @class RootAthenaDumpReader
///
/// @brief Reader for measurements and spacepoints from an Athena
///        object dumper.
///        Specifically written for the input ntuple for GNN
///        See:
///        https://gitlab.cern.ch/atlas/athena/-/blob/main/InnerDetector/InDetGNNTracking/src/DumpObjects.cxx
class RootAthenaDumpReader : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    // Name of tree
    std::string treename;
    // Name of inputfile
    std::vector<std::string> inputfiles;
    // name of the output measurements
    std::string outputMeasurements = "athena_measurements";
    // name of the output pixel space points
    std::string outputPixelSpacePoints = "athena_pixel_spacepoints";
    // name of the output strip space points
    std::string outputStripSpacePoints = "athena_strip_spacepoints";
    // name of the output space points
    std::string outputSpacePoints = "athena_spacepoints";
    // name of the output clusters
    std::string outputClusters = "athena_clusters";
    // name of the output particles
    std::string outputParticles = "athena_particles";
    // name of the measurements -> particles map
    std::string outputMeasurementParticlesMap = "athena_meas_parts_map";
    // name of the particles -> measurements map
    std::string outputParticleMeasurementsMap = "athena_parts_meas_map";
    // name of the track parameters (fitted by athena?)
    std::string outputTrackParameters = "athena_track_parameters";

    /// Only extract spacepoints
    bool onlySpacepoints = false;

    /// Skip truth data
    bool noTruth = false;

    /// Only extract particles that passed the tracking requirements, for
    /// details see:
    /// https://gitlab.cern.ch/atlas/athena/-/blob/main/InnerDetector/InDetGNNTracking/src/DumpObjects.cxx?ref_type=heads#L1363
    bool onlyPassedParticles = false;

    /// Skip spacepoints with phi overlap
    bool skipOverlapSPsPhi = false;

    /// Skip spacepoints with eta overlap
    bool skipOverlapSPsEta = false;

    /// A map that provides a mapping between ACTS and Athena surface
    /// identifiers
    std::shared_ptr<ActsExamples::GeometryIdMapActsAthena> geometryIdMap =
        nullptr;

    /// Tracking Geometry that contains the surfaces where we project
    /// the measurements on
    std::shared_ptr<Acts::TrackingGeometry> trackingGeometry = nullptr;

    /// When projecting measurements on ACTS surfaces, which euclidean boundary
    /// tolerance should be allowed. If a value above zero is needed, this
    /// indicates that the ACTS surfaces do not 100% include the athena surfaces
    double absBoundaryTolerance = 0.0;
  };

  RootAthenaDumpReader(const RootAthenaDumpReader &) = delete;
  RootAthenaDumpReader(const RootAthenaDumpReader &&) = delete;

  // Constructor
  /// @param config The configuration struct
  RootAthenaDumpReader(const Config &config, Acts::Logging::Level level);

  std::string name() const override { return "RootAthenaDumpReader"; }

  /// Return the available events range.
  std::pair<std::size_t, std::size_t> availableEvents() const override {
    return {0u, m_events};
  }

  /// Read out data from the input stream
  ///
  /// @param context The algorithm context
  ProcessCode read(const ActsExamples::AlgorithmContext &ctx) override;

  /// Readonly access to the config
  const Config &config() const { return m_cfg; }

 private:
  /// Particles with barcodes larger then this value are considered to be
  /// secondary particles
  /// https://gitlab.cern.ch/atlas/athena/-/blob/main/InnerDetector/InDetGNNTracking/src/DumpObjects.h?ref_type=heads#L101
  constexpr static int s_maxBarcodeForPrimary = 200000;

  /// Private access to the logging instance
  const Acts::Logger &logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  /// Helper method to read particles
  SimParticleContainer readParticles() const;

  /// Helper method to read measurements
  std::tuple<ClusterContainer, MeasurementContainer,
             IndexMultimap<ActsFatras::Barcode>,
             std::unordered_map<int, std::size_t>>
  readMeasurements(SimParticleContainer &particles,
                   const Acts::GeometryContext &gctx) const;

  /// Helper method to read spacepoints
  /// @param imIdxMap optional remapping of indices. Since the measurement
  /// index must be continuous, we need to remap the measurements indices
  /// if we skip measurements in the first place
  std::tuple<SimSpacePointContainer, SimSpacePointContainer,
             SimSpacePointContainer>
  readSpacepoints(const std::optional<std::unordered_map<int, std::size_t>>
                      &imIdxMap) const;

  /// Helper method to reprocess particle ids
  std::pair<SimParticleContainer, IndexMultimap<ActsFatras::Barcode>>
  reprocessParticles(
      const SimParticleContainer &particles,
      const IndexMultimap<ActsFatras::Barcode> &measPartMap) const;

  /// Write handlers
  WriteDataHandle<SimSpacePointContainer> m_outputPixelSpacePoints{
      this, "outputPixelSpacepoints"};
  WriteDataHandle<SimSpacePointContainer> m_outputStripSpacePoints{
      this, "outputStripSpacepoints"};
  WriteDataHandle<SimSpacePointContainer> m_outputSpacePoints{
      this, "output_spacepoints"};
  WriteDataHandle<ClusterContainer> m_outputClusters{this, "output_clusters"};
  WriteDataHandle<SimParticleContainer> m_outputParticles{this,
                                                          "output_particles"};
  WriteDataHandle<MeasurementContainer> m_outputMeasurements{
      this, "output_measurements"};
  WriteDataHandle<IndexMultimap<ActsFatras::Barcode>> m_outputMeasParticleMap{
      this, "output_meas_part_map"};
  WriteDataHandle<InverseMultimap<ActsFatras::Barcode>> m_outputParticleMeasMap{
      this, "output_part_meas_map"};

  std::unique_ptr<const Acts::Logger> m_logger;
  std::mutex m_read_mutex;

  /// Vector of {eventNr, entryMin, entryMax}
  std::vector<std::tuple<std::uint32_t, std::size_t, std::size_t>> m_eventMap;
  std::shared_ptr<TChain> m_inputchain;
  std::size_t m_events;
  bool m_haveStripFeatures = true;

  static constexpr unsigned int maxCL = 1500000;
  static constexpr unsigned int maxSP = 1500000;
  static constexpr unsigned int maxDTT = 1500000;
  static constexpr unsigned int maxTRK = 1500000;
  static constexpr unsigned int maxPart = 1500000;

  // Declaration of leaf types
  unsigned int run_number = 0;
  ULong64_t event_number = 0;
  int nSE = 0;
  int SEID[4] = {};  //[nSE]
  int nCL = 0;
  int CLindex[maxCL] = {};  //[nCL]

  // Clusters
  std::vector<std::string> *CLhardware{};
  Double_t CLx[maxCL] = {};           //[nCL]
  Double_t CLy[maxCL] = {};           //[nCL]
  Double_t CLz[maxCL] = {};           //[nCL]
  Int_t CLbarrel_endcap[maxCL] = {};  //[nCL]
  Int_t CLlayer_disk[maxCL] = {};     //[nCL]
  Int_t CLeta_module[maxCL] = {};     //[nCL]
  Int_t CLphi_module[maxCL] = {};     //[nCL]
  Int_t CLside[maxCL] = {};           //[nCL]
  ULong64_t CLmoduleID[maxCL] = {};   //[nCL]
  std::vector<std::vector<int>> *CLparticleLink_eventIndex{};
  std::vector<std::vector<int>> *CLparticleLink_barcode{};
  std::vector<std::vector<bool>> *CLbarcodesLinked{};
  std::vector<std::vector<float>> *CLparticle_charge{};
  std::vector<std::vector<int>> *CLphis{};
  std::vector<std::vector<int>> *CLetas{};
  std::vector<std::vector<int>> *CLtots{};
  Double_t CLloc_direction1[maxCL] = {};      //[nCL]
  Double_t CLloc_direction2[maxCL] = {};      //[nCL]
  Double_t CLloc_direction3[maxCL] = {};      //[nCL]
  Double_t CLJan_loc_direction1[maxCL] = {};  //[nCL]
  Double_t CLJan_loc_direction2[maxCL] = {};  //[nCL]
  Double_t CLJan_loc_direction3[maxCL] = {};  //[nCL]
  Int_t CLpixel_count[maxCL] = {};            //[nCL]
  Float_t CLcharge_count[maxCL] = {};         //[nCL]
  Float_t CLloc_eta[maxCL] = {};              //[nCL]
  Float_t CLloc_phi[maxCL] = {};              //[nCL]
  Float_t CLglob_eta[maxCL] = {};             //[nCL]
  Float_t CLglob_phi[maxCL] = {};             //[nCL]
  Double_t CLeta_angle[maxCL] = {};           //[nCL]
  Double_t CLphi_angle[maxCL] = {};           //[nCL]
  Float_t CLnorm_x[maxCL] = {};               //[nCL]
  Float_t CLnorm_y[maxCL] = {};               //[nCL]
  Float_t CLnorm_z[maxCL] = {};               //[nCL]
  std::vector<std::vector<double>> *CLlocal_cov{};

  // Particles
  Int_t nPartEVT = 0;
  Int_t Part_event_number[maxPart] = {};  //[nPartEVT]
  Int_t Part_barcode[maxPart] = {};       //[nPartEVT]
  Float_t Part_px[maxPart] = {};          //[nPartEVT]
  Float_t Part_py[maxPart] = {};          //[nPartEVT]
  Float_t Part_pz[maxPart] = {};          //[nPartEVT]
  Float_t Part_pt[maxPart] = {};          //[nPartEVT]
  Float_t Part_eta[maxPart] = {};         //[nPartEVT]
  Float_t Part_vx[maxPart] = {};          //[nPartEVT]
  Float_t Part_vy[maxPart] = {};          //[nPartEVT]
  Float_t Part_vz[maxPart] = {};          //[nPartEVT]
  Float_t Part_radius[maxPart] = {};      //[nPartEVT]
  Float_t Part_status[maxPart] = {};      //[nPartEVT]
  Float_t Part_charge[maxPart] = {};      //[nPartEVT]
  Int_t Part_pdg_id[maxPart] = {};        //[nPartEVT]
  Int_t Part_passed[maxPart] = {};        //[nPartEVT]
  Int_t Part_vProdNin[maxPart] = {};      //[nPartEVT]
  Int_t Part_vProdNout[maxPart] = {};     //[nPartEVT]
  Int_t Part_vProdStatus[maxPart] = {};   //[nPartEVT]
  Int_t Part_vProdBarcode[maxPart] = {};  //[nPartEVT]
  std::vector<std::vector<int>> *Part_vParentID{};
  std::vector<std::vector<int>> *Part_vParentBarcode{};

  // Spacepoints
  Int_t nSP = 0;
  Int_t SPindex[maxSP] = {};        //[nSP]
  Double_t SPx[maxSP] = {};         //[nSP]
  Double_t SPy[maxSP] = {};         //[nSP]
  Double_t SPz[maxSP] = {};         //[nSP]
  Int_t SPCL1_index[maxSP] = {};    //[nSP]
  Int_t SPCL2_index[maxSP] = {};    //[nSP]
  Int_t SPisOverlap[maxSP] = {};    //[nSP]
  double SPradius[maxSP] = {};      //[nSP]
  double SPcovr[maxSP] = {};        //[nSP]
  double SPcovz[maxSP] = {};        //[nSP]
  float SPhl_topstrip[maxSP] = {};  //[nSP]
  float SPhl_botstrip[maxSP] = {};  //[nSP]
  std::vector<std::vector<float>> *SPtopStripDirection{};
  std::vector<std::vector<float>> *SPbottomStripDirection{};
  std::vector<std::vector<float>> *SPstripCenterDistance{};
  std::vector<std::vector<float>> *SPtopStripCenterPosition{};

  // Those fields are not used currently
  // Keep the code though, since it is annoying to write
  /*
  // Tracks
  Int_t nTRK = 0;
  Int_t TRKindex[maxTRK] = {};                //[nTRK]
  Int_t TRKtrack_fitter[maxTRK] = {};         //[nTRK]
  Int_t TRKparticle_hypothesis[maxTRK] = {};  //[nTRK]
  std::vector<std::vector<int>> *TRKproperties{};
  std::vector<std::vector<int>> *TRKpattern{};
  Int_t TRKndof[maxTRK] = {};     //[nTRK]
  Int_t TRKmot[maxTRK] = {};      //[nTRK]
  Int_t TRKoot[maxTRK] = {};      //[nTRK]
  Float_t TRKchiSq[maxTRK] = {};  //[nTRK]
  std::vector<std::vector<int>> *TRKmeasurementsOnTrack_pixcl_sctcl_index{};
  std::vector<std::vector<int>> *TRKoutliersOnTrack_pixcl_sctcl_index{};
  Int_t TRKcharge[maxTRK] = {};  //[nTRK]
  std::vector<std::vector<double>> *TRKperigee_position{};
  std::vector<std::vector<double>> *TRKperigee_momentum{};
  Int_t TTCindex[maxTRK] = {};          //[nTRK]
  Int_t TTCevent_index[maxTRK] = {};    //[nTRK]
  Int_t TTCparticle_link[maxTRK] = {};  //[nTRK]
  Float_t TTCprobability[maxTRK] = {};  //[nTRK]

  // DDT
  Int_t nDTT = 0;
  Int_t DTTindex[maxDTT] = {};  //[nDTT]
  Int_t DTTsize[maxDTT] = {};   //[nDTT]
  std::vector<std::vector<int>> *DTTtrajectory_eventindex{};
  std::vector<std::vector<int>> *DTTtrajectory_barcode{};
  std::vector<std::vector<int>> *DTTstTruth_subDetType{};
  std::vector<std::vector<int>> *DTTstTrack_subDetType{};
  std::vector<std::vector<int>> *DTTstCommon_subDetType{};
  */
};
}  // namespace ActsExamples
