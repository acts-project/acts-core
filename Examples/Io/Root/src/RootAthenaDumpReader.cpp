// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootAthenaDumpReader.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include <ActsExamples/Digitization/MeasurementCreation.hpp>

#include <TChain.h>
#include <boost/container/static_vector.hpp>

using namespace Acts::UnitLiterals;

namespace {
std::uint64_t concatInts(int a, int b) {
  auto va = static_cast<std::uint32_t>(a);
  auto vb = static_cast<std::uint32_t>(b);
  std::uint64_t value = (static_cast<std::uint64_t>(va) << 32) | vb;
  return value;
}

std::pair<std::uint32_t, std::uint32_t> splitInt(std::uint64_t v) {
  return {static_cast<std::uint32_t>((v & 0xFFFFFFFF00000000LL) >> 32),
          static_cast<std::uint32_t>(v & 0xFFFFFFFFLL)};
}

/// In cases when there is built up a particle collection in an iterative way it
/// can be way faster to build up a vector and afterwards use a special
/// constructor to speed up the set creation.
inline auto particleVectorToSet(
    std::vector<ActsExamples::SimParticle>& particles) {
  using namespace ActsExamples;
  auto cmp = [](const auto& a, const auto& b) {
    return a.particleId().value() == b.particleId().value();
  };

  std::sort(particles.begin(), particles.end(), detail::CompareParticleId{});
  particles.erase(std::unique(particles.begin(), particles.end(), cmp),
                  particles.end());

  return SimParticleContainer(boost::container::ordered_unique_range_t{},
                              particles.begin(), particles.end());
}

}  // namespace

enum SpacePointType { ePixel = 1, eStrip = 2 };

namespace ActsExamples {

RootAthenaDumpReader::RootAthenaDumpReader(
    const RootAthenaDumpReader::Config& config, Acts::Logging::Level level)
    : IReader(),
      m_cfg(config),
      m_logger(Acts::getDefaultLogger(name(), level)) {
  if (m_cfg.inputfiles.empty()) {
    throw std::invalid_argument("Empty input file list");
  }
  if (m_cfg.treename.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_inputchain = std::make_shared<TChain>(m_cfg.treename.c_str());

  m_outputPixelSpacePoints.initialize(m_cfg.outputPixelSpacePoints);
  m_outputStripSpacePoints.initialize(m_cfg.outputStripSpacePoints);
  m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);
  if (!m_cfg.onlySpacepoints) {
    m_outputMeasurements.initialize(m_cfg.outputMeasurements);
    if (!m_cfg.noTruth) {
      m_outputClusters.initialize(m_cfg.outputClusters);
      m_outputParticles.initialize(m_cfg.outputParticles);
      m_outputMeasParticleMap.initialize(m_cfg.outputMeasurementParticlesMap);
      m_outputParticleMeasMap.initialize(m_cfg.outputParticleMeasurementsMap);
    }
  }

  if (m_inputchain->GetBranch("SPtopStripDirection") == nullptr) {
    ACTS_WARNING("Additional SP strip features not available");
    m_haveStripFeatures = false;
  }

  // Set the branches
  m_inputchain->SetBranchAddress("run_number", &run_number);
  m_inputchain->SetBranchAddress("event_number", &event_number);
  m_inputchain->SetBranchAddress("nSE", &nSE);
  m_inputchain->SetBranchAddress("SEID", SEID);
  m_inputchain->SetBranchAddress("nCL", &nCL);
  m_inputchain->SetBranchAddress("CLindex", CLindex);
  m_inputchain->SetBranchAddress("CLhardware", &CLhardware);
  m_inputchain->SetBranchAddress("CLx", CLx);
  m_inputchain->SetBranchAddress("CLy", CLy);
  m_inputchain->SetBranchAddress("CLz", CLz);
  m_inputchain->SetBranchAddress("CLbarrel_endcap", CLbarrel_endcap);
  m_inputchain->SetBranchAddress("CLlayer_disk", CLlayer_disk);
  m_inputchain->SetBranchAddress("CLeta_module", CLeta_module);
  m_inputchain->SetBranchAddress("CLphi_module", CLphi_module);
  m_inputchain->SetBranchAddress("CLside", CLside);
  m_inputchain->SetBranchAddress("CLmoduleID", CLmoduleID);
  m_inputchain->SetBranchAddress("CLphis", &CLphis);
  m_inputchain->SetBranchAddress("CLetas", &CLetas);
  m_inputchain->SetBranchAddress("CLtots", &CLtots);
  m_inputchain->SetBranchAddress("CLloc_direction1", CLloc_direction1);
  m_inputchain->SetBranchAddress("CLloc_direction2", CLloc_direction2);
  m_inputchain->SetBranchAddress("CLloc_direction3", CLloc_direction3);
  m_inputchain->SetBranchAddress("CLJan_loc_direction1", CLJan_loc_direction1);
  m_inputchain->SetBranchAddress("CLJan_loc_direction2", CLJan_loc_direction2);
  m_inputchain->SetBranchAddress("CLJan_loc_direction3", CLJan_loc_direction3);
  m_inputchain->SetBranchAddress("CLpixel_count", CLpixel_count);
  m_inputchain->SetBranchAddress("CLcharge_count", CLcharge_count);
  m_inputchain->SetBranchAddress("CLloc_eta", CLloc_eta);
  m_inputchain->SetBranchAddress("CLloc_phi", CLloc_phi);
  m_inputchain->SetBranchAddress("CLglob_eta", CLglob_eta);
  m_inputchain->SetBranchAddress("CLglob_phi", CLglob_phi);
  m_inputchain->SetBranchAddress("CLeta_angle", CLeta_angle);
  m_inputchain->SetBranchAddress("CLphi_angle", CLphi_angle);
  m_inputchain->SetBranchAddress("CLnorm_x", CLnorm_x);
  m_inputchain->SetBranchAddress("CLnorm_y", CLnorm_y);
  m_inputchain->SetBranchAddress("CLnorm_z", CLnorm_z);
  m_inputchain->SetBranchAddress("CLlocal_cov", &CLlocal_cov);
  if (!m_cfg.noTruth) {
    m_inputchain->SetBranchAddress("CLparticleLink_eventIndex",
                                   &CLparticleLink_eventIndex);
    m_inputchain->SetBranchAddress("CLparticleLink_barcode",
                                   &CLparticleLink_barcode);
    m_inputchain->SetBranchAddress("CLbarcodesLinked", &CLbarcodesLinked);
    m_inputchain->SetBranchAddress("CLparticle_charge", &CLparticle_charge);
  }

  if (!m_cfg.noTruth) {
    m_inputchain->SetBranchAddress("nPartEVT", &nPartEVT);
    m_inputchain->SetBranchAddress("Part_event_number", Part_event_number);
    m_inputchain->SetBranchAddress("Part_barcode", Part_barcode);
    m_inputchain->SetBranchAddress("Part_px", Part_px);
    m_inputchain->SetBranchAddress("Part_py", Part_py);
    m_inputchain->SetBranchAddress("Part_pz", Part_pz);
    m_inputchain->SetBranchAddress("Part_pt", Part_pt);
    m_inputchain->SetBranchAddress("Part_eta", Part_eta);
    m_inputchain->SetBranchAddress("Part_vx", Part_vx);
    m_inputchain->SetBranchAddress("Part_vy", Part_vy);
    m_inputchain->SetBranchAddress("Part_vz", Part_vz);
    m_inputchain->SetBranchAddress("Part_radius", Part_radius);
    m_inputchain->SetBranchAddress("Part_status", Part_status);
    m_inputchain->SetBranchAddress("Part_charge", Part_charge);
    m_inputchain->SetBranchAddress("Part_pdg_id", Part_pdg_id);
    m_inputchain->SetBranchAddress("Part_passed", Part_passed);
    m_inputchain->SetBranchAddress("Part_vProdNin", Part_vProdNin);
    m_inputchain->SetBranchAddress("Part_vProdNout", Part_vProdNout);
    m_inputchain->SetBranchAddress("Part_vProdStatus", Part_vProdStatus);
    m_inputchain->SetBranchAddress("Part_vProdBarcode", Part_vProdBarcode);
    m_inputchain->SetBranchAddress("Part_vParentID", &Part_vParentID);
    m_inputchain->SetBranchAddress("Part_vParentBarcode", &Part_vParentBarcode);
  }

  m_inputchain->SetBranchAddress("nSP", &nSP);
  m_inputchain->SetBranchAddress("SPindex", SPindex);
  m_inputchain->SetBranchAddress("SPx", SPx);
  m_inputchain->SetBranchAddress("SPy", SPy);
  m_inputchain->SetBranchAddress("SPz", SPz);
  m_inputchain->SetBranchAddress("SPCL1_index", SPCL1_index);
  m_inputchain->SetBranchAddress("SPCL2_index", SPCL2_index);
  m_inputchain->SetBranchAddress("SPisOverlap", SPisOverlap);

  if (m_haveStripFeatures) {
    m_inputchain->SetBranchAddress("SPradius", SPradius);
    m_inputchain->SetBranchAddress("SPcovr", SPcovr);
    m_inputchain->SetBranchAddress("SPcovz", SPcovz);
    m_inputchain->SetBranchAddress("SPhl_topstrip", SPhl_topstrip);
    m_inputchain->SetBranchAddress("SPhl_botstrip", SPhl_botstrip);
    m_inputchain->SetBranchAddress("SPtopStripDirection", &SPtopStripDirection);
    m_inputchain->SetBranchAddress("SPbottomStripDirection",
                                   &SPbottomStripDirection);
    m_inputchain->SetBranchAddress("SPstripCenterDistance",
                                   &SPstripCenterDistance);
    m_inputchain->SetBranchAddress("SPtopStripCenterPosition",
                                   &SPtopStripCenterPosition);
  }

  // These quantities are not used currently and thus commented out
  // I would like to keep the code, since it is always a pain to write it
  /*
  m_inputchain->SetBranchAddress("nTRK", &nTRK);
  m_inputchain->SetBranchAddress("TRKindex", TRKindex);
  m_inputchain->SetBranchAddress("TRKtrack_fitter", TRKtrack_fitter);
  m_inputchain->SetBranchAddress("TRKparticle_hypothesis",
                                 TRKparticle_hypothesis);
  m_inputchain->SetBranchAddress("TRKproperties", &TRKproperties);
  m_inputchain->SetBranchAddress("TRKpattern", &TRKpattern);
  m_inputchain->SetBranchAddress("TRKndof", TRKndof);
  m_inputchain->SetBranchAddress("TRKmot", TRKmot);
  m_inputchain->SetBranchAddress("TRKoot", TRKoot);
  m_inputchain->SetBranchAddress("TRKchiSq", TRKchiSq);
  m_inputchain->SetBranchAddress("TRKmeasurementsOnTrack_pixcl_sctcl_index",
                                 &TRKmeasurementsOnTrack_pixcl_sctcl_index);
  m_inputchain->SetBranchAddress("TRKoutliersOnTrack_pixcl_sctcl_index",
                                 &TRKoutliersOnTrack_pixcl_sctcl_index);
  m_inputchain->SetBranchAddress("TRKcharge", TRKcharge);
  m_inputchain->SetBranchAddress("TRKperigee_position", &TRKperigee_position);
  m_inputchain->SetBranchAddress("TRKperigee_momentum", &TRKperigee_momentum);
  m_inputchain->SetBranchAddress("TTCindex", TTCindex);
  m_inputchain->SetBranchAddress("TTCevent_index", TTCevent_index);
  m_inputchain->SetBranchAddress("TTCparticle_link", TTCparticle_link);
  m_inputchain->SetBranchAddress("TTCprobability", TTCprobability);
  m_inputchain->SetBranchAddress("nDTT", &nDTT);
  m_inputchain->SetBranchAddress("DTTindex", DTTindex);
  m_inputchain->SetBranchAddress("DTTsize", DTTsize);
  m_inputchain->SetBranchAddress("DTTtrajectory_eventindex",
                                 &DTTtrajectory_eventindex);
  m_inputchain->SetBranchAddress("DTTtrajectory_barcode",
                                 &DTTtrajectory_barcode);
  m_inputchain->SetBranchAddress("DTTstTruth_subDetType",
                                 &DTTstTruth_subDetType);
  m_inputchain->SetBranchAddress("DTTstTrack_subDetType",
                                 &DTTstTrack_subDetType);
  m_inputchain->SetBranchAddress("DTTstCommon_subDetType",
                                 &DTTstCommon_subDetType);
  */

  for (const auto& file : m_cfg.inputfiles) {
    m_inputchain->Add(file.c_str());
    ACTS_DEBUG("Adding file '" << file << "' to tree " << m_cfg.treename);
  }

  m_events = m_inputchain->GetEntries();

  ACTS_DEBUG("End of constructor. In total available events=" << m_events);
}  // constructor

SimParticleContainer RootAthenaDumpReader::readParticles() const {
  std::vector<SimParticle> particles;
  particles.reserve(nPartEVT);

  for (auto ip = 0; ip < nPartEVT; ++ip) {
    if (m_cfg.onlyPassedParticles && !static_cast<bool>(Part_passed[ip])) {
      continue;
    }

    auto dummyBarcode = concatInts(Part_barcode[ip], Part_event_number[ip]);
    SimParticleState particle(dummyBarcode,
                              static_cast<Acts::PdgParticle>(Part_pdg_id[ip]));

    Acts::Vector3 p = Acts::Vector3{Part_px[ip], Part_py[ip], Part_pz[ip]} *
                      Acts::UnitConstants::MeV;
    particle.setAbsoluteMomentum(p.norm());

    particle.setDirection(p.normalized());

    auto x = Acts::Vector4{Part_vx[ip], Part_vy[ip], Part_vz[ip], 0.0};
    particle.setPosition4(x);

    particles.push_back(SimParticle(particle, particle));
  }

  ACTS_DEBUG("Created " << particles.size() << " particles");
  auto before = particles.size();

  auto particlesSet = particleVectorToSet(particles);

  if (particlesSet.size() < before) {
    ACTS_WARNING("Particle IDs not unique for " << before - particles.size()
                                                << " particles!");
  }

  return particlesSet;
}

std::tuple<ClusterContainer, MeasurementContainer,
           IndexMultimap<ActsFatras::Barcode>,
           std::unordered_map<int, std::size_t>>
RootAthenaDumpReader::readMeasurements(
    SimParticleContainer& particles, const Acts::GeometryContext& gctx) const {
  ClusterContainer clusters(nCL);
  clusters.reserve(nCL);

  MeasurementContainer measurements;
  measurements.reserve(nCL);

  std::size_t nTotalTotZero = 0;

  const auto prevParticlesSize = particles.size();
  IndexMultimap<ActsFatras::Barcode> measPartMap;

  // We cannot use im for the index since we might skip measurements
  std::unordered_map<int, std::size_t> imIdxMap;

  for (int im = 0; im < nCL; im++) {
    if (!(CLhardware->at(im) == "PIXEL" || CLhardware->at(im) == "STRIP")) {
      ACTS_ERROR("hardware is neither 'PIXEL' or 'STRIP', skip particle");
      continue;
    }
    ACTS_VERBOSE("Cluster " << im << ": " << CLhardware->at(im));

    auto type = (CLhardware->at(im) == "PIXEL") ? ePixel : eStrip;

    // Make cluster
    // TODO refactor Cluster class so it is not so tedious
    Cluster cluster;

    const auto& etas = CLetas->at(im);
    const auto& phis = CLetas->at(im);
    const auto& tots = CLtots->at(im);

    const auto totalTot = std::accumulate(tots.begin(), tots.end(), 0);

    const auto [minEta, maxEta] = std::minmax_element(etas.begin(), etas.end());
    const auto [minPhi, maxPhi] = std::minmax_element(phis.begin(), phis.end());

    cluster.sizeLoc0 = *maxEta - *minEta;
    cluster.sizeLoc1 = *maxPhi - *minPhi;

    if (totalTot == 0.0) {
      ACTS_VERBOSE("total time over threshold is 0, set all activations to 0");
      nTotalTotZero++;
    }

    for (const auto& [eta, phi, tot] : Acts::zip(etas, phis, tots)) {
      // Make best out of what we have:
      // Weight the overall collected charge corresponding to the
      // time-over-threshold of each cell Use this as activation (does this make
      // sense?)
      auto activation =
          (totalTot != 0.0) ? CLcharge_count[im] * tot / totalTot : 0.0;

      // This bases every cluster at zero, but shouldn't matter right now
      ActsFatras::Segmentizer::Bin2D bin;
      bin[0] = eta - *minEta;
      bin[1] = phi - *minPhi;

      // Of course we have no Segment2D because this is not Fatras
      cluster.channels.emplace_back(bin, ActsFatras::Segmentizer::Segment2D{},
                                    activation);
    }

    cluster.globalPosition = {CLx[im], CLy[im], CLz[im]};
    cluster.localDirection = {CLloc_direction1[im], CLloc_direction2[im],
                              CLloc_direction3[im]};
    cluster.lengthDirection = {CLJan_loc_direction1[im],
                               CLJan_loc_direction2[im],
                               CLJan_loc_direction3[im]};
    cluster.localEta = CLloc_eta[im];
    cluster.localPhi = CLloc_phi[im];
    cluster.globalEta = CLglob_eta[im];
    cluster.globalPhi = CLglob_phi[im];
    cluster.etaAngle = CLeta_angle[im];
    cluster.phiAngle = CLphi_angle[im];

    ACTS_VERBOSE("CL shape: " << cluster.channels.size()
                              << "cells, dimensions: " << cluster.sizeLoc0
                              << ", " << cluster.sizeLoc1);

    clusters[im] = cluster;

    // Measurement creation
    ACTS_VERBOSE("CL loc dims:" << CLloc_direction1[im] << ", "
                                << CLloc_direction2[im] << ", "
                                << CLloc_direction3[im]);
    const auto& locCov = CLlocal_cov->at(im);

    Acts::GeometryIdentifier geoId;
    std::vector<double> localParams;
    if (m_cfg.geometryIdMap && m_cfg.trackingGeometry) {
      const auto& geoIdMap = m_cfg.geometryIdMap->left;
      if (geoIdMap.find(CLmoduleID[im]) == geoIdMap.end()) {
        ACTS_WARNING("Missing geo id for " << CLmoduleID[im] << ", skip hit");
        continue;
      }

      geoId = m_cfg.geometryIdMap->left.at(CLmoduleID[im]);

      auto surface = m_cfg.trackingGeometry->findSurface(geoId);
      if (surface == nullptr) {
        ACTS_WARNING("Did not find " << geoId
                                     << " in tracking geometry, skip hit");
        continue;
      }

      bool inside =
          surface->isOnSurface(gctx, cluster.globalPosition, {},
                               Acts::BoundaryTolerance::AbsoluteEuclidean(
                                   m_cfg.absBoundaryTolerance),
                               std::numeric_limits<double>::max());

      if (!inside) {
        const Acts::Vector3 v =
            surface->transform(gctx).inverse() * cluster.globalPosition;
        ACTS_WARNING("Projected position is not in surface bounds for "
                     << surface->geometryId() << ", skip hit");
        ACTS_WARNING("Position in local coordinates: " << v.transpose());
        ACTS_WARNING("Surface details:\n" << surface->toStream(gctx));
        continue;
      }

      const double tol = (type == ePixel) ? Acts::s_onSurfaceTolerance : 1.3_mm;
      auto loc = surface->globalToLocal(gctx, cluster.globalPosition, {}, tol);

      if (!loc.ok()) {
        const Acts::Vector3 v =
            surface->transform(gctx).inverse() * cluster.globalPosition;
        ACTS_WARNING("Global-to-local fit failed on "
                     << geoId << " (z dist: " << v[2]
                     << ", projected on surface: " << std::boolalpha << inside
                     << ") , skip hit");
        continue;
      }

      // TODO is this in strip coordinates or in polar coordinates for annulus
      // bounds?
      localParams = std::vector<double>(loc->begin(), loc->end());
    } else {
      geoId = Acts::GeometryIdentifier(CLmoduleID[im]);
      localParams = {CLloc_direction1[im], CLloc_direction2[im]};
    }

    DigitizedParameters digiPars;
    if (type == ePixel) {
      digiPars.indices = {Acts::eBoundLoc0, Acts::eBoundLoc1};
      assert(locCov.size() == 4);
      digiPars.variances = {locCov[0], locCov[3]};
      digiPars.values = localParams;
    } else {
      digiPars.indices = {Acts::eBoundLoc0};
      assert(!locCov.empty());
      digiPars.variances = {locCov[0]};
      digiPars.values = {localParams[0]};
    }

    std::size_t measIndex = measurements.size();
    imIdxMap.emplace(im, measIndex);
    createMeasurement(measurements, geoId, digiPars);

    if (!m_cfg.noTruth) {
      // Create measurement particles map and particles container
      for (const auto& [subevt, barcode] :
           Acts::zip(CLparticleLink_eventIndex->at(im),
                     CLparticleLink_barcode->at(im))) {
        auto dummyBarcode = concatInts(barcode, subevt);
        // If we don't find the particle, create one with default values
        if (particles.find(dummyBarcode) == particles.end()) {
          ACTS_VERBOSE("Particle with subevt "
                       << subevt << ", barcode " << barcode
                       << "not found, create dummy one");
          particles.emplace(dummyBarcode, Acts::PdgParticle::eInvalid);
        }
        measPartMap.insert(
            std::pair<Index, ActsFatras::Barcode>{measIndex, dummyBarcode});
      }
    }
  }

  if (measurements.size() < static_cast<std::size_t>(nCL)) {
    ACTS_WARNING("Could not convert " << nCL - measurements.size() << " / "
                                      << nCL << " measurements");
  }

  if (particles.size() - prevParticlesSize > 0) {
    ACTS_DEBUG("Created " << particles.size() - prevParticlesSize
                          << " dummy particles");
  }

  if (nTotalTotZero > 0) {
    ACTS_DEBUG(nTotalTotZero << " / " << nCL
                             << " clusters have zero time-over-threshold");
  }

  return {clusters, measurements, measPartMap, imIdxMap};
}

std::tuple<SimSpacePointContainer, SimSpacePointContainer,
           SimSpacePointContainer>
RootAthenaDumpReader::readSpacepoints(
    const std::optional<std::unordered_map<int, std::size_t>>& imIdxMap) const {
  SimSpacePointContainer pixelSpacePoints;
  pixelSpacePoints.reserve(nSP);

  SimSpacePointContainer stripSpacePoints;
  stripSpacePoints.reserve(nSP);

  SimSpacePointContainer spacePoints;
  spacePoints.reserve(nSP);

  // Loop on space points
  std::size_t skippedSpacePoints = 0;
  for (int isp = 0; isp < nSP; isp++) {
    auto isPhiOverlap = (SPisOverlap[isp] == 2) || (SPisOverlap[isp] == 3);
    auto isEtaOverlap = (SPisOverlap[isp] == 1) || (SPisOverlap[isp] == 3);
    if (m_cfg.skipOverlapSPsPhi && isPhiOverlap) {
      ++skippedSpacePoints;
      continue;
    }
    if (m_cfg.skipOverlapSPsEta && isEtaOverlap) {
      ++skippedSpacePoints;
      continue;
    }

    const Acts::Vector3 globalPos{SPx[isp], SPy[isp], SPz[isp]};
    const double spCovr = SPcovr[isp];
    const double spCovz = SPcovz[isp];

    // PIX=1  STRIP = 2
    auto type = SPCL2_index[isp] == -1 ? ePixel : eStrip;

    ACTS_VERBOSE("SP:: " << type << " [" << globalPos.transpose() << "] "
                         << spCovr << " " << spCovz);

    boost::container::static_vector<Acts::SourceLink, 2> sLinks;

    const auto cl1Index = SPCL1_index[isp];
    assert(cl1Index >= 0 && cl1Index < nCL);

    auto getGeoId =
        [&](auto athenaId) -> std::optional<Acts::GeometryIdentifier> {
      if (m_cfg.geometryIdMap == nullptr) {
        return Acts::GeometryIdentifier{athenaId};
      }
      if (m_cfg.geometryIdMap->left.find(athenaId) ==
          m_cfg.geometryIdMap->left.end()) {
        return std::nullopt;
      }
      return m_cfg.geometryIdMap->left.at(athenaId);
    };

    auto cl1GeoId = getGeoId(CLmoduleID[cl1Index]);
    if (!cl1GeoId) {
      ACTS_WARNING("Could not find geoId for spacepoint cluster 1");
      continue;
    }

    if (imIdxMap && !imIdxMap->contains(cl1Index)) {
      ACTS_WARNING("Measurement 1 for spacepoint " << isp << " not created");
      continue;
    }

    IndexSourceLink first(*cl1GeoId,
                          imIdxMap ? imIdxMap->at(cl1Index) : cl1Index);
    sLinks.emplace_back(first);

    // First create pixel spacepoint here, later maybe overwrite with strip
    // spacepoint
    SimSpacePoint sp(globalPos, std::nullopt, spCovr, spCovz, std::nullopt,
                     sLinks);

    if (type == ePixel) {
      pixelSpacePoints.push_back(sp);
    } else {
      const auto cl2Index = SPCL2_index[isp];
      assert(cl2Index >= 0 && cl2Index < nCL);

      auto cl2GeoId = getGeoId(CLmoduleID[cl1Index]);
      if (!cl2GeoId) {
        ACTS_WARNING("Could not find geoId for spacepoint cluster 2");
        continue;
      }

      if (imIdxMap && !imIdxMap->contains(cl2Index)) {
        ACTS_WARNING("Measurement 2 for spacepoint " << isp << " not created");
        continue;
      }

      IndexSourceLink second(*cl2GeoId,
                             imIdxMap ? imIdxMap->at(cl2Index) : cl2Index);
      sLinks.emplace_back(second);

      using Vector3f = Eigen::Matrix<float, 3, 1>;
      Vector3f topStripDirection = Vector3f::Zero();
      Vector3f bottomStripDirection = Vector3f::Zero();
      Vector3f stripCenterDistance = Vector3f::Zero();
      Vector3f topStripCenterPosition = Vector3f::Zero();

      if (m_haveStripFeatures) {
        topStripDirection = {SPtopStripDirection->at(isp).at(0),
                             SPtopStripDirection->at(isp).at(1),
                             SPtopStripDirection->at(isp).at(2)};
        bottomStripDirection = {SPbottomStripDirection->at(isp).at(0),
                                SPbottomStripDirection->at(isp).at(1),
                                SPbottomStripDirection->at(isp).at(2)};
        stripCenterDistance = {SPstripCenterDistance->at(isp).at(0),
                               SPstripCenterDistance->at(isp).at(1),
                               SPstripCenterDistance->at(isp).at(2)};
        topStripCenterPosition = {SPtopStripCenterPosition->at(isp).at(0),
                                  SPtopStripCenterPosition->at(isp).at(1),
                                  SPtopStripCenterPosition->at(isp).at(2)};
      }
      sp = SimSpacePoint(globalPos, std::nullopt, spCovr, spCovz, std::nullopt,
                         sLinks, SPhl_topstrip[isp], SPhl_botstrip[isp],
                         topStripDirection.cast<double>(),
                         bottomStripDirection.cast<double>(),
                         stripCenterDistance.cast<double>(),
                         topStripCenterPosition.cast<double>());

      stripSpacePoints.push_back(sp);
    }

    spacePoints.push_back(sp);
  }

  if (m_cfg.skipOverlapSPsEta || m_cfg.skipOverlapSPsPhi) {
    ACTS_DEBUG("Skipped " << skippedSpacePoints
                          << " because of eta/phi overlaps");
  }
  if (spacePoints.size() <
      (static_cast<std::size_t>(nSP) - skippedSpacePoints)) {
    ACTS_WARNING("Could not convert " << nSP - spacePoints.size() << " of "
                                      << nSP << " spacepoints");
  }

  ACTS_DEBUG("Created " << spacePoints.size() << " overall space points");
  ACTS_DEBUG("Created " << pixelSpacePoints.size() << " "
                        << " pixel space points");

  return {spacePoints, pixelSpacePoints, stripSpacePoints};
}

std::pair<SimParticleContainer, IndexMultimap<ActsFatras::Barcode>>
RootAthenaDumpReader::reprocessParticles(
    const SimParticleContainer& particles,
    const IndexMultimap<ActsFatras::Barcode>& measPartMap) const {
  std::vector<ActsExamples::SimParticle> newParticles;
  newParticles.reserve(particles.size());
  IndexMultimap<ActsFatras::Barcode> newMeasPartMap;

  const auto partMeasMap = invertIndexMultimap(measPartMap);

  std::uint16_t primaryCount = 0;
  std::uint16_t secondaryCount = 0;

  for (const auto& particle : particles) {
    const auto [begin, end] = partMeasMap.equal_range(particle.particleId());

    if (begin == end) {
      ACTS_VERBOSE("Particle " << particle.particleId().value()
                               << " has no measurements");
      continue;
    }

    auto [athBarcode, athSubevent] = splitInt(particle.particleId().value());
    auto primary = (athBarcode < s_maxBarcodeForPrimary);

    ActsFatras::Barcode fatrasBarcode;

    // vertex primary shouldn't be zero for a valid particle
    fatrasBarcode.setVertexPrimary(1);
    if (primary) {
      fatrasBarcode.setVertexSecondary(0);
      fatrasBarcode.setParticle(primaryCount);
      assert(primaryCount < std::numeric_limits<std::uint16_t>::max());
      primaryCount++;
    } else {
      fatrasBarcode.setVertexSecondary(1);
      fatrasBarcode.setParticle(secondaryCount);
      assert(primaryCount < std::numeric_limits<std::uint16_t>::max());
      secondaryCount++;
    }

    auto newParticle = particle.withParticleId(fatrasBarcode);
    newParticle.final().setNumberOfHits(std::distance(begin, end));
    newParticles.push_back(newParticle);

    for (auto it = begin; it != end; ++it) {
      newMeasPartMap.insert(
          std::pair<Index, ActsFatras::Barcode>{it->second, fatrasBarcode});
    }
  }

  ACTS_DEBUG("After reprocessing particles " << newParticles.size() << " of "
                                             << particles.size() << " remain");
  return {particleVectorToSet(newParticles), newMeasPartMap};
}

ProcessCode RootAthenaDumpReader::read(const AlgorithmContext& ctx) {
  ACTS_DEBUG("Reading event " << ctx.eventNumber);
  auto entry = ctx.eventNumber;
  if (entry >= m_events) {
    ACTS_ERROR("event out of bounds");
    return ProcessCode::ABORT;
  }

  std::lock_guard<std::mutex> lock(m_read_mutex);

  m_inputchain->GetEntry(entry);

  std::optional<std::unordered_map<int, std::size_t>> optImIdxMap;

  if (!m_cfg.onlySpacepoints) {
    SimParticleContainer candidateParticles;

    if (!m_cfg.noTruth) {
      candidateParticles = readParticles();
    }

    auto [clusters, measurements, candidateMeasPartMap, imIdxMap] =
        readMeasurements(candidateParticles, ctx.geoContext);
    optImIdxMap.emplace(std::move(imIdxMap));

    m_outputClusters(ctx, std::move(clusters));
    m_outputMeasurements(ctx, std::move(measurements));

    if (!m_cfg.noTruth) {
      auto [particles, measPartMap] =
          reprocessParticles(candidateParticles, candidateMeasPartMap);

      m_outputParticles(ctx, std::move(particles));
      m_outputParticleMeasMap(ctx, invertIndexMultimap(measPartMap));
      m_outputMeasParticleMap(ctx, std::move(measPartMap));
    }
  }

  auto [spacePoints, pixelSpacePoints, stripSpacePoints] =
      readSpacepoints(optImIdxMap);

  m_outputPixelSpacePoints(ctx, std::move(pixelSpacePoints));
  m_outputStripSpacePoints(ctx, std::move(stripSpacePoints));
  m_outputSpacePoints(ctx, std::move(spacePoints));

  return ProcessCode::SUCCESS;
}
}  // namespace ActsExamples
