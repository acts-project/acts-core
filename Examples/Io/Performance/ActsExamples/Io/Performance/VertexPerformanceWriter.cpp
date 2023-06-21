// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Performance/VertexPerformanceWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/SingleBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/MultiIndex.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ios>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <stdexcept>
#include <utility>

#include <TFile.h>
#include <TTree.h>

namespace ActsExamples {
struct AlgorithmContext;
}  // namespace ActsExamples

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

ActsExamples::VertexPerformanceWriter::VertexPerformanceWriter(
    const ActsExamples::VertexPerformanceWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT(config.inputVertices, "VertexPerformanceWriter", level),
      m_cfg(config) {
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }
  if (m_cfg.inputAllTruthParticles.empty()) {
    throw std::invalid_argument("Collection with all truth particles missing");
  }
  if (m_cfg.inputSelectedTruthParticles.empty()) {
    throw std::invalid_argument(
        "Collection with selected truth particles missing");
  }

  m_inputAllTruthParticles.initialize(m_cfg.inputAllTruthParticles);
  m_inputSelectedTruthParticles.initialize(m_cfg.inputSelectedTruthParticles);

  if (!m_cfg.inputAssociatedTruthParticles.empty()) {
    m_inputAssociatedTruthParticles.initialize(
        m_cfg.inputAssociatedTruthParticles);
    if (!m_cfg.inputTrackParameters.empty()) {
      m_inputTrackParameters.initialize(m_cfg.inputTrackParameters);
    } else {
      m_inputTrajectories.initialize(m_cfg.inputTrajectories);
    }
  } else {
    m_inputMeasurementParticlesMap.initialize(
        m_cfg.inputMeasurementParticlesMap);
    m_inputTrajectories.initialize(m_cfg.inputTrajectories);
  }

  // Setup ROOT I/O
  auto path = m_cfg.filePath;
  m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + path);
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  } else {
    // I/O parameters
    m_outputTree->Branch("truthX", &m_truthX);
    m_outputTree->Branch("truthY", &m_truthY);
    m_outputTree->Branch("truthZ", &m_truthZ);
    m_outputTree->Branch("truthT", &m_truthT);
    m_outputTree->Branch("truthPhi", &m_truthPhi);
    m_outputTree->Branch("truthTheta", &m_truthTheta);
    m_outputTree->Branch("truthQOverP", &m_truthQOverP);

    m_outputTree->Branch("recoX", &m_recoX);
    m_outputTree->Branch("recoY", &m_recoY);
    m_outputTree->Branch("recoZ", &m_recoZ);
    m_outputTree->Branch("recoT", &m_recoT);
    m_outputTree->Branch("recoPhi", &m_recoPhi);
    m_outputTree->Branch("recoPhiFitted", &m_recoPhiFitted);
    m_outputTree->Branch("recoTheta", &m_recoTheta);
    m_outputTree->Branch("recoThetaFitted", &m_recoThetaFitted);
    m_outputTree->Branch("recoQOverP", &m_recoQOverP);
    m_outputTree->Branch("recoQOverPFitted", &m_recoQOverPFitted);

    m_outputTree->Branch("resX", &m_resX);
    m_outputTree->Branch("resY", &m_resY);
    m_outputTree->Branch("resZ", &m_resZ);
    m_outputTree->Branch("resT", &m_resT);
    m_outputTree->Branch("resPhi", &m_resPhi);
    m_outputTree->Branch("resPhiFitted", &m_resPhiFitted);
    m_outputTree->Branch("resTheta", &m_resTheta);
    m_outputTree->Branch("resThetaFitted", &m_resThetaFitted);
    m_outputTree->Branch("resQOverP", &m_resQOverP);
    m_outputTree->Branch("resQOverPFitted", &m_resQOverPFitted);
    m_outputTree->Branch("momOverlap", &m_momOverlap);
    m_outputTree->Branch("momOverlapFitted", &m_momOverlapFitted);

    m_outputTree->Branch("pullX", &m_pullX);
    m_outputTree->Branch("pullY", &m_pullY);
    m_outputTree->Branch("pullZ", &m_pullZ);
    m_outputTree->Branch("pullT", &m_pullT);
    m_outputTree->Branch("pullPhi", &m_pullPhi);
    m_outputTree->Branch("pullPhiFitted", &m_pullPhiFitted);
    m_outputTree->Branch("pullTheta", &m_pullTheta);
    m_outputTree->Branch("pullThetaFitted", &m_pullThetaFitted);
    m_outputTree->Branch("pullQOverP", &m_pullQOverP);
    m_outputTree->Branch("pullQOverPFitted", &m_pullQOverPFitted);

    m_outputTree->Branch("covXX", &m_covXX);
    m_outputTree->Branch("covYY", &m_covYY);
    m_outputTree->Branch("covZZ", &m_covZZ);
    m_outputTree->Branch("covTT", &m_covTT);
    m_outputTree->Branch("covXY", &m_covXY);
    m_outputTree->Branch("covXZ", &m_covXZ);
    m_outputTree->Branch("covXT", &m_covXT);
    m_outputTree->Branch("covYZ", &m_covYZ);
    m_outputTree->Branch("covYT", &m_covYT);
    m_outputTree->Branch("covZT", &m_covZT);

    m_outputTree->Branch("nTracksTruthVtx", &m_nTracksOnTruthVertex);
    m_outputTree->Branch("nTracksRecoVtx", &m_nTracksOnRecoVertex);

    m_outputTree->Branch("trkVtxMatch", &m_trackVtxMatchFraction);

    m_outputTree->Branch("nRecoVtx", &m_nRecoVtx);
    m_outputTree->Branch("nTrueVtx", &m_nTrueVtx);
    m_outputTree->Branch("nVtxDetectorAcceptance", &m_nVtxDetAcceptance);
    m_outputTree->Branch("nVtxReconstructable", &m_nVtxReconstructable);
  }
}

ActsExamples::VertexPerformanceWriter::~VertexPerformanceWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::VertexPerformanceWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  return ProcessCode::SUCCESS;
}

int ActsExamples::VertexPerformanceWriter::getNumberOfReconstructableVertices(
    const SimParticleContainer& collection) const {
  // map for finding frequency
  std::map<int, int> fmap;

  std::vector<int> reconstructableTruthVertices;

  // traverse the array for frequency
  for (const auto& p : collection) {
    int secVtxId = p.particleId().vertexSecondary();
    if (secVtxId != 0) {
      // truthparticle from secondary vtx
      continue;
    }
    int priVtxId = p.particleId().vertexPrimary();
    fmap[priVtxId]++;
  }

  // iterate over the map
  for (auto it : fmap) {
    // Require at least 2 tracks
    if (it.second > 1) {
      reconstructableTruthVertices.push_back(it.first);
    }
  }

  return reconstructableTruthVertices.size();
}

int ActsExamples::VertexPerformanceWriter::getNumberOfTruePriVertices(
    const SimParticleContainer& collection) const {
  // Vector to store indices of all primary vertices
  std::set<int> allPriVtxIds;
  for (const auto& p : collection) {
    int priVtxId = p.particleId().vertexPrimary();
    int secVtxId = p.particleId().vertexSecondary();
    if (secVtxId != 0) {
      // truthparticle from secondary vtx
      continue;
    }
    // Insert to set, removing duplicates
    allPriVtxIds.insert(priVtxId);
  }
  // Size of set corresponds to total number of primary vertices
  return allPriVtxIds.size();
}

ActsExamples::ProcessCode ActsExamples::VertexPerformanceWriter::writeT(
    const AlgorithmContext& ctx,
    const std::vector<Acts::Vertex<Acts::BoundTrackParameters>>& vertices) {
  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  m_nRecoVtx = vertices.size();

  ACTS_DEBUG("Number of reco vertices in event: " << m_nRecoVtx);

  // Read truth particle input collection
  const auto& allTruthParticles = m_inputAllTruthParticles(ctx);
  // Get number of generated true primary vertices
  m_nTrueVtx = getNumberOfTruePriVertices(allTruthParticles);

  ACTS_VERBOSE("Total number of generated truth particles in event : "
               << allTruthParticles.size());
  ACTS_VERBOSE(
      "Total number of generated truth primary vertices : " << m_nTrueVtx);

  // Read selected truth particle input collection
  const auto& selectedTruthParticles = m_inputSelectedTruthParticles(ctx);
  // Get number of detector-accepted true primary vertices
  m_nVtxDetAcceptance = getNumberOfTruePriVertices(selectedTruthParticles);

  ACTS_VERBOSE("Total number of selected truth particles in event : "
               << selectedTruthParticles.size());
  ACTS_VERBOSE("Total number of detector-accepted truth primary vertices : "
               << m_nVtxDetAcceptance);

  std::vector<Acts::BoundTrackParameters> trackParameters;
  std::vector<SimParticle> associatedTruthParticles;

  // The i-th entry in associatedTruthParticles corresponds to the i-th entry in
  // trackParameters. If we know the truth particles associated to the track
  // parameters a priori:
  if (!m_cfg.inputAssociatedTruthParticles.empty()) {
    if (!m_cfg.inputTrackParameters.empty()) {
      trackParameters = m_inputTrackParameters(ctx);
    } else {
      const auto& inputTrajectories = m_inputTrajectories(ctx);

      for (const auto& trajectories : inputTrajectories) {
        for (auto tip : trajectories.tips()) {
          if (!trajectories.hasTrackParameters(tip)) {
            continue;
          }
          trackParameters.push_back(trajectories.trackParameters(tip));
        }
      }
    }

    // Read track-associated truth particle input collection
    associatedTruthParticles =
        std::vector<SimParticle>(m_inputAssociatedTruthParticles(ctx).begin(),
                                 m_inputAssociatedTruthParticles(ctx).end());

    auto mismatchMsg = [&](auto level, const auto& extra) {
      ACTS_LOG(level,
               "Number of fitted tracks and associated truth particles do not "
               "match. ("
                   << trackParameters.size()
                   << " != " << associatedTruthParticles.size()
                   << ") Not able to match fitted tracks at reconstructed "
                      "vertex to truth vertex."
                   << extra);
    };

    if (associatedTruthParticles.size() < trackParameters.size()) {
      mismatchMsg(Acts::Logging::ERROR, " Switch to hit based truth matching.");
    } else if (associatedTruthParticles.size() > trackParameters.size()) {
      mismatchMsg(Acts::Logging::INFO,
                  " This is likely due to track efficiency < 1");
    }
  }
  // If we don't know which truth particle corresponds to which track a priori,
  // we check how many hits particles and tracks share. We match the particle
  // to the track if a fraction of more than truthMatchProbMin of hits that
  // contribute to the track come from the particle.
  // Note that not all tracksatVertex have matching parameters in
  // trackParameters in this case. Equivalently, one could say that not all
  // tracksAtVertex will be assigned to a truth particle.
  else {
    // get active tips
    const auto& inputTrajectories = m_inputTrajectories(ctx);

    std::vector<ParticleHitCount> particleHitCounts;

    const auto& hitParticlesMap = m_inputMeasurementParticlesMap(ctx);

    for (const auto& trajectories : inputTrajectories) {
      for (auto tip : trajectories.tips()) {
        if (!trajectories.hasTrackParameters(tip)) {
          continue;
        }

        identifyContributingParticles(hitParticlesMap, trajectories, tip,
                                      particleHitCounts);
        ActsFatras::Barcode majorityParticleId =
            particleHitCounts.front().particleId;
        size_t nMajorityHits = particleHitCounts.front().hitCount;

        auto trajState = Acts::MultiTrajectoryHelpers::trajectoryState(
            trajectories.multiTrajectory(), tip);

        if (nMajorityHits / trajState.nMeasurements < m_cfg.truthMatchProbMin) {
          continue;
        }

        auto it = std::find_if(allTruthParticles.begin(),
                               allTruthParticles.end(), [&](const auto& tp) {
                                 return tp.particleId() == majorityParticleId;
                               });

        if (it == allTruthParticles.end()) {
          continue;
        }

        const auto& majorityParticle = *it;
        trackParameters.push_back(trajectories.trackParameters(tip));
        associatedTruthParticles.push_back(majorityParticle);
      }
    }
  }

  // Get number of track-associated true primary vertices
  m_nVtxReconstructable =
      getNumberOfReconstructableVertices(SimParticleContainer(
          associatedTruthParticles.begin(), associatedTruthParticles.end()));

  ACTS_INFO(
      "Total number of reconstructed tracks : " << trackParameters.size());
  ACTS_INFO("Total number of reco track-associated truth particles in event : "
            << associatedTruthParticles.size());
  ACTS_INFO("Total number of reco track-associated truth primary vertices : "
            << m_nVtxReconstructable);

  // Loop over all reco vertices and find associated truth particles
  std::vector<SimParticleContainer> truthParticlesAtVtxContainer;

  // We compare the reconstructed momenta to the true momenta at the vertex. For
  // this, we propagate the reconstructed tracks to the PCA of the true vertex
  // position. Setting up propagator:
  Acts::EigenStepper<> stepper(m_cfg.bField);
  using Propagator = Acts::Propagator<Acts::EigenStepper<>>;
  auto propagator = std::make_shared<Propagator>(stepper);

  // Loop over reconstructed vertices and see if they can be matched to a true
  // vertex.
  for (const auto& vtx : vertices) {
    // Reconstructed tracks that contribute to the reconstructed vertex
    const auto tracksAtVtx = vtx.tracks();

    // Containers for storing truth particles and truth vertices that contribute
    // to the reconstructed vertex
    SimParticleContainer particleAtVtx;
    std::vector<int> contributingTruthVertices;

    for (const auto& trk : tracksAtVtx) {
      // Track parameters before the vertex fit
      Acts::BoundTrackParameters origTrack = *(trk.originalParams);

      // Finding the matching parameters in the container of all track
      // parameters. This allows us to identify the corresponding particle,
      // since we expect trackParameters and associatedTruthParticles to align.
      bool foundMatchingParams = false;
      for (std::size_t i = 0; i < trackParameters.size(); ++i) {
        const auto& params = trackParameters[i].parameters();

        if (origTrack.parameters() == params) {
          // We expect that the i-th associated truth particle corresponds to
          // the i-th track parameters
          const auto& particle = associatedTruthParticles[i];
          particleAtVtx.insert(particle);
          int priVtxId = particle.particleId().vertexPrimary();
          contributingTruthVertices.push_back(priVtxId);
          foundMatchingParams = true;
          break;
        }
      }
      if (!foundMatchingParams) {
        ACTS_VERBOSE("Track has no matching truth particle.");
      }
    }  // end loop tracksAtVtx

    // Find true vertex that contributes most to the reconstructed vertex
    std::map<int, int> fmap;
    for (int priVtxId : contributingTruthVertices) {
      fmap[priVtxId]++;
    }
    int maxOccurrence = -1;
    int maxOccurrenceId = -1;
    for (auto it : fmap) {
      if (it.second > maxOccurrence) {
        maxOccurrenceId = it.first;
        maxOccurrence = it.second;
      }
    }

    // Count number of reconstructible tracks on truth vertex
    int nTracksOnTruthVertex = 0;
    for (const auto& particle : associatedTruthParticles) {
      int priVtxId = particle.particleId().vertexPrimary();
      if (priVtxId == maxOccurrenceId) {
        ++nTracksOnTruthVertex;
      }
    }

    // Match reconstructed and truth vertex if the tracks of the truth vertex
    // make up at least minTrackVtxMatchFraction of the tracks at the
    // reconstructed vertex.
    double trackVtxMatchFraction =
        (double)fmap[maxOccurrenceId] / tracksAtVtx.size();
    if (trackVtxMatchFraction > m_cfg.minTrackVtxMatchFraction) {
      int count = 0;
      for (std::size_t j = 0; j < associatedTruthParticles.size(); ++j) {
        const auto& particle = associatedTruthParticles[j];
        int priVtxId = particle.particleId().vertexPrimary();
        int secVtxId = particle.particleId().vertexSecondary();

        if (secVtxId != 0) {
          // truthparticle from secondary vtx
          continue;
        }
        if (priVtxId == maxOccurrenceId) {
          // Vertex found, fill variables
          const auto& truePos = particle.fourPosition();

          auto pull = [&](int i, const auto& covariance, const auto& recoVec,
                          const auto& trueVec) {
            double var = covariance(i, i);
            if (var < 0) {
              ACTS_WARNING("var(" << i << ") = " << var << " < 0");
              return std::numeric_limits<double>::quiet_NaN();
            }
            double std = std::sqrt(var);
            if (std == 0) {
              ACTS_WARNING("std(" << i << ") = 0");
              return std::numeric_limits<double>::quiet_NaN();
            }
            return (recoVec[i] - trueVec[i]) / std;
          };

          // Save reconstructed/true vertex position only in the first iteration
          // to avoid duplicates
          if (count == 0) {
            m_truthX.push_back(truePos[Acts::FreeIndices::eFreePos0]);
            m_truthY.push_back(truePos[Acts::FreeIndices::eFreePos1]);
            m_truthZ.push_back(truePos[Acts::FreeIndices::eFreePos2]);
            m_truthT.push_back(truePos[Acts::FreeIndices::eFreeTime]);

            m_recoX.push_back(vtx.fullPosition()[Acts::FreeIndices::eFreePos0]);
            m_recoY.push_back(vtx.fullPosition()[Acts::FreeIndices::eFreePos1]);
            m_recoZ.push_back(vtx.fullPosition()[Acts::FreeIndices::eFreePos2]);
            m_recoT.push_back(vtx.fullPosition()[Acts::FreeIndices::eFreeTime]);

            m_resX.push_back(vtx.fullPosition()[Acts::FreeIndices::eFreePos0] -
                             truePos[Acts::FreeIndices::eFreePos0]);
            m_resY.push_back(vtx.fullPosition()[Acts::FreeIndices::eFreePos1] -
                             truePos[Acts::FreeIndices::eFreePos1]);
            m_resZ.push_back(vtx.fullPosition()[Acts::FreeIndices::eFreePos2] -
                             truePos[Acts::FreeIndices::eFreePos2]);
            m_resT.push_back(vtx.fullPosition()[Acts::FreeIndices::eFreeTime] -
                             truePos[Acts::FreeIndices::eFreeTime]);

            m_pullX.push_back(pull(Acts::FreeIndices::eFreePos0,
                                   vtx.fullCovariance(), vtx.fullPosition(),
                                   truePos));
            m_pullY.push_back(pull(Acts::FreeIndices::eFreePos1,
                                   vtx.fullCovariance(), vtx.fullPosition(),
                                   truePos));
            m_pullZ.push_back(pull(Acts::FreeIndices::eFreePos2,
                                   vtx.fullCovariance(), vtx.fullPosition(),
                                   truePos));
            m_pullT.push_back(pull(Acts::FreeIndices::eFreeTime,
                                   vtx.fullCovariance(), vtx.fullPosition(),
                                   truePos));

            m_covXX.push_back(vtx.fullCovariance()(
                Acts::FreeIndices::eFreePos0, Acts::FreeIndices::eFreePos0));
            m_covYY.push_back(vtx.fullCovariance()(
                Acts::FreeIndices::eFreePos1, Acts::FreeIndices::eFreePos1));
            m_covZZ.push_back(vtx.fullCovariance()(
                Acts::FreeIndices::eFreePos2, Acts::FreeIndices::eFreePos2));
            m_covTT.push_back(vtx.fullCovariance()(
                Acts::FreeIndices::eFreeTime, Acts::FreeIndices::eFreeTime));
            m_covXY.push_back(vtx.fullCovariance()(
                Acts::FreeIndices::eFreePos0, Acts::FreeIndices::eFreePos1));
            m_covXZ.push_back(vtx.fullCovariance()(
                Acts::FreeIndices::eFreePos0, Acts::FreeIndices::eFreePos2));
            m_covXT.push_back(vtx.fullCovariance()(
                Acts::FreeIndices::eFreePos0, Acts::FreeIndices::eFreeTime));
            m_covYZ.push_back(vtx.fullCovariance()(
                Acts::FreeIndices::eFreePos1, Acts::FreeIndices::eFreePos2));
            m_covYT.push_back(vtx.fullCovariance()(
                Acts::FreeIndices::eFreePos1, Acts::FreeIndices::eFreeTime));
            m_covZT.push_back(vtx.fullCovariance()(
                Acts::FreeIndices::eFreePos2, Acts::FreeIndices::eFreeTime));

            m_nTracksOnTruthVertex.push_back(nTracksOnTruthVertex);
            m_nTracksOnRecoVertex.push_back(tracksAtVtx.size());

            m_trackVtxMatchFraction.push_back(trackVtxMatchFraction);
          }

          // Saving the reconstructed/truth momenta. The reconstructed momenta
          // are taken at the PCA to the truth vertex position -> we need to
          // perform a propagation.

          // Perigee at the true vertex position
          const std::shared_ptr<Acts::PerigeeSurface> perigeeSurface =
              Acts::Surface::makeShared<Acts::PerigeeSurface>(truePos.head(3));
          // Setting the geometry/magnetic field context context for the event
          Acts::PropagatorOptions pOptions(ctx.geoContext, ctx.magFieldContext);
          // Lambda for propagating the tracks to the PCA
          auto propagateToVtx = [&](const auto& params) {
            auto intersection = perigeeSurface->intersect(
                ctx.geoContext, params.position(ctx.geoContext),
                params.unitDirection(), false);
            pOptions.direction = Acts::Direction::fromScalarZeroAsPositive(
                intersection.intersection.pathLength);

            auto result = propagator->propagate(
                params, *perigeeSurface,
                pOptions);  // TODO: not checking if result is ok - is that a
                            // problem?
            auto& paramsAtVtx = *result->endParameters;
            return paramsAtVtx;
          };
          // Get the reconstructed track parameters corresponding to the
          // particle
          const auto& params = trackParameters[j].parameters();
          // Check if they correspond to a track that contributed to the vertex.
          // We save the momenta if we find a match.
          for (const auto& trk : tracksAtVtx) {
            if (trk.originalParams->parameters() == params) {
              const auto& trueUnitDir = particle.unitDirection();
              Acts::ActsVector<3> trueMom;
              trueMom.head(2) =
                  Acts::makePhiThetaFromDirectionUnit(trueUnitDir);
              trueMom[2] = particle.qop();
              m_truthPhi.push_back(trueMom[0]);
              m_truthTheta.push_back(trueMom[1]);
              m_truthQOverP.push_back(trueMom[2]);

              // Track parameters before the vertex fit
              auto paramsAtVtx = propagateToVtx(*(trk.originalParams));
              Acts::ActsVector<3> recoMom =
                  paramsAtVtx.parameters().segment(Acts::eBoundPhi, 3);
              const Acts::ActsMatrix<3, 3>& momCov =
                  paramsAtVtx.covariance()->block<3, 3>(Acts::eBoundPhi,
                                                        Acts::eBoundPhi);
              m_recoPhi.push_back(recoMom[0]);
              m_recoTheta.push_back(recoMom[1]);
              m_recoQOverP.push_back(recoMom[2]);

              // TODO: subtract angles correctly
              m_resPhi.push_back(recoMom[0] - trueMom[0]);
              m_resTheta.push_back(recoMom[1] - trueMom[1]);
              m_resQOverP.push_back(recoMom[2] - trueMom[2]);

              m_pullPhi.push_back(pull(0, momCov, recoMom, trueMom));
              m_pullTheta.push_back(pull(1, momCov, recoMom, trueMom));
              m_pullQOverP.push_back(pull(2, momCov, recoMom, trueMom));

              const auto& recoUnitDir = paramsAtVtx.unitDirection();
              double overlap = trueUnitDir.dot(recoUnitDir);
              m_momOverlap.push_back(overlap);

              // Track parameters after the vertex fit
              auto paramsAtVtxFitted = propagateToVtx(trk.fittedParams);
              Acts::ActsVector<3> recoMomFitted =
                  paramsAtVtxFitted.parameters().segment(Acts::eBoundPhi, 3);
              const Acts::ActsMatrix<3, 3>& momCovFitted =
                  paramsAtVtxFitted.covariance()->block<3, 3>(Acts::eBoundPhi,
                                                              Acts::eBoundPhi);
              m_recoPhiFitted.push_back(recoMomFitted[0]);
              m_recoThetaFitted.push_back(recoMomFitted[1]);
              m_recoQOverPFitted.push_back(recoMomFitted[2]);

              // TODO: subtract angles correctly
              m_resPhiFitted.push_back(recoMomFitted[0] - trueMom[0]);
              m_resThetaFitted.push_back(recoMomFitted[1] - trueMom[1]);
              m_resQOverPFitted.push_back(recoMomFitted[2] - trueMom[2]);

              m_pullPhiFitted.push_back(
                  pull(0, momCovFitted, recoMomFitted, trueMom));
              m_pullThetaFitted.push_back(
                  pull(1, momCovFitted, recoMomFitted, trueMom));
              m_pullQOverPFitted.push_back(
                  pull(2, momCovFitted, recoMomFitted, trueMom));

              const auto& recoUnitDirFitted = paramsAtVtxFitted.unitDirection();
              double overlapFitted = trueUnitDir.dot(recoUnitDirFitted);
              m_momOverlapFitted.push_back(overlapFitted);
            }
          }
          count++;
        }
      }
    }
  }  // end loop vertices

  // fill the variables
  m_outputTree->Fill();

  m_truthX.clear();
  m_truthY.clear();
  m_truthZ.clear();
  m_truthT.clear();
  m_truthPhi.clear();
  m_truthTheta.clear();
  m_truthQOverP.clear();
  m_recoX.clear();
  m_recoY.clear();
  m_recoZ.clear();
  m_recoT.clear();
  m_recoPhi.clear();
  m_recoPhiFitted.clear();
  m_recoTheta.clear();
  m_recoThetaFitted.clear();
  m_recoQOverP.clear();
  m_recoQOverPFitted.clear();
  m_resX.clear();
  m_resY.clear();
  m_resZ.clear();
  m_resT.clear();
  m_resPhi.clear();
  m_resPhiFitted.clear();
  m_resTheta.clear();
  m_resThetaFitted.clear();
  m_resQOverP.clear();
  m_resQOverPFitted.clear();
  m_momOverlap.clear();
  m_momOverlapFitted.clear();
  m_pullX.clear();
  m_pullY.clear();
  m_pullZ.clear();
  m_pullT.clear();
  m_pullPhi.clear();
  m_pullPhiFitted.clear();
  m_pullTheta.clear();
  m_pullThetaFitted.clear();
  m_pullQOverP.clear();
  m_pullQOverPFitted.clear();
  m_covXX.clear();
  m_covYY.clear();
  m_covZZ.clear();
  m_covTT.clear();
  m_covXY.clear();
  m_covXZ.clear();
  m_covXT.clear();
  m_covYZ.clear();
  m_covYT.clear();
  m_covZT.clear();

  m_nTracksOnTruthVertex.clear();
  m_nTracksOnRecoVertex.clear();

  m_trackVtxMatchFraction.clear();

  return ProcessCode::SUCCESS;
}
