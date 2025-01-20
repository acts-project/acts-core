// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// SeedFinderGbts.ipp
// TODO: update to C++17 style

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Seeding/SeedFinder.hpp"
#include "Acts/Seeding/SeedFinderGbtsConfig.hpp"
#include "Acts/Seeding/SeedFinderUtils.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <list>
#include <numbers>
#include <numeric>
#include <type_traits>
#include <vector>

// core so in ACTS namespace

namespace Acts {

template <typename external_spacepoint_t>
SeedFinderGbts<external_spacepoint_t>::SeedFinderGbts(
    const SeedFinderGbtsConfig<external_spacepoint_t>& config,
    const GbtsGeometry<external_spacepoint_t>& gbtsGeo,
    std::unique_ptr<const Acts::Logger> logger)
    : m_config(config),
      m_storage(
          std::make_unique<GbtsDataStorage<external_spacepoint_t>>(gbtsGeo)),
      m_logger(std::move(logger)) {}

// define loadspace points function
template <typename external_spacepoint_t>
void SeedFinderGbts<external_spacepoint_t>::loadSpacePoints(
    const std::vector<GbtsSP<external_spacepoint_t>>& gbtsSPvect) {
  ACTS_VERBOSE("Loading space points");
  for (const auto& gbtssp : gbtsSPvect) {
    bool is_Pixel = gbtssp.isPixel();
    if (!is_Pixel) {
      continue;
    }
    m_storage->addSpacePoint(gbtssp, (m_config.m_useClusterWidth > 0));
  }

  m_config.m_phiSliceWidth = 2 * std::numbers::pi / m_config.m_nMaxPhiSlice;

  m_storage->sortByPhi();

  m_storage->generatePhiIndexing(1.5 * m_config.m_phiSliceWidth);
}

template <typename external_spacepoint_t>
void SeedFinderGbts<external_spacepoint_t>::runGbts_TrackFinder(
    std::vector<GbtsTrigTracklet<external_spacepoint_t>>& vTracks,
    const Acts::RoiDescriptor& roi,
    const Acts::GbtsGeometry<external_spacepoint_t>& gbtsGeo) {
  ACTS_VERBOSE("Running GBTS Track Finder");
  const float min_z0 = roi.zedMinus();
  const float max_z0 = roi.zedPlus();
  const float cut_zMinU = min_z0 + m_config.maxOuterRadius * roi.dzdrMinus();
  const float cut_zMaxU = max_z0 + m_config.maxOuterRadius * roi.dzdrPlus();
  float m_minR_squ = m_config.m_tripletPtMin * m_config.m_tripletPtMin /
                     std::pow(m_config.ptCoeff, 2);  // from athena
  float m_maxCurv = m_config.ptCoeff / m_config.m_tripletPtMin;
  const float maxKappa_high_eta = 0.8 / m_minR_squ;
  const float maxKappa_low_eta = 0.6 / m_minR_squ;

  // 1. loop over stages

  int currentStage = 0;

  const Acts::GbtsConnector& connector = *(gbtsGeo.connector());

  std::vector<Acts::GbtsEdge<external_spacepoint_t>> edgeStorage;

  edgeStorage.reserve(m_config.MaxEdges);

  int nEdges = 0;

  for (std::map<int, std::vector<GbtsConnector::LayerGroup>>::const_iterator
           it = connector.m_layerGroups.begin();
       it != connector.m_layerGroups.end(); ++it, currentStage++) {
    // loop over L1 layers for the current stage

    for (const auto& layerGroup : (*it).second) {
      unsigned int dst = layerGroup.m_dst;  // n1 : inner nodes

      const GbtsLayer<external_spacepoint_t>* pL1 =
          gbtsGeo.getGbtsLayerByKey(dst);

      if (pL1 == nullptr) {
        continue;
      }

      for (const auto& conn :
           layerGroup.m_sources) {  // loop over L2(L1) for the current stage

        unsigned int src = conn->m_src;  // n2 : the new connectors

        const GbtsLayer<external_spacepoint_t>* pL2 =
            gbtsGeo.getGbtsLayerByKey(src);

        if (pL2 == nullptr) {
          continue;
        }
        int nDstBins = pL1->m_bins.size();
        int nSrcBins = pL2->m_bins.size();

        for (int b1 = 0; b1 < nDstBins; b1++) {  // loop over bins in Layer 1

          const GbtsEtaBin<external_spacepoint_t>& B1 =
              m_storage->getEtaBin(pL1->m_bins.at(b1));

          if (B1.empty()) {
            continue;
          }

          float rb1 = pL1->getMinBinRadius(b1);

          // 3. loops over source eta-bins

          for (int b2 = 0; b2 < nSrcBins; b2++) {  // loop over bins in Layer 2

            if (m_config.m_useEtaBinning && (nSrcBins + nDstBins > 2)) {
              if (conn->m_binTable[b1 + b2 * nDstBins] != 1) {
                continue;  // using precomputed LUT
              }
            }

            const GbtsEtaBin<external_spacepoint_t>& B2 =
                m_storage->getEtaBin(pL2->m_bins.at(b2));

            if (B2.empty()) {
              continue;
            }

            float rb2 = pL2->getMaxBinRadius(b2);

            // calculated delta Phi for rb1 ---> rb2 extrapolation

            float deltaPhi =
                0.5f *
                m_config
                    .m_phiSliceWidth;  // the default sliding window along phi

            if (m_config.m_useEtaBinning) {
              deltaPhi = 0.001f + m_maxCurv * std::abs(rb2 - rb1);
            }

            unsigned int first_it = 0;
            for (const auto& n1 : B1.m_vn) {  // loop over nodes in Layer 1

              if (n1->m_in.size() >= MAX_SEG_PER_NODE) {
                continue;
              }

              float r1 = n1->m_spGbts.r();
              float x1 = n1->m_spGbts.SP->x();
              float y1 = n1->m_spGbts.SP->y();
              float z1 = n1->m_spGbts.SP->z();
              float phi1 = n1->m_spGbts.phi();

              float minPhi = phi1 - deltaPhi;
              float maxPhi = phi1 + deltaPhi;

              for (unsigned int n2PhiIdx = first_it;
                   n2PhiIdx < B2.m_vPhiNodes.size();
                   n2PhiIdx++) {  // sliding window over nodes in Layer 2

                float phi2 = B2.m_vPhiNodes.at(n2PhiIdx).first;

                if (phi2 < minPhi) {
                  first_it = n2PhiIdx;
                  continue;
                }
                if (phi2 > maxPhi) {
                  break;
                }

                GbtsNode<external_spacepoint_t>* n2 =
                    B2.m_vn.at(B2.m_vPhiNodes.at(n2PhiIdx).second).get();

                if (n2->m_out.size() >= MAX_SEG_PER_NODE) {
                  continue;
                }
                if (n2->isFull()) {
                  continue;
                }

                float r2 = n2->m_spGbts.r();

                float dr = r2 - r1;

                if (dr < m_config.m_minDeltaRadius) {
                  continue;
                }

                float z2 = n2->m_spGbts.SP->z();

                float dz = z2 - z1;
                float tau = dz / dr;
                float ftau = std::abs(tau);
                if (ftau > 36.0) {
                  continue;
                }

                if (ftau < n1->m_minCutOnTau) {
                  continue;
                }
                if (ftau < n2->m_minCutOnTau) {
                  continue;
                }
                if (ftau > n1->m_maxCutOnTau) {
                  continue;
                }
                if (ftau > n2->m_maxCutOnTau) {
                  continue;
                }

                if (m_config.m_doubletFilterRZ) {
                  float z0 = z1 - r1 * tau;

                  if (z0 < min_z0 || z0 > max_z0) {
                    continue;
                  }

                  float zouter = z0 + m_config.maxOuterRadius * tau;

                  if (zouter < cut_zMinU || zouter > cut_zMaxU) {
                    continue;
                  }
                }

                float dx = n2->m_spGbts.SP->x() - x1;
                float dy = n2->m_spGbts.SP->y() - y1;

                float L2 = 1 / (dx * dx + dy * dy);

                float D =
                    (n2->m_spGbts.SP->y() * x1 - y1 * n2->m_spGbts.SP->x()) /
                    (r1 * r2);

                float kappa = D * D * L2;

                if (ftau < 4.0) {  // eta = 2.1
                  if (kappa > maxKappa_low_eta) {
                    continue;
                  }

                } else {
                  if (kappa > maxKappa_high_eta) {
                    continue;
                  }
                }

                // match edge candidate against edges incoming to n2

                float exp_eta = std::sqrt(1 + tau * tau) - tau;

                bool isGood =
                    n2->m_in.size() <=
                    2;  // we must have enough incoming edges to decide

                if (!isGood) {
                  float uat_1 = 1.0f / exp_eta;

                  for (const auto& n2_in_idx : n2->m_in) {
                    float tau2 = edgeStorage.at(n2_in_idx).m_p[0];
                    float tau_ratio = tau2 * uat_1 - 1.0f;

                    // bad match
                    if (std::abs(tau_ratio) > m_config.cut_tau_ratio_max) {
                      continue;
                    }

                    // good match found
                    isGood = true;
                    break;
                  }
                }
                if (!isGood) {
                  continue;  // no match found, skip creating [n1 <- n2] edge
                }

                float curv = D * std::sqrt(L2);  // signed curvature
                float dPhi2 = std::asin(curv * r2);
                float dPhi1 = std::asin(curv * r1);

                if (nEdges < m_config.MaxEdges) {
                  edgeStorage.emplace_back(n1.get(), n2, exp_eta, curv,
                                           phi1 + dPhi1, phi2 + dPhi2);

                  n1->addIn(nEdges);
                  n2->addOut(nEdges);

                  nEdges++;
                }
              }  // loop over n2 (outer) nodes
            }  // loop over n1 (inner) nodes
          }  // loop over source eta bins
        }  // loop over dst eta bins
      }  // loop over L2(L1) layers
    }  // loop over dst layers
  }  // loop over the stages of doublet making

  std::vector<const GbtsNode<external_spacepoint_t>*> vNodes;

  m_storage->getConnectingNodes(vNodes);

  if (vNodes.empty()) {
    ACTS_VERBOSE("No nodes");
    return;
  }

  int nNodes = vNodes.size();

  for (int nodeIdx = 0; nodeIdx < nNodes; nodeIdx++) {
    const GbtsNode<external_spacepoint_t>& pN = *vNodes.at(nodeIdx);

    std::vector<std::pair<float, int>> in_sort, out_sort;
    in_sort.resize(pN.m_in.size());
    out_sort.resize(pN.m_out.size());

    for (int inIdx = 0; inIdx < static_cast<int>(pN.m_in.size()); inIdx++) {
      int inEdgeIdx = pN.m_in.at(inIdx);
      Acts::GbtsEdge<external_spacepoint_t>* pS = &(edgeStorage.at(inEdgeIdx));
      in_sort[inIdx].second = inEdgeIdx;
      in_sort[inIdx].first = pS->m_p[0];
    }
    for (int outIdx = 0; outIdx < static_cast<int>(pN.m_out.size()); outIdx++) {
      int outEdgeIdx = pN.m_out.at(outIdx);
      Acts::GbtsEdge<external_spacepoint_t>* pS = &(edgeStorage.at(outEdgeIdx));
      out_sort[outIdx].second = outEdgeIdx;
      out_sort[outIdx].first = pS->m_p[0];
    }

    std::ranges::sort(in_sort);
    std::ranges::sort(out_sort);

    unsigned int last_out = 0;

    for (unsigned int in_idx = 0; in_idx < in_sort.size();
         in_idx++) {  // loop over incoming edges

      int inEdgeIdx = in_sort[in_idx].second;

      Acts::GbtsEdge<external_spacepoint_t>* pS = &(edgeStorage.at(inEdgeIdx));

      pS->m_nNei = 0;
      float tau1 = pS->m_p[0];
      float uat_1 = 1.0f / tau1;
      float curv1 = pS->m_p[1];
      float Phi1 = pS->m_p[2];

      for (unsigned int out_idx = last_out; out_idx < out_sort.size();
           out_idx++) {
        int outEdgeIdx = out_sort[out_idx].second;

        Acts::GbtsEdge<external_spacepoint_t>* pNS =
            &(edgeStorage.at(outEdgeIdx));

        float tau2 = pNS->m_p[0];
        float tau_ratio = tau2 * uat_1 - 1.0f;

        if (tau_ratio < -m_config.cut_tau_ratio_max) {
          last_out = out_idx;
          continue;
        }
        if (tau_ratio > m_config.cut_tau_ratio_max) {
          break;
        }

        float dPhi = pNS->m_p[3] - Phi1;

        if (dPhi < -std::numbers::pi_v<float>) {
          dPhi += static_cast<float>(2 * std::numbers::pi);
        } else if (dPhi > std::numbers::pi_v<float>) {
          dPhi -= static_cast<float>(2 * std::numbers::pi);
        }

        if (dPhi < -m_config.cut_dphi_max || dPhi > m_config.cut_dphi_max) {
          continue;
        }

        float curv2 = pNS->m_p[1];
        float dcurv = curv2 - curv1;

        if (dcurv < -m_config.cut_dcurv_max || dcurv > m_config.cut_dcurv_max) {
          continue;
        }

        pS->m_vNei[pS->m_nNei++] = outEdgeIdx;
        if (pS->m_nNei >= N_SEG_CONNS) {
          break;
        }
      }
    }
  }

  const int maxIter = 15;

  int maxLevel = 0;

  int iter = 0;

  std::vector<Acts::GbtsEdge<external_spacepoint_t>*> v_old;

  for (int edgeIndex = 0; edgeIndex < nEdges; edgeIndex++) {
    Acts::GbtsEdge<external_spacepoint_t>* pS = &(edgeStorage.at(edgeIndex));
    if (pS->m_nNei == 0) {
      continue;
    }
    v_old.push_back(pS);  // TO-DO: increment level for segments as they already
                          // have at least one neighbour
  }

  for (; iter < maxIter; iter++) {
    // generate proposals
    std::vector<Acts::GbtsEdge<external_spacepoint_t>*> v_new;
    v_new.clear();

    for (auto pS : v_old) {
      int next_level = pS->m_level;

      for (int nIdx = 0; nIdx < pS->m_nNei; nIdx++) {
        unsigned int nextEdgeIdx = pS->m_vNei[nIdx];

        Acts::GbtsEdge<external_spacepoint_t>* pN =
            &(edgeStorage.at(nextEdgeIdx));

        if (pS->m_level == pN->m_level) {
          next_level = pS->m_level + 1;
          v_new.push_back(pS);
          break;
        }
      }

      pS->m_next = next_level;  // proposal
    }
    // update

    int nChanges = 0;

    for (auto pS : v_new) {
      if (pS->m_next != pS->m_level) {
        nChanges++;
        pS->m_level = pS->m_next;
        if (maxLevel < pS->m_level) {
          maxLevel = pS->m_level;
        }
      }
    }

    if (nChanges == 0) {
      break;
    }

    v_old = std::move(v_new);
    v_new.clear();
  }

  int minLevel = 3;  // a triplet + 2 confirmation

  std::vector<Acts::GbtsEdge<external_spacepoint_t>*> vSeeds;

  vSeeds.reserve(m_config.MaxEdges / 2);

  for (int edgeIndex = 0; edgeIndex < nEdges; edgeIndex++) {
    Acts::GbtsEdge<external_spacepoint_t>* pS = &(edgeStorage.at(edgeIndex));

    if (pS->m_level < minLevel) {
      continue;
    }

    vSeeds.push_back(pS);
  }

  m_triplets.clear();

  std::ranges::sort(
      vSeeds, typename Acts::GbtsEdge<external_spacepoint_t>::CompareLevel());

  if (vSeeds.empty()) {
    return;
  }

  // backtracking

  GbtsTrackingFilter<external_spacepoint_t> tFilter(
      m_config.m_layerGeometry, edgeStorage,
      logger().cloneWithSuffix("GbtsFilter"));

  for (auto pS : vSeeds) {
    if (pS->m_level == -1) {
      continue;
    }

    GbtsEdgeState<external_spacepoint_t> rs(false);

    tFilter.followTrack(pS, rs);

    if (!rs.m_initialized) {
      continue;
    }

    if (static_cast<int>(rs.m_vs.size()) < minLevel) {
      continue;
    }

    std::vector<const GbtsSP<external_spacepoint_t>*> vSP;

    for (typename std::vector<
             Acts::GbtsEdge<external_spacepoint_t>*>::reverse_iterator sIt =
             rs.m_vs.rbegin();
         sIt != rs.m_vs.rend(); ++sIt) {
      (*sIt)->m_level = -1;  // mark as collected

      if (sIt == rs.m_vs.rbegin()) {
        vSP.push_back(&(*sIt)->m_n1->m_spGbts);
      }
      vSP.push_back(&(*sIt)->m_n2->m_spGbts);
    }

    if (vSP.size() < 3) {
      continue;
    }

    // making triplets

    unsigned int nTriplets = 0;

    std::vector<TrigInDetTriplet<external_spacepoint_t>> output;

    for (unsigned int idx_m = 1; idx_m < vSP.size() - 1; idx_m++) {
      const GbtsSP<external_spacepoint_t>& spM = *vSP.at(idx_m);
      const double pS_r = spM.r();
      const double pS_x = spM.SP->x();
      const double pS_y = spM.SP->y();
      const double cosA = pS_x / pS_r;
      const double sinA = pS_y / pS_r;

      for (unsigned int idx_o = idx_m + 1; idx_o < vSP.size(); idx_o++) {
        const GbtsSP<external_spacepoint_t>& spO = *vSP.at(idx_o);

        double dx = spO.SP->x() - pS_x;
        double dy = spO.SP->y() - pS_y;
        double R2inv = 1.0 / (dx * dx + dy * dy);
        double xn = dx * cosA + dy * sinA;
        double yn = -dx * sinA + dy * cosA;

        const double uo = xn * R2inv;
        const double vo = yn * R2inv;

        for (unsigned int idx_i = 0; idx_i < idx_m; idx_i++) {
          const GbtsSP<external_spacepoint_t>& spI = *vSP.at(idx_i);

          dx = spI.SP->x() - pS_x;
          dy = spI.SP->y() - pS_y;
          R2inv = 1.0 / (dx * dx + dy * dy);

          xn = dx * cosA + dy * sinA;
          yn = -dx * sinA + dy * cosA;

          const double ui = xn * R2inv;
          const double vi = yn * R2inv;

          // 1. pT estimate

          const double du = uo - ui;
          if (du == 0.0) {
            continue;
          }
          const double A = (vo - vi) / du;
          const double B = vi - A * ui;

          if ((1 + A * A) < (B * B) * m_minR_squ) {
            continue;
          }

          // 2. d0 cut

          const double fabs_d0 = std::abs(pS_r * (B * pS_r - A));

          if (fabs_d0 > m_config.m_tripletD0Max) {
            continue;
          }

          // 3. phi0 cut

          // if (!roi.isFullscan()) {
          //   const double uc = 2 * B * pS_r - A;
          //   // const double phi0 = std::atan2(sinA - uc * cosA, cosA + uc *
          //   sinA);
          //   //           if ( !RoiUtil::containsPhi( *roiDescriptor, phi0 ) )
          //   {
          //   //             continue;
          //   //           }
          // }

          // 4. add new triplet

          const double Q = fabs_d0 * fabs_d0;

          output.emplace_back(spI, spM, spO, Q);

          nTriplets++;

          if (nTriplets >= m_config.m_maxTripletBufferLength) {
            break;
          }
        }
        if (nTriplets >= m_config.m_maxTripletBufferLength) {
          break;
        }
      }
      if (nTriplets >= m_config.m_maxTripletBufferLength) {
        break;
      }
    }

    if (output.empty()) {
      continue;
    }

    vTracks.emplace_back(vSP, output);
  }
}

template <typename external_spacepoint_t>
template <typename output_container_t>
void SeedFinderGbts<external_spacepoint_t>::createSeeds(
    const Acts::RoiDescriptor& roi,
    const Acts::GbtsGeometry<external_spacepoint_t>& gbtsGeo,
    output_container_t& out_cont) {
  ACTS_VERBOSE("Creating seeds");
  std::vector<GbtsTrigTracklet<external_spacepoint_t>>
      vTracks;  // make empty vector

  vTracks.reserve(5000);

  runGbts_TrackFinder(vTracks, roi, gbtsGeo);  // returns filled vector

  if (vTracks.empty()) {
    return;
  }

  m_triplets.clear();  // member of class , saying not declared, maybe public?

  for (auto& track : vTracks) {
    for (auto& seed : track.m_seeds) {  // access member of GbtsTrigTracklet

      float newQ = seed.Q();  // function of TrigInDetTriplet
      if (m_config.m_LRTmode) {
        // In LRT mode penalize pixels in Triplets
        if (seed.s1().isPixel()) {
          newQ += 1000;  // functions of TrigSiSpacePointBase
        }
        if (seed.s2().isPixel()) {
          newQ += 1000;
        }
        if (seed.s3().isPixel()) {
          newQ += 1000;
        }
      } else {
        // In normal (non LRT) mode penalise SSS by 1000, PSS (if enabled) and
        // PPS by 10000
        if (seed.s3().isStrip()) {
          newQ += seed.s1().isStrip() ? 1000.0 : 10000.0;
        }
      }
      seed.Q(newQ);
      m_triplets.emplace_back(seed);
    }
  }
  vTracks.clear();

  for (auto& triplet : m_triplets) {
    const external_spacepoint_t* S1 =
        triplet.s1().SP;  // triplet-> GbtsSP-> simspacepoint
    const external_spacepoint_t* S2 = triplet.s2().SP;
    const external_spacepoint_t* S3 = triplet.s3().SP;

    // input to seed
    float Vertex = 0;
    float Quality = triplet.Q();
    // make a new seed, add to vector of seeds
    out_cont.emplace_back(*S1, *S2, *S3);
    out_cont.back().setVertexZ(Vertex);
    out_cont.back().setQuality(Quality);
  }
}

// outer called in alg
template <typename external_spacepoint_t>
std::vector<Seed<external_spacepoint_t>>
SeedFinderGbts<external_spacepoint_t>::createSeeds(
    const Acts::RoiDescriptor& roi,
    const Acts::GbtsGeometry<external_spacepoint_t>& gbtsGeo) {
  std::vector<seed_t> r;
  createSeeds(roi, gbtsGeo, r);
  return r;
}

}  // namespace Acts
