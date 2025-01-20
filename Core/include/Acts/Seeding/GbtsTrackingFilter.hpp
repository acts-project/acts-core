// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// TODO: update to C++17 style
#include "Acts/Seeding/GbtsDataStorage.hpp"  //includes geo which has trigindetsilayer, may move this to trigbase
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <list>
#include <vector>

template <typename external_spacepoint_t>
struct GbtsEdgeState {
 public:
  struct Compare {
    bool operator()(const struct GbtsEdgeState* s1,
                    const struct GbtsEdgeState* s2) {
      return s1->m_J > s2->m_J;
    }
  };

  GbtsEdgeState() = default;

  GbtsEdgeState(bool f) : m_initialized(f) {}

  void initialize(Acts::GbtsEdge<external_spacepoint_t>* pS) {
    m_initialized = true;

    m_J = 0.0;
    m_vs.clear();

    // n2->n1

    float dx = pS->m_n1->m_spGbts.SP->x() - pS->m_n2->m_spGbts.SP->x();
    float dy = pS->m_n1->m_spGbts.SP->y() - pS->m_n2->m_spGbts.SP->y();
    float L = std::sqrt(dx * dx + dy * dy);

    m_s = dy / L;
    m_c = dx / L;

    // transform for extrapolation and update
    //  x' =  x*m_c + y*m_s
    //  y' = -x*m_s + y*m_c

    m_refY = pS->m_n2->m_spGbts.r();
    m_refX =
        pS->m_n2->m_spGbts.SP->x() * m_c + pS->m_n2->m_spGbts.SP->y() * m_s;

    // X-state: y, dy/dx, d2y/dx2

    m_X[0] =
        -pS->m_n2->m_spGbts.SP->x() * m_s + pS->m_n2->m_spGbts.SP->y() * m_c;
    m_X[1] = 0.0;
    m_X[2] = 0.0;

    // Y-state: z, dz/dr

    m_Y[0] = pS->m_n2->m_spGbts.SP->z();
    m_Y[1] = (pS->m_n1->m_spGbts.SP->z() - pS->m_n2->m_spGbts.SP->z()) /
             (pS->m_n1->m_spGbts.r() - pS->m_n2->m_spGbts.r());

    memset(&m_Cx[0][0], 0, sizeof(m_Cx));
    memset(&m_Cy[0][0], 0, sizeof(m_Cy));

    m_Cx[0][0] = 0.25;
    m_Cx[1][1] = 0.001;
    m_Cx[2][2] = 0.001;

    m_Cy[0][0] = 1.5;
    m_Cy[1][1] = 0.001;
  }

  void clone(const struct GbtsEdgeState& st) {
    memcpy(&m_X[0], &st.m_X[0], sizeof(m_X));
    memcpy(&m_Y[0], &st.m_Y[0], sizeof(m_Y));
    memcpy(&m_Cx[0][0], &st.m_Cx[0][0], sizeof(m_Cx));
    memcpy(&m_Cy[0][0], &st.m_Cy[0][0], sizeof(m_Cy));
    m_refX = st.m_refX;
    m_refY = st.m_refY;
    m_c = st.m_c;
    m_s = st.m_s;
    m_J = st.m_J;
    m_vs.clear();
    m_vs.reserve(st.m_vs.size());
    std::copy(st.m_vs.begin(), st.m_vs.end(), std::back_inserter(m_vs));

    m_initialized = true;
  }

  float m_J{};

  std::vector<Acts::GbtsEdge<external_spacepoint_t>*> m_vs;

  float m_X[3]{}, m_Y[2]{}, m_Cx[3][3]{}, m_Cy[2][2]{};
  float m_refX{}, m_refY{}, m_c{}, m_s{};

  bool m_initialized{false};
};

#define MAX_EDGE_STATE 2500

template <typename external_spacepoint_t>
class GbtsTrackingFilter {
 public:
  GbtsTrackingFilter(const std::vector<Acts::TrigInDetSiLayer>& g,
                     std::vector<Acts::GbtsEdge<external_spacepoint_t>>& sb,
                     std::unique_ptr<const Acts::Logger> logger =
                         Acts::getDefaultLogger("Filter",
                                                Acts::Logging::Level::INFO))
      : m_geo(g), m_segStore(sb), m_logger(std::move(logger)) {}

  void followTrack(Acts::GbtsEdge<external_spacepoint_t>* pS,
                   GbtsEdgeState<external_spacepoint_t>& output) {
    if (pS->m_level == -1) {
      return;  // already collected
    }
    m_globalStateCounter = 0;

    // create track state

    GbtsEdgeState<external_spacepoint_t>* pInitState =
        &m_stateStore[m_globalStateCounter++];

    pInitState->initialize(pS);

    m_stateVec.clear();

    // recursive branching and propagation

    propagate(pS, *pInitState);

    if (m_stateVec.empty()) {
      return;
    }
    std::ranges::sort(m_stateVec,
                      typename GbtsEdgeState<external_spacepoint_t>::Compare());

    GbtsEdgeState<external_spacepoint_t>* best = (*m_stateVec.begin());

    output.clone(*best);

    m_globalStateCounter = 0;
  }

 protected:
  void propagate(Acts::GbtsEdge<external_spacepoint_t>* pS,
                 GbtsEdgeState<external_spacepoint_t>& ts) {
    if (m_globalStateCounter >= MAX_EDGE_STATE) {
      return;
    }
    GbtsEdgeState<external_spacepoint_t>* p_new_ts =
        &m_stateStore[m_globalStateCounter++];

    GbtsEdgeState<external_spacepoint_t>& new_ts = *p_new_ts;
    new_ts.clone(ts);

    new_ts.m_vs.push_back(pS);

    bool accepted = update(pS, new_ts);  // update using n1 of the segment

    if (!accepted) {
      return;  // stop further propagation
    }
    int level = pS->m_level;

    std::list<Acts::GbtsEdge<external_spacepoint_t>*> lCont;

    for (int nIdx = 0; nIdx < pS->m_nNei;
         nIdx++) {  // loop over the neighbours of this segment
      unsigned int nextSegmentIdx = pS->m_vNei[nIdx];

      Acts::GbtsEdge<external_spacepoint_t>* pN =
          &(m_segStore.at(nextSegmentIdx));

      if (pN->m_level == -1) {
        continue;  // already collected
      }
      if (pN->m_level == level - 1) {
        lCont.push_back(pN);
      }
    }
    if (lCont.empty()) {  // the end of chain

      // store in the vector
      if (m_globalStateCounter < MAX_EDGE_STATE) {
        if (m_stateVec.empty()) {  // add the first segment state
          GbtsEdgeState<external_spacepoint_t>* p =
              &m_stateStore[m_globalStateCounter++];
          p->clone(new_ts);
          m_stateVec.push_back(p);
        } else {  // compare with the best and add
          float best_so_far = (*m_stateVec.begin())->m_J;
          if (new_ts.m_J > best_so_far) {
            GbtsEdgeState<external_spacepoint_t>* p =
                &m_stateStore[m_globalStateCounter++];
            p->clone(new_ts);
            m_stateVec.push_back(p);
          }
        }
      }
    } else {  // branching
      int nBranches = 0;
      for (typename std::list<Acts::GbtsEdge<external_spacepoint_t>*>::iterator
               sIt = lCont.begin();
           sIt != lCont.end(); ++sIt, nBranches++) {
        propagate((*sIt), new_ts);  // recursive call
      }
    }
  }

  bool update(Acts::GbtsEdge<external_spacepoint_t>* pS,
              GbtsEdgeState<external_spacepoint_t>& ts) {
    const float sigma_t = 0.0003;
    const float sigma_w = 0.00009;

    const float sigmaMS = 0.016;

    const float sigma_x = 0.25;  // was 0.22
    const float sigma_y = 2.5;   // was 1.7

    const float weight_x = 0.5;
    const float weight_y = 0.5;

    const float maxDChi2_x = 60.0;  // was 35.0;
    const float maxDChi2_y = 60.0;  // was 31.0;

    const float add_hit = 14.0;

    if (ts.m_Cx[2][2] < 0.0 || ts.m_Cx[1][1] < 0.0 || ts.m_Cx[0][0] < 0.0) {
      ACTS_WARNING("Negative covariance detected in X components: "
                   << "cov[2][2]=" << ts.m_Cx[2][2] << " cov[1][1]="
                   << ts.m_Cx[1][1] << " cov[0][0]=" << ts.m_Cx[0][0]);
    }

    if (ts.m_Cy[1][1] < 0.0 || ts.m_Cy[0][0] < 0.0) {
      ACTS_WARNING("Negative covariance detected in Y components: "
                   << "cov[1][1]=" << ts.m_Cy[1][1]
                   << " cov[0][0]=" << ts.m_Cy[0][0]);
    }

    // add ms.

    ts.m_Cx[2][2] += sigma_w * sigma_w;
    ts.m_Cx[1][1] += sigma_t * sigma_t;

    int type1 = getLayerType(pS->m_n1->m_spGbts.combined_ID);

    float t2 = type1 == 0 ? 1.0 + ts.m_Y[1] * ts.m_Y[1]
                          : 1.0 + 1.0 / (ts.m_Y[1] * ts.m_Y[1]);
    float s1 = sigmaMS * t2;
    float s2 = s1 * s1;

    s2 *= std::sqrt(t2);

    ts.m_Cy[1][1] += s2;

    // extrapolation

    float X[3], Y[2];
    float Cx[3][3], Cy[2][2];

    float refX{}, refY{}, mx{}, my{};

    float x{}, y{}, z{}, r{};

    x = pS->m_n1->m_spGbts.SP->x();
    y = pS->m_n1->m_spGbts.SP->y();
    z = pS->m_n1->m_spGbts.SP->z();
    r = pS->m_n1->m_spGbts.r();

    refX = x * ts.m_c + y * ts.m_s;
    mx = -x * ts.m_s + y * ts.m_c;  // measured X[0]
    refY = r;
    my = z;  // measured Y[0]

    float A = refX - ts.m_refX;
    float B = 0.5 * A * A;
    float dr = refY - ts.m_refY;

    X[0] = ts.m_X[0] + ts.m_X[1] * A + ts.m_X[2] * B;
    X[1] = ts.m_X[1] + ts.m_X[2] * A;
    X[2] = ts.m_X[2];

    Cx[0][0] = ts.m_Cx[0][0] + 2 * ts.m_Cx[0][1] * A + 2 * ts.m_Cx[0][2] * B +
               A * A * ts.m_Cx[1][1] + 2 * A * B * ts.m_Cx[1][2] +
               B * B * ts.m_Cx[2][2];
    Cx[0][1] = Cx[1][0] = ts.m_Cx[0][1] + ts.m_Cx[1][1] * A +
                          ts.m_Cx[1][2] * B + ts.m_Cx[0][2] * A +
                          A * A * ts.m_Cx[1][2] + A * B * ts.m_Cx[2][2];
    Cx[0][2] = Cx[2][0] = ts.m_Cx[0][2] + ts.m_Cx[1][2] * A + ts.m_Cx[2][2] * B;

    Cx[1][1] = ts.m_Cx[1][1] + 2 * A * ts.m_Cx[1][2] + A * A * ts.m_Cx[2][2];
    Cx[1][2] = Cx[2][1] = ts.m_Cx[1][2] + ts.m_Cx[2][2] * A;

    Cx[2][2] = ts.m_Cx[2][2];

    Y[0] = ts.m_Y[0] + ts.m_Y[1] * dr;
    Y[1] = ts.m_Y[1];

    Cy[0][0] = ts.m_Cy[0][0] + 2 * ts.m_Cy[0][1] * dr + dr * dr * ts.m_Cy[1][1];
    Cy[0][1] = Cy[1][0] = ts.m_Cy[0][1] + dr * ts.m_Cy[1][1];
    Cy[1][1] = ts.m_Cy[1][1];

    float resid_x = mx - X[0];
    float resid_y = my - Y[0];

    float CHx[3] = {Cx[0][0], Cx[0][1], Cx[0][2]};
    float CHy[2] = {Cy[0][0], Cy[0][1]};

    float sigma_rz = 0.0;

    int type = getLayerType(pS->m_n1->m_spGbts.combined_ID);

    if (type == 0) {  // barrel TODO: split into barrel Pixel and barrel SCT
      sigma_rz = sigma_y * sigma_y;
    } else {
      sigma_rz = sigma_y * ts.m_Y[1];
      sigma_rz = sigma_rz * sigma_rz;
    }

    float Dx = 1.0 / (Cx[0][0] + sigma_x * sigma_x);

    float Dy = 1.0 / (Cy[0][0] + sigma_rz);

    float dchi2_x = resid_x * resid_x * Dx;
    float dchi2_y = resid_y * resid_y * Dy;

    if (dchi2_x > maxDChi2_x || dchi2_y > maxDChi2_y) {
      return false;
    }

    ts.m_J += add_hit - dchi2_x * weight_x - dchi2_y * weight_y;

    // state update
    float Kx[3] = {Dx * Cx[0][0], Dx * Cx[0][1], Dx * Cx[0][2]};
    float Ky[2] = {Dy * Cy[0][0], Dy * Cy[0][1]};

    for (int i = 0; i < 3; i++) {
      ts.m_X[i] = X[i] + Kx[i] * resid_x;
    }
    for (int i = 0; i < 2; i++) {
      ts.m_Y[i] = Y[i] + Ky[i] * resid_y;
    }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        ts.m_Cx[i][j] = Cx[i][j] - Kx[i] * CHx[j];
      }
    }

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        ts.m_Cy[i][j] = Cy[i][j] - Ky[i] * CHy[j];
      }
    }

    ts.m_refX = refX;
    ts.m_refY = refY;

    return true;
  }

  int getLayerType(int l) {
    auto iterator = find_if(m_geo.begin(), m_geo.end(), [l](auto n) {
      return n.m_subdet == l;
    });  // iterator to vector member with this id
    int index = std::distance(m_geo.begin(), iterator);

    return m_geo.at(index).m_type;  // needs to be 0, 2, or -2
  }

  const std::vector<Acts::TrigInDetSiLayer>& m_geo;

  std::vector<Acts::GbtsEdge<external_spacepoint_t>>& m_segStore;

  std::vector<GbtsEdgeState<external_spacepoint_t>*> m_stateVec;

  GbtsEdgeState<external_spacepoint_t> m_stateStore[MAX_EDGE_STATE];

  int m_globalStateCounter{0};

  const Acts::Logger& logger() const { return *m_logger; }
  std::unique_ptr<const Acts::Logger> m_logger{nullptr};
};
