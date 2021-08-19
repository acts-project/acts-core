// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedingOrthogonalAlgorithm.hpp"

#include "Acts/Utilities/KDTree.hpp"
#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SeedfinderUtils.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>

using SP = Acts::InternalSpacePoint<ActsExamples::SimSpacePoint>;
using Protoseed = std::tuple<const SP *, const SP *, const SP *>;

namespace {
  void transformCoordinates(
      std::vector<const SP*>& vec,
      const SP& spM, bool bottom,
      std::vector<Acts::LinCircle>& linCircleVec) {
    float xM = spM.x();
    float yM = spM.y();
    float zM = spM.z();
    float rM = spM.radius();
    float varianceZM = spM.varianceZ();
    float varianceRM = spM.varianceR();
    float cosPhiM = xM / rM;
    float sinPhiM = yM / rM;
    for (auto sp : vec) {
      float deltaX = sp->x() - xM;
      float deltaY = sp->y() - yM;
      float deltaZ = sp->z() - zM;
      // calculate projection fraction of spM->sp vector pointing in same
      // direction as
      // vector origin->spM (x) and projection fraction of spM->sp vector pointing
      // orthogonal to origin->spM (y)
      float x = deltaX * cosPhiM + deltaY * sinPhiM;
      float y = deltaY * cosPhiM - deltaX * sinPhiM;
      // 1/(length of M -> SP)
      float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
      float iDeltaR = std::sqrt(iDeltaR2);
      //
      int bottomFactor = 1 * (int(!bottom)) - 1 * (int(bottom));
      // cot_theta = (deltaZ/deltaR)
      float cot_theta = deltaZ * iDeltaR * bottomFactor;
      // VERY frequent (SP^3) access
      Acts::LinCircle l;
      l.cotTheta = cot_theta;
      // location on z-axis of this SP-duplet
      l.Zo = zM - rM * cot_theta;
      l.iDeltaR = iDeltaR;
      // transformation of circle equation (x,y) into linear equation (u,v)
      // x^2 + y^2 - 2x_0*x - 2y_0*y = 0
      // is transformed into
      // 1 - 2x_0*u - 2y_0*v = 0
      // using the following m_U and m_V
      // (u = A + B*v); A and B are created later on
      l.U = x * iDeltaR2;
      l.V = y * iDeltaR2;
      // error term for sp-pair without correlation of middle space point
      l.Er = ((varianceZM + sp->varianceZ()) +
              (cot_theta * cot_theta) * (varianceRM + sp->varianceR())) *
            iDeltaR2;
      linCircleVec.push_back(l);
    }
  }
}

ActsExamples::SeedingOrthogonalAlgorithm::SeedingOrthogonalAlgorithm(
    ActsExamples::SeedingOrthogonalAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("SeedingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point input collections");
  }
  for (const auto& i : m_cfg.inputSpacePoints) {
    if (i.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing proto tracks output collection");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds output collection");
  }

  m_gridCfg.bFieldInZ = m_cfg.bFieldInZ;
  m_gridCfg.minPt = m_cfg.minPt;
  m_gridCfg.rMax = m_cfg.rMax;
  m_gridCfg.zMax = m_cfg.zMax;
  m_gridCfg.zMin = m_cfg.zMin;
  m_gridCfg.deltaRMax = m_cfg.deltaRMax;
  m_gridCfg.cotThetaMax = m_cfg.cotThetaMax;

  // construct seed filter
  Acts::SeedFilterConfig filterCfg;
  filterCfg.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
  m_finderCfg.seedFilter = std::make_unique<Acts::SeedFilter<SimSpacePoint>>(
      Acts::SeedFilter<SimSpacePoint>(filterCfg));

  m_finderCfg.rMax = m_cfg.rMax;
  m_finderCfg.deltaRMin = m_cfg.deltaRMin;
  m_finderCfg.deltaRMax = m_cfg.deltaRMax;
  m_finderCfg.collisionRegionMin = m_cfg.collisionRegionMin;
  m_finderCfg.collisionRegionMax = m_cfg.collisionRegionMax;
  m_finderCfg.zMin = m_cfg.zMin;
  m_finderCfg.zMax = m_cfg.zMax;
  m_finderCfg.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
  m_finderCfg.cotThetaMax = m_cfg.cotThetaMax;
  m_finderCfg.sigmaScattering = m_cfg.sigmaScattering;
  m_finderCfg.radLengthPerSeed = m_cfg.radLengthPerSeed;
  m_finderCfg.minPt = m_cfg.minPt;
  m_finderCfg.bFieldInZ = m_cfg.bFieldInZ;
  m_finderCfg.beamPos = Acts::Vector2(m_cfg.beamPosX, m_cfg.beamPosY);
  m_finderCfg.impactMax = m_cfg.impactMax;

  auto _config = m_finderCfg.toInternalUnits();
  m_finderCfg = _config;

  // calculation of scattering using the highland formula
  // convert pT to p once theta angle is known
  m_finderCfg.highland = 13.6 * std::sqrt(m_finderCfg.radLengthPerSeed) *
                      (1 + 0.038 * std::log(m_finderCfg.radLengthPerSeed));
  float maxScatteringAngle = m_finderCfg.highland / m_finderCfg.minPt;
  m_finderCfg.maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;
  // helix radius in homogeneous magnetic field. Units are Kilotesla, MeV and
  // millimeter
  // TODO: change using ACTS units
  m_finderCfg.pTPerHelixRadius = 300. * m_finderCfg.bFieldInZ;
  m_finderCfg.minHelixDiameter2 =
      std::pow(m_finderCfg.minPt * 2 / m_finderCfg.pTPerHelixRadius, 2);
  std::cout << "Min helix radius = " << (0.5 * std::sqrt(m_finderCfg.minHelixDiameter2)) << std::endl;
  m_finderCfg.pT2perRadius =
      std::pow(m_finderCfg.highland / m_finderCfg.pTPerHelixRadius, 2);
}

std::optional<std::pair<double, double>> ActsExamples::SeedingOrthogonalAlgorithm::validTriple(
  const Acts::InternalSpacePoint<ActsExamples::SimSpacePoint> & b,
  const Acts::InternalSpacePoint<ActsExamples::SimSpacePoint> & m,
  const Acts::InternalSpacePoint<ActsExamples::SimSpacePoint> & t
) const {
  const double maxScatteringAngle = std::sqrt(m_finderCfg.maxScatteringAngle2);
  const double pTperRadius = std::sqrt(m_finderCfg.pT2perRadius);
  std::complex<double> rotPhiM = std::polar(1.0, static_cast<double>(std::asin(-m.y() / m.radius())));
  std::complex<double> mc(m.x(), m.y());

  Acts::LinCircle lb = lcBotMid(b, m);

  std::complex<double> bc(b.x(), b.y());

  std::complex<double> qc = (mc - bc);
  std::complex<double> tf = std::complex(1.0, 0.0) / qc;

  float csc = std::sqrt(1. + lb.cotTheta * lb.cotTheta);

  std::complex<double> tc(t.x(), t.y());

  std::complex<double> w = (tc - bc) * tf;
  std::complex<double> mp = (w - std::norm(w)) / (w - std::conj(w));
  std::complex<double> M = qc * mp + bc;
  double norm = std::norm(bc - M);

  if (4 * norm < m_finderCfg.minHelixDiameter2) {
    return std::nullopt;
  }

  double radius = std::sqrt(norm);

  std::complex<double> I = M - M * radius / std::abs(M);

  if (std::norm(I) > std::pow(m_finderCfg.impactMax, 2)) {
    return std::nullopt;
  }

  std::complex<double> deltalb = bc - mc;
  std::complex<double> deltalt = tc - mc;

  double ltCotTheta = (t.z() - m.z()) / std::abs(deltalt);

  double ltEr = ((m.varianceZ() + t.varianceZ()) +
    (ltCotTheta * ltCotTheta) * (m.varianceR() + t.varianceR())) / std::norm(deltalt);




  // add errors of spB-spM and spM-spT pairs and add the correlation term
  // for errors on spM
  double error2 = ltEr + lb.Er +
                  2 * (lb.cotTheta * ltCotTheta * m.varianceR() + m.varianceZ()) *
                      lb.iDeltaR / std::abs(tc - mc);

  double deltaCotTheta = std::abs(lb.cotTheta - ltCotTheta);

  // if the error is larger than the difference in theta, no need to
  // compare with scattering
  if (deltaCotTheta * deltaCotTheta > error2) {
    double pT2scatter;

    if (m_finderCfg.pTPerHelixRadius * radius > m_finderCfg.maxPtScattering) {
      pT2scatter = m_finderCfg.highland / m_finderCfg.maxPtScattering;
    } else {
      pT2scatter = pTperRadius / radius;
    }

    double max_scat = std::min(maxScatteringAngle, pT2scatter);

    if (deltaCotTheta - std::sqrt(error2) > m_finderCfg.sigmaScattering * csc * max_scat) {
      return std::nullopt;
    }
  }

  std::complex<double> dlb = deltalb / std::complex(std::norm(deltalb), 0.0);
  std::complex<double> dlt = deltalt / std::complex(std::norm(deltalt), 0.0);
  std::complex<double> ltc = dlt * rotPhiM;
  std::complex<double> lbc = dlb * rotPhiM;
  std::complex<double> myd = ltc - lbc;
  double secMyd = std::abs(myd) / myd.real();
  double tanMyd = myd.imag() / myd.real();
  double sectanMyd = secMyd * tanMyd;

  return std::pair<double, double>(std::abs(I), lbc.imag() - sectanMyd * lbc.real());
}

std::vector<std::pair<float, std::unique_ptr<const Acts::InternalSeed<ActsExamples::SimSpacePoint>>>>
ActsExamples::SeedingOrthogonalAlgorithm::helper(
  const SP & middle,
  std::vector<const SP *> & bottom,
  std::vector<const SP *> & top
) const {
  float rM = middle.radius();
  float varianceRM = middle.varianceR();
  float varianceZM = middle.varianceZ();

  std::vector<const SP *> top_valid;
  std::vector<float> curvatures;
  std::vector<float> impactParameters;

  std::vector<std::pair<float, std::unique_ptr<const Acts::InternalSeed<ActsExamples::SimSpacePoint>>>> protoseeds;
  // contains parameters required to calculate circle with linear equation
  // ...for bottom-middle
  std::vector<Acts::LinCircle> linCircleBottom;
  linCircleBottom.reserve(bottom.size());
  // ...for middle-top
  std::vector<Acts::LinCircle> linCircleTop;
  linCircleTop.reserve(top.size());
  transformCoordinates(bottom, middle, true, linCircleBottom);
  transformCoordinates(top, middle, false, linCircleTop);

  size_t numBotSP = bottom.size();
  size_t numTopSP = top.size();

  for (size_t b = 0; b < numBotSP; b++) {
    auto lb = linCircleBottom[b];
    float Zob = lb.Zo;
    float cotThetaB = lb.cotTheta;
    float Vb = lb.V;
    float Ub = lb.U;
    float ErB = lb.Er;
    float iDeltaRB = lb.iDeltaR;

    // 1+(cot^2(theta)) = 1/sin^2(theta)
    float iSinTheta2 = (1. + cotThetaB * cotThetaB);
    // calculate max scattering for min momentum at the seed's theta angle
    // scaling scatteringAngle^2 by sin^2(theta) to convert pT^2 to p^2
    // accurate would be taking 1/atan(thetaBottom)-1/atan(thetaTop) <
    // scattering
    // but to avoid trig functions we approximate cot by scaling by
    // 1/sin^4(theta)
    // resolving with pT to p scaling --> only divide by sin^2(theta)
    // max approximation error for allowed scattering angles of 0.04 rad at
    // eta=infinity: ~8.5%
    float scatteringInRegion2 = m_finderCfg.maxScatteringAngle2 * iSinTheta2;
    // multiply the squared sigma onto the squared scattering
    scatteringInRegion2 *=
        m_finderCfg.sigmaScattering * m_finderCfg.sigmaScattering;

    // clear all vectors used in each inner for loop
    top_valid.clear();
    curvatures.clear();
    impactParameters.clear();
    for (size_t t = 0; t < numTopSP; t++) {
      auto lt = linCircleTop[t];

      // add errors of spB-spM and spM-spT pairs and add the correlation term
      // for errors on spM
      float error2 = lt.Er + ErB +
                      2 * (cotThetaB * lt.cotTheta * varianceRM + varianceZM) *
                          iDeltaRB * lt.iDeltaR;

      float deltaCotTheta = cotThetaB - lt.cotTheta;
      float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;
      float error;
      float dCotThetaMinusError2;
      // if the error is larger than the difference in theta, no need to
      // compare with scattering
      if (deltaCotTheta2 - error2 > 0) {
        deltaCotTheta = std::abs(deltaCotTheta);
        // if deltaTheta larger than the scattering for the lower pT cut, skip
        error = std::sqrt(error2);
        dCotThetaMinusError2 =
            deltaCotTheta2 + error2 - 2 * deltaCotTheta * error;
        // avoid taking root of scatteringInRegion
        // if left side of ">" is positive, both sides of unequality can be
        // squared
        // (scattering is always positive)

        if (dCotThetaMinusError2 > scatteringInRegion2) {
          continue;
        }
      }

      // protects against division by 0
      float dU = lt.U - Ub;
      if (dU == 0.) {
        continue;
      }
      // A and B are evaluated as a function of the circumference parameters
      // x_0 and y_0
      float A = (lt.V - Vb) / dU;
      float S2 = 1. + A * A;
      float B = Vb - A * Ub;
      float B2 = B * B;
      // sqrt(S2)/B = 2 * helixradius
      // calculated radius must not be smaller than minimum radius
      if (S2 < B2 * m_finderCfg.minHelixDiameter2) {
        continue;
      }
      // 1/helixradius: (B/sqrt(S2))*2 (we leave everything squared)
      float iHelixDiameter2 = B2 / S2;
      // calculate scattering for p(T) calculated from seed curvature
      float pT2scatter = 4 * iHelixDiameter2 * m_finderCfg.pT2perRadius;
      // if pT > maxPtScattering, calculate allowed scattering angle using
      // maxPtScattering instead of pt.
      float pT = m_finderCfg.pTPerHelixRadius * std::sqrt(S2 / B2) / 2.;
      if (pT > m_finderCfg.maxPtScattering) {
        float pTscatter = m_finderCfg.highland / m_finderCfg.maxPtScattering;
        pT2scatter = pTscatter * pTscatter;
      }
      // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
      // from rad to deltaCotTheta
      float p2scatter = pT2scatter * iSinTheta2;
      // if deltaTheta larger than allowed scattering for calculated pT, skip
      if ((deltaCotTheta2 - error2 > 0) &&
          (dCotThetaMinusError2 >
            p2scatter * m_finderCfg.sigmaScattering * m_finderCfg.sigmaScattering)) {
        continue;
      }
      // A and B allow calculation of impact params in U/V plane with linear
      // function
      // (in contrast to having to solve a quadratic function in x/y plane)
      float Im = std::abs((A - B * rM) * rM);

      if (Im <= m_finderCfg.impactMax) {
        top_valid.push_back(top[t]);
        // inverse diameter is signed depending if the curvature is
        // positive/negative in phi
        curvatures.push_back(B / std::sqrt(S2));
        impactParameters.push_back(Im);
      }
    }
    if (!top_valid.empty()) {
      m_finderCfg.seedFilter->filterSeeds_2SpFixed(
          *bottom[b], middle, top_valid, curvatures, impactParameters,
          Zob, std::back_inserter(protoseeds));
    }
  }

  return protoseeds;
}

template<std::size_t NDims, typename tree_t>
ActsExamples::SimSeedContainer ActsExamples::SeedingOrthogonalAlgorithm::strategy1(
  const std::vector<const Acts::InternalSpacePoint<ActsExamples::SimSpacePoint> *> spacePoints,
  const tree_t & tree
) const {
  using range_t = typename tree_t::range_t;

  ActsExamples::SimSeedContainer seeds;
  std::size_t timer_0 = 0, timer_1 = 0, timer_2 = 0, timer_3 = 0, timer_4 = 0;

  std::vector<const SP *> bottom_lh_v, bottom_hl_v, top_lh_v, top_hl_v;

  std::vector<std::pair<float, std::unique_ptr<const Acts::InternalSeed<ActsExamples::SimSpacePoint>>>> protoseeds;

  /***************************************************************************
   * STRATEGY 1
   ***************************************************************************/

  for (const SP * middle : spacePoints) {
    std::chrono::high_resolution_clock::time_point tt0 = std::chrono::high_resolution_clock::now();


    if (!validMiddle(*middle, m_finderCfg)) {
      continue;
    }

    range_t bottom_r = validMidToLowOrthoRange<NDims>(*middle, m_finderCfg);
    range_t top_r = validTupleOrthoRange<NDims>(*middle, m_finderCfg);

    double myCotTheta = std::max(std::abs(middle->z() / middle->radius()), m_finderCfg.cotThetaMax);

    std::complex<double> rotPhiM = std::polar(1.0, static_cast<double>(std::asin(-middle->y() / middle->radius())));

    double deltaRMaxTop = top_r[1].max() - middle->radius();
    double deltaRMaxBottom = middle->radius() - bottom_r[1].min();

    range_t bottom_lh_r = bottom_r;
    bottom_lh_r[2].shrink_min(middle->z() - myCotTheta * deltaRMaxBottom);
    bottom_lh_r[2].shrink_max(middle->z());
    range_t top_lh_r = top_r;
    top_lh_r[2].shrink_min(middle->z());
    top_lh_r[2].shrink_max(middle->z() + myCotTheta * deltaRMaxTop);

    range_t bottom_hl_r = bottom_r;
    bottom_hl_r[2].shrink_min(middle->z());
    bottom_hl_r[2].shrink_max(middle->z() + myCotTheta * deltaRMaxBottom);
    range_t top_hl_r = top_r;
    top_hl_r[2].shrink_max(middle->z());
    top_hl_r[2].shrink_min(middle->z() - myCotTheta * deltaRMaxTop);

    // std::cout << "BotLH " << bottom_lh_r.str() << std::endl;
    // std::cout << "BotHL " << bottom_hl_r.str() << std::endl;
    // std::cout << "TopLH " << top_lh_r.str() << std::endl;
    // std::cout << "TopHL " << top_hl_r.str() << std::endl;

    std::complex<double> mc(middle->x(), middle->y());

    std::chrono::high_resolution_clock::time_point tt1 = std::chrono::high_resolution_clock::now();

    if (!bottom_lh_r.degenerate() && !top_lh_r.degenerate()) {
      bottom_lh_v.clear();
      for (const SP * bottom : tree.range_search(bottom_lh_r)) {
        if (validTuple(bottom, middle, m_finderCfg)) {
          bottom_lh_v.push_back(bottom);
        }
      }
    }

    if (!bottom_hl_r.degenerate() && !top_hl_r.degenerate()) {
      bottom_hl_v.clear();
      for (const SP * bottom : tree.range_search(bottom_hl_r)) {
        if (validTuple(bottom, middle, m_finderCfg)) {
          bottom_hl_v.push_back(bottom);
        }
      }
    }

    std::chrono::high_resolution_clock::time_point tt2 = std::chrono::high_resolution_clock::now();

    if (!bottom_lh_v.empty()) {
      top_lh_v.clear();
      for (const SP * top : tree.range_search(top_lh_r)) {
        if (validTuple(middle, top, m_finderCfg)) {
          top_lh_v.push_back(top);
        }
      }
    }
    if (!bottom_hl_v.empty()) {
      top_hl_v.clear();
      for (const SP * top : tree.range_search(top_hl_r)) {
        if (validTuple(middle, top, m_finderCfg)) {
          top_hl_v.push_back(top);
        }
      }
    }

    std::chrono::high_resolution_clock::time_point tt3 = std::chrono::high_resolution_clock::now();

    protoseeds.clear();

    if (!bottom_lh_v.empty() && !top_lh_v.empty()) {
      std::vector<std::pair<float, std::unique_ptr<const Acts::InternalSeed<ActsExamples::SimSpacePoint>>>> lh = helper(*middle, bottom_lh_v, top_lh_v);
      protoseeds.insert(protoseeds.end(), std::make_move_iterator(lh.begin()), std::make_move_iterator(lh.end()));
    }

    if (!bottom_hl_v.empty() && !top_hl_v.empty()) {
      std::vector<std::pair<float, std::unique_ptr<const Acts::InternalSeed<ActsExamples::SimSpacePoint>>>> hl = helper(*middle, bottom_hl_v, top_hl_v);
      protoseeds.insert(protoseeds.end(), std::make_move_iterator(hl.begin()), std::make_move_iterator(hl.end()));
    }

    m_finderCfg.seedFilter->filterSeeds_1SpFixed(protoseeds, std::back_inserter(seeds));

    std::chrono::high_resolution_clock::time_point tt4 = std::chrono::high_resolution_clock::now();

    timer_0 += std::chrono::duration_cast<std::chrono::microseconds>(tt1 - tt0).count();
    timer_1 += std::chrono::duration_cast<std::chrono::microseconds>(tt2 - tt1).count();
    timer_2 += std::chrono::duration_cast<std::chrono::microseconds>(tt3 - tt2).count();
    timer_3 += std::chrono::duration_cast<std::chrono::microseconds>(tt4 - tt3).count();
  }

  std::cout << "Time taken = " <<
    (timer_0 / 1000.0) << "ms, " <<
    (timer_1 / 1000.0) << "ms, " <<
    (timer_2 / 1000.0) << "ms, " <<
    (timer_3 / 1000.0) << "ms, " <<
    (timer_4 / 1000.0) << "ms" <<
    std::endl;

  return seeds;
}

template<std::size_t NDims, typename tree_t>
ActsExamples::SimSeedContainer ActsExamples::SeedingOrthogonalAlgorithm::strategy2(
  const std::vector<const Acts::InternalSpacePoint<ActsExamples::SimSpacePoint> *> spacePoints,
  const tree_t & tree
) const {
  using range_t = typename tree_t::range_t;

  const double maxScatteringAngle = std::sqrt(m_finderCfg.maxScatteringAngle2);
  const double pTperRadius = std::sqrt(m_finderCfg.pT2perRadius);

  ActsExamples::SimSeedContainer seeds;
  std::size_t timer_0 = 0, timer_1 = 0, timer_2 = 0, timer_3 = 0, timer_4 = 0;

  for (const SP * middle : spacePoints) {
    if (!validMiddle(*middle, m_finderCfg)) {
      continue;
    }

    range_t bottom_r = validMidToLowOrthoRange<NDims>(*middle, m_finderCfg);
    range_t top_r = validTupleOrthoRange<NDims>(*middle, m_finderCfg);

    std::complex<double> rotPhiM = std::polar(1.0, static_cast<double>(std::asin(-middle->y() / middle->radius())));

    double myCotTheta = std::max(std::abs(middle->z() / middle->radius()), m_finderCfg.cotThetaMax);

    bottom_r[2].shrink_min(middle->z() - myCotTheta * m_finderCfg.deltaRMax);
    bottom_r[2].shrink_max(middle->z() + myCotTheta * m_finderCfg.deltaRMax);

    std::complex<double> mc(middle->x(), middle->y());

    std::chrono::high_resolution_clock::time_point tt1 = std::chrono::high_resolution_clock::now();

    std::vector<const SP *> bottoms;

    for (const SP * bottom : tree.range_search(bottom_r)) {
      if (validTupleOrtho<NDims>(bottom, middle, m_finderCfg) && validTuple(bottom, middle, m_finderCfg)) {
        bottoms.push_back(bottom);
      }
    }

    std::chrono::high_resolution_clock::time_point tt2 = std::chrono::high_resolution_clock::now();

    std::vector<std::pair<
          float, std::unique_ptr<const Acts::InternalSeed<ActsExamples::SimSpacePoint>>>> protoseeds;

    for (const SP * bottom : bottoms) {
      range_t top_nr = top_r;

      if (middle->z() > bottom->z()) {
        top_nr[2].shrink_min(middle->z());
      } else {
        top_nr[2].shrink_max(middle->z());
      }

      if constexpr (NDims >= 6) {
        double v1 = (middle->z() - bottom->z()) / (middle->radius() - bottom->radius());
        double v2 = middle->z() / middle->radius();

        top_nr[5].shrink_min(std::min(v1, v2));
        top_nr[5].shrink_max(std::max(v1, v2));
      }

      // std::cout << top_nr.str() << std::endl;

      std::vector<const SP *> tops = tree.range_search(top_nr);

      Acts::LinCircle lb = lcBotMid(*bottom, *middle);

      std::complex<double> bc(bottom->x(), bottom->y());

      std::complex<double> qc = (mc - bc);
      std::complex<double> tf = std::complex(1.0, 0.0) / qc;

      float csc = std::sqrt(1. + lb.cotTheta * lb.cotTheta);

      std::vector<const SP *> top_valid;
      std::vector<float> curvatures;
      std::vector<float> impactParameters;

      // std::cout << "We have " << tops.size() << " tops" << std::endl;

      for (const SP * top : tops) {
        if (!validTuple(middle, top, m_finderCfg)) {
          continue;
        }
        std::complex<double> tc(top->x(), top->y());

        std::complex<double> w = (tc - bc) * tf;
        std::complex<double> m = (w - std::norm(w)) / (w - std::conj(w));
        std::complex<double> M = qc * m + bc;
        double norm = std::norm(bc - M);

        if (4 * norm < m_finderCfg.minHelixDiameter2) {
          continue;
        }

        double radius = std::sqrt(norm);

        std::complex<double> I = M - M * radius / std::abs(M);

        if (std::norm(I) > std::pow(m_finderCfg.impactMax, 2)) {
          continue;
        }

        std::complex<double> deltalb = bc - mc;
        std::complex<double> deltalt = tc - mc;

        double ltCotTheta = (top->z() - middle->z()) / std::abs(deltalt);

        double ltEr = ((middle->varianceZ() + top->varianceZ()) +
          (ltCotTheta * ltCotTheta) * (middle->varianceR() + top->varianceR())) / std::norm(deltalt);




        // add errors of spB-spM and spM-spT pairs and add the correlation term
        // for errors on spM
        double error2 = ltEr + lb.Er +
                        2 * (lb.cotTheta * ltCotTheta * middle->varianceR() + middle->varianceZ()) *
                            lb.iDeltaR / std::abs(tc - mc);

        double deltaCotTheta = std::abs(lb.cotTheta - ltCotTheta);

        // if the error is larger than the difference in theta, no need to
        // compare with scattering
        if (deltaCotTheta * deltaCotTheta > error2) {
          double pT2scatter;

          if (m_finderCfg.pTPerHelixRadius * radius > m_finderCfg.maxPtScattering) {
            pT2scatter = m_finderCfg.highland / m_finderCfg.maxPtScattering;
          } else {
            pT2scatter = pTperRadius / radius;
          }

          double max_scat = std::min(maxScatteringAngle, pT2scatter);

          if (deltaCotTheta - std::sqrt(error2) > m_finderCfg.sigmaScattering * csc * max_scat) {
            continue;
          }
        }

        std::complex<double> dlb = deltalb / std::complex(std::norm(deltalb), 0.0);
        std::complex<double> dlt = deltalt / std::complex(std::norm(deltalt), 0.0);
        std::complex<double> ltc = dlt * rotPhiM;
        std::complex<double> lbc = dlb * rotPhiM;
        std::complex<double> myd = ltc - lbc;
        double secMyd = std::abs(myd) / myd.real();
        double tanMyd = myd.imag() / myd.real();
        double sectanMyd = secMyd * tanMyd;

        top_valid.push_back(top);
        impactParameters.push_back(std::abs(I));
        double curvature = lbc.imag() - sectanMyd * lbc.real();
        curvatures.push_back(curvature);
      }

      m_finderCfg.seedFilter->filterSeeds_2SpFixed(
              *bottom, *middle, top_valid, curvatures, impactParameters,
              lb.Zo, std::back_inserter(protoseeds));
    }

    m_finderCfg.seedFilter->filterSeeds_1SpFixed(protoseeds, std::back_inserter(seeds));

    std::chrono::high_resolution_clock::time_point tt3 = std::chrono::high_resolution_clock::now();

    timer_1 += std::chrono::duration_cast<std::chrono::microseconds>(tt2 - tt1).count();
    timer_2 += std::chrono::duration_cast<std::chrono::microseconds>(tt3 - tt2).count();
  }

  std::cout << "Time taken = " <<
    (timer_0 / 1000.0) << "ms, " <<
    (timer_1 / 1000.0) << "ms, " <<
    (timer_2 / 1000.0) << "ms, " <<
    (timer_3 / 1000.0) << "ms, " <<
    (timer_4 / 1000.0) << "ms" <<
    std::endl;

  return seeds;
}

template<std::size_t NDims, typename tree_t>
ActsExamples::SimSeedContainer ActsExamples::SeedingOrthogonalAlgorithm::strategy3(
  const std::vector<const Acts::InternalSpacePoint<ActsExamples::SimSpacePoint> *> spacePoints,
  const tree_t & tree
) const {
  using range_t = typename tree_t::range_t;

  ActsExamples::SimSeedContainer seeds;
  std::size_t timer_0 = 0, timer_1 = 0, timer_2 = 0, timer_3 = 0, timer_4 = 0;

  for (const SP * middle : spacePoints) {
    if (!validMiddle(*middle, m_finderCfg)) {
      continue;
    }

    range_t bottom_r = validMidToLowOrthoRange<NDims>(*middle, m_finderCfg);
    range_t top_r = validTupleOrthoRange<NDims>(*middle, m_finderCfg);

    double myCotTheta = std::max(std::abs(middle->z() / middle->radius()), m_finderCfg.cotThetaMax);


    bottom_r[2].shrink_min(middle->z() - myCotTheta * m_finderCfg.deltaRMax);
    bottom_r[2].shrink_max(middle->z() + myCotTheta * m_finderCfg.deltaRMax);

    top_r[2].shrink_min(middle->z() - myCotTheta * m_finderCfg.deltaRMax);
    top_r[2].shrink_max(middle->z() + myCotTheta * m_finderCfg.deltaRMax);



    std::chrono::high_resolution_clock::time_point tt1 = std::chrono::high_resolution_clock::now();

    std::vector<const SP *> tops;

    for (const SP * top : tree.range_search(top_r)) {
      if (validTupleOrtho<NDims>(middle, top, m_finderCfg) && validTuple(middle, top, m_finderCfg)) {
        tops.push_back(top);
      }
    }

    // std::cout << "Top size " << tops.size() << std::endl;

    std::chrono::high_resolution_clock::time_point tt2 = std::chrono::high_resolution_clock::now();

    std::vector<std::pair<
          float, std::unique_ptr<const Acts::InternalSeed<ActsExamples::SimSpacePoint>>>> protoseeds;

    std::map<const SP *, std::set<const SP *>> bottom_top_map;

    for (const SP * top : tops) {
      range_t bottom_nr = bottom_r;

      if (middle->z() > top->z()) {
        bottom_nr[2].shrink_min(middle->z());
      } else {
        bottom_nr[2].shrink_max(middle->z());
      }

      double slope = (top->radius() - middle->radius()) / (top->z() - middle->z());
      double offset = top->radius() - slope * top->z();
      double zc = (bottom_nr[1].min() - offset) / slope;

      // std::cout << slope << ", " << offset << ", " << zc << std::endl;

      if (zc > middle->z()) {
        bottom_nr[2].shrink_max(zc);
      } else {
        bottom_nr[2].shrink_min(zc);
      }

      if constexpr (NDims >= 6) {
        // double cotTheta = middle->z() / middle->radius();
        // bottom_nr[5].shrink_max(std::max(zc / bottom_nr[1].min(), cotTheta));
        // bottom_nr[5].shrink_min(std::min(zc / bottom_nr[1].min(), cotTheta));
      }

      std::vector<const SP *> bottoms = tree.range_search(bottom_nr);

      // std::cout << bottom_nr.str() << std::endl;
      // std::cout << "Bottom size " << bottoms.size() << std::endl;

      for (const SP * bottom : bottoms) {
        if (!validTupleOrtho<NDims>(bottom, middle, m_finderCfg) || !validTuple(bottom, middle, m_finderCfg)) {
          continue;
        }

        std::optional<std::pair<double, double>> res = validTriple(*bottom, *middle, *top);

        if (!res) {
          continue;
        }

        bottom_top_map.try_emplace(bottom);
        bottom_top_map.at(bottom).insert(top);
      }
    }

    std::chrono::high_resolution_clock::time_point tt3 = std::chrono::high_resolution_clock::now();

    for (auto & p : bottom_top_map) {
      std::vector<float> curvatures;
      std::vector<float> impactParameters;
      std::vector<const SP *> mytops;
      Acts::LinCircle lb = lcBotMid(*p.first, *middle);

      for (const SP * t : p.second) {
        std::optional<std::pair<double, double>> res = validTriple(*p.first, *middle, *t);

        impactParameters.push_back(res->first);
        curvatures.push_back(res->second);
        mytops.push_back(t);
      }

      m_finderCfg.seedFilter->filterSeeds_2SpFixed(
              *p.first, *middle, mytops, curvatures, impactParameters,
              lb.Zo, std::back_inserter(protoseeds));
    }

    m_finderCfg.seedFilter->filterSeeds_1SpFixed(protoseeds, std::back_inserter(seeds));

    std::chrono::high_resolution_clock::time_point tt4 = std::chrono::high_resolution_clock::now();

    timer_0 += std::chrono::duration_cast<std::chrono::microseconds>(tt2 - tt1).count();
    timer_1 += std::chrono::duration_cast<std::chrono::microseconds>(tt3 - tt2).count();
    timer_2 += std::chrono::duration_cast<std::chrono::microseconds>(tt4 - tt3).count();
  }

  std::cout << "Time taken = " <<
    (timer_0 / 1000.0) << "ms, " <<
    (timer_1 / 1000.0) << "ms, " <<
    (timer_2 / 1000.0) << "ms, " <<
    (timer_3 / 1000.0) << "ms, " <<
    (timer_4 / 1000.0) << "ms" <<
    std::endl;

  return seeds;
}


ActsExamples::ProcessCode ActsExamples::SeedingOrthogonalAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
  // construct the combined input container of space point pointers from all
  // configured input sources.
  // pre-compute the total size required so we only need to allocate once
  size_t nSpacePoints = 0;
  for (const auto& isp : m_cfg.inputSpacePoints) {
    nSpacePoints += ctx.eventStore.get<SimSpacePointContainer>(isp).size();
  }
  std::vector<const Acts::InternalSpacePoint<SimSpacePoint> *> spacePoints;
  spacePoints.reserve(nSpacePoints);
  for (const auto& isp : m_cfg.inputSpacePoints) {
    for (const auto& spacePoint :
         ctx.eventStore.get<SimSpacePointContainer>(isp)) {
      // since the event store owns the space points, their pointers should be
      // stable and we do not need to create local copies.
      spacePoints.push_back(new Acts::InternalSpacePoint<SimSpacePoint>(
        spacePoint,
        {spacePoint.x(), spacePoint.y(), spacePoint.z()},
        {0.0, 0.0},
        {spacePoint.varianceR(), spacePoint.varianceZ()}
      ));
    }
  }

  constexpr std::size_t NDims = 3;

  std::vector<std::pair<std::array<double, NDims>, const SP *>> points;

  enum class Dim {
    Phi = 0,
    Radius = 1,
    Z = 2,
    Y = 3,
    X = 4,
    CotTheta = 5
  };

  for (auto sp : spacePoints) {
    std::array<double, NDims> point;

    if constexpr (NDims >= 6) {
      point[5] = sp->z() / sp->radius();
    }

    if constexpr (NDims >= 5) {
      point[3] = sp->y();
      point[4] = sp->x();
    }

    if constexpr (NDims >= 3) {
      point[0] = sp->phi();
      point[1] = sp->radius();
      point[2] = sp->z();
    }

    points.emplace_back(point, sp);

    // if (point[0] <= -3.14159265 + 0.1) {
    //   std::array<double, NDims> point_dup = point;
    //   point_dup[0] += 2 * 3.14159265;
    //   points.emplace_back(point_dup, sp);
    // }

    // if (point[0] >= 3.14159265 - 0.1) {
    //   std::array<double, NDims> point_dup = point;
    //   point_dup[0] -= 2 * 3.14159265;
    //   points.emplace_back(point_dup, sp);
    // }
  }

  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  Acts::KDTree<NDims, const SP *, Acts::ActsScalar, std::array, 4> tree(std::move(points));

  // tree.print();

  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

  // run the seeding
  SimSeedContainer seeds = strategy1<NDims>(spacePoints, tree);


  std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();

  // extract proto tracks, i.e. groups of measurement indices, from tracks seeds
  size_t nSeeds = seeds.size();
  ProtoTrackContainer protoTracks;
  protoTracks.reserve(nSeeds);
  for (const auto& seed : seeds) {
    ProtoTrack protoTrack;
    protoTrack.reserve(seed.sp().size());
    for (auto spacePointPtr : seed.sp()) {
      protoTrack.push_back(spacePointPtr->measurementIndex());
    }
    protoTracks.push_back(std::move(protoTrack));
  }

  std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();

  std::cout << "Preprocessing:     " << (std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0) << "ms" << std::endl;
  std::cout << "k-d tree creation: " << (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.0) << "ms" << std::endl;
  std::cout << "Seed finding:      " << (std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count() / 1000.0) << "ms" << std::endl;
  std::cout << "Postprocessing:    " << (std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count() / 1000.0) << "ms" << std::endl;

  ACTS_DEBUG("Created " << seeds.size() << " track seeds from "
                        << spacePoints.size() << " space points");

  ctx.eventStore.add(m_cfg.outputSeeds, std::move(seeds));
  ctx.eventStore.add(m_cfg.outputProtoTracks, std::move(protoTracks));
  return ActsExamples::ProcessCode::SUCCESS;
}
