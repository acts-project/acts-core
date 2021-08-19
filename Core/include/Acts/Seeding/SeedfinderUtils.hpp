#pragma once

#include "Acts/Utilities/KDRange.hpp"

template<typename SP>
Acts::LinCircle _lincircle(
  const SP & middle,
  const SP & other,
  bool bottom
) {
  float xM = middle.x();
  float yM = middle.y();
  float zM = middle.z();
  float rM = middle.radius();
  float varianceZM = middle.varianceZ();
  float varianceRM = middle.varianceR();
  float cosPhiM = xM / rM;
  float sinPhiM = yM / rM;

  float deltaX = other.x() - xM;
  float deltaY = other.y() - yM;
  float deltaZ = other.z() - zM;

  float x = deltaX * cosPhiM + deltaY * sinPhiM;
  float y = deltaY * cosPhiM - deltaX * sinPhiM;

  float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
  float iDeltaR = std::sqrt(iDeltaR2);

  int bottomFactor = 1 * (int(!bottom)) - 1 * (int(bottom));

  float cot_theta = deltaZ * iDeltaR * bottomFactor;

  Acts::LinCircle l;
  l.cotTheta = cot_theta;

  l.Zo = zM - rM * cot_theta;
  l.iDeltaR = iDeltaR;

  l.U = x * iDeltaR2;
  l.V = y * iDeltaR2;

  l.Er = ((varianceZM + other.varianceZ()) +
          (cot_theta * cot_theta) * (varianceRM + other.varianceR())) *
          iDeltaR2;
  return l;
}

template<typename SP>
Acts::LinCircle lcBotMid(
  const SP & bottom,
  const SP & middle
) {
  return _lincircle(middle, bottom, true);
}

template<typename SP>
Acts::LinCircle lcMidTop(
  const SP & middle,
  const SP & top
) {
  return _lincircle(middle, top, false);
}

template<typename SP>
bool validMiddle(
  const SP & s,
  const Acts::SeedfinderConfig<ActsExamples::SimSpacePoint> & config
) {
  if (s.radius() > 120) {
    return false;
  }

  if (s.radius() < 60) {
    return false;
  }

  if (std::abs(s.z() / s.radius()) > config.cotThetaMax) {
    return false;
  }

  return true;
}

template<std::size_t N, typename SP>
Acts::KDRange<N, double> validTupleRangeBase(
  const Acts::SeedfinderConfig<ActsExamples::SimSpacePoint> & config
) {
  Acts::KDRange<N, double> res;

  res[0].shrink_min(config.phiMin);
  res[0].shrink_max(config.phiMax);

  // res[1].shrink_min(config.rMin);
  res[1].shrink_max(config.rMax);

  res[2].shrink_min(config.zMin);
  res[2].shrink_max(config.zMax);

  // if constexpr(N >= 6) {
  //   res[5].shrink_min(-config.cotThetaMax);
  //   res[5].shrink_max(config.cotThetaMax);
  // }

  return res;
}

template<std::size_t N, typename SP>
Acts::KDRange<N, double> validTupleOrthoRange(
  const SP & low,
  const Acts::SeedfinderConfig<ActsExamples::SimSpacePoint> & config
) {
  double pL = low.phi();
  double rL = low.radius();
  double zL = low.z();

  Acts::KDRange<N, double> res = validTupleRangeBase<N, SP>(config);

  res[1].shrink_min(rL + config.deltaRMin);
  res[1].shrink_max(rL + config.deltaRMax);

  double zMax = (res[1].max() / rL) * (zL - config.collisionRegionMin) + config.collisionRegionMin;
  double zMin = config.collisionRegionMax - (res[1].max() / rL) * (config.collisionRegionMax - zL);

  if (zL > config.collisionRegionMin) {
    res[2].shrink_max(zMax);
  }

  if (zL < config.collisionRegionMax) {
    res[2].shrink_min(zMin);
  }

  res[2].shrink_min(zL - config.cotThetaMax * (res[1].max() - rL));
  res[2].shrink_max(zL + config.cotThetaMax * (res[1].max() - rL));

  double delta_phi = 0.065;

  res[0].shrink_min(pL - delta_phi);
  res[0].shrink_max(pL + delta_phi);

  return res;
}

template<std::size_t N, typename SP>
Acts::KDRange<N, double> validMidToLowOrthoRange(
  const SP & mid,
  const Acts::SeedfinderConfig<ActsExamples::SimSpacePoint> & config
) {
  float pM = mid.phi();
  float rM = mid.radius();

  Acts::KDRange<N, double> res = validTupleRangeBase<N, SP>(config);

  res[1].shrink_min(rM - config.deltaRMax);
  res[1].shrink_max(rM - config.deltaRMin);

  double frac_r = res[1].min() / rM;

  float zMin = (mid.z() - config.collisionRegionMin) * frac_r + config.collisionRegionMin;
  float zMax = (mid.z() - config.collisionRegionMax) * frac_r + config.collisionRegionMax;

  res[2].shrink_min(std::min(zMin, mid.z()));
  res[2].shrink_max(std::max(zMax, mid.z()));

  double delta_phi = 0.065;

  res[0].shrink_min(pM - delta_phi);
  res[0].shrink_max(pM + delta_phi);

  return res;
}

template<std::size_t N, typename SP>
bool validTupleOrtho(
  const SP * low,
  const SP * high,
  const Acts::SeedfinderConfig<ActsExamples::SimSpacePoint> & config
) {
  Acts::KDRange<N, double> r = validTupleOrthoRange<N>(*low, config);
  typename Acts::KDRange<N, double>::coordinate_t p = {high->phi(), high->radius(), high->z()};

  return r.contains(p);
}

template<typename SP>
bool validTuple(
  const SP * low,
  const SP * high,
  const Acts::SeedfinderConfig<ActsExamples::SimSpacePoint> & config
) {
  float rL = low->radius();
  float rH = high->radius();

  float zL = low->z();
  float zH = high->z();

  float deltaR = rH - rL;

  // ratio Z/R (forward angle) of space point duplet
  // float cotTheta2 = (zH - zL) / config.deltaRMax;
  // if (std::fabs(cotTheta2) > config.cotThetaMax) {
  //   std::cout << "Tuple discriminator 0.1" << std::endl;
  //   return false;
  // }

  // ratio Z/R (forward angle) of space point duplet
  float cotTheta = (zH - zL) / deltaR;
  if (std::fabs(cotTheta) > config.cotThetaMax) {
    // std::cout << "Tuple discriminator 1" << std::endl;
    return false;
  }
  // check if duplet origin on z axis within collision region
  float zOrigin = zL - rL * cotTheta;
  if (zOrigin < config.collisionRegionMin ||
      zOrigin > config.collisionRegionMax) {
    // std::cout << "Tuple discriminator 2" << std::endl;
    return false;
  }

  return true;
}
