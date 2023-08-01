// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "Acts/Vertexing/VertexingError.hpp"

#include <algorithm>

template <int trkGridSize>
typename Acts::AdaptiveGridTrackDensity<trkGridSize>::DensityMap::const_iterator Acts::AdaptiveGridTrackDensity<trkGridSize>::highestDensityEntry(const DensityMap& densityMap) const {
  auto maxEntry = std::max_element
  (
      std::begin(densityMap), std::end(densityMap),
      [] (const auto& densityEntry1, const auto& densityEntry2) {
          return densityEntry1.second < densityEntry2.second;
      }
  );
  return maxEntry;
}
/// TODO use .at() function to retrieve/modify map elems

template <int trkGridSize>
Acts::Result<float>
Acts::AdaptiveGridTrackDensity<trkGridSize>::getMaxZPosition(
    std::vector<float>& mainGridDensity,
    const std::vector<int>& mainGridZValues) const {

  DensityMap mainDensityMap;
  for (unsigned int i = 0; i < mainGridDensity.size(); i++) {
    mainDensityMap[mainGridZValues[i]] = mainGridDensity[i];
  }

  if (mainDensityMap.empty()) {
    return VertexingError::EmptyInput;
  }

  int zBin = -1;
  if (!m_cfg.useHighestSumZPosition) {
    zBin = highestDensityEntry(mainDensityMap)->first;
    //std::cout << "\n mainGridZValues[zGridPos], maxEntry:\n" << mainGridZValues[zGridPos] << " " << maxEntry->first << "\n";
  } else {
    // Get z position with highest density sum
    // of surrounding bins
    zBin = highestDensitySumBin(mainDensityMap);
  }

  // Derive corresponding z value
  int sign = (zBin > 0) ? +1 : -1;
  return (zBin + sign * 0.5f) * m_cfg.binSize;
}

template <int trkGridSize>
Acts::Result<float>
Acts::AdaptiveGridTrackDensity<trkGridSize>::getMaxZPosition(DensityMap& densityMap) const {

  if (densityMap.empty()) {
    return VertexingError::EmptyInput;
  }

  int zBin = -1;
  if (!m_cfg.useHighestSumZPosition) {
    zBin = highestDensityEntry(densityMap)->first;
    //std::cout << "\n mainGridZValues[zGridPos], maxEntry:\n" << mainGridZValues[zGridPos] << " " << maxEntry->first << "\n";
  } else {
    // Get z position with highest density sum
    // of surrounding bins
    zBin = highestDensitySumBin(densityMap);
  }

  // Derive corresponding z value
  int sign = (zBin > 0) ? +1 : -1;
  return (zBin + sign * 0.5f) * m_cfg.binSize;
}

template <int trkGridSize>
Acts::Result<std::pair<float, float>>
Acts::AdaptiveGridTrackDensity<trkGridSize>::getMaxZPositionAndWidth(
    std::vector<float>& mainGridDensity,
    const std::vector<int>& mainGridZValues) const {
  // Get z maximum value
  auto maxZRes = getMaxZPosition(mainGridDensity, mainGridZValues);
  if (not maxZRes.ok()) {
    return maxZRes.error();
  }
  float maxZ = *maxZRes;

  // Get seed width estimate
  auto widthRes = estimateSeedWidth(mainGridDensity, mainGridZValues, maxZ);
  if (not widthRes.ok()) {
    return widthRes.error();
  }
  float width = *widthRes;
  std::pair<float, float> returnPair{maxZ, width};
  return returnPair;
}

template <int trkGridSize>
Acts::Result<std::pair<float, float>>
Acts::AdaptiveGridTrackDensity<trkGridSize>::getMaxZPositionAndWidth(DensityMap& densityMap) const {
  // Get z maximum value
  auto maxZRes = getMaxZPosition(densityMap);
  if (not maxZRes.ok()) {
    return maxZRes.error();
  }
  float maxZ = *maxZRes;

  // Get seed width estimate
  auto widthRes = estimateSeedWidth(densityMap, maxZ);
  if (not widthRes.ok()) {
    return widthRes.error();
  }
  float width = *widthRes;
  std::pair<float, float> returnPair{maxZ, width};
  return returnPair;
}

template <int trkGridSize>
std::pair<int,
          typename Acts::AdaptiveGridTrackDensity<trkGridSize>::TrackGridVector>
Acts::AdaptiveGridTrackDensity<trkGridSize>::addTrack(
    const Acts::BoundTrackParameters& trk, std::vector<float>& mainGridDensity,
    std::vector<int>& mainGridZValues) const {
  SymMatrix2 cov = trk.covariance().value().block<2, 2>(0, 0);
  float d0 = trk.parameters()[0];
  float z0 = trk.parameters()[1];

  DensityMap mainDensityMap;
  for (unsigned int i = 0; i < mainGridDensity.size(); i++) {
    mainDensityMap[mainGridZValues[i]] = mainGridDensity[i];
  }

  // Calculate offset in d direction to central bin at z-axis
  int dOffset = static_cast<int>(std::floor(d0 / m_cfg.binSize - 0.5) + 1);
  // Check if current track affects grid density
  // in central bins at z-axis
  if (std::abs(dOffset) > (trkGridSize - 1) / 2.) {
    return {0, TrackGridVector::Zero()};
  }
  // Calculate bin in z
  int zBin = int(z0 / m_cfg.binSize);

  // Calculate the positions of the bin centers
  float binCtrD = dOffset * m_cfg.binSize;

  int sign = (z0 > 0) ? +1 : -1;
  float binCtrZ = (zBin + sign * 0.5f) * m_cfg.binSize;

  // Calculate the distance between IP values and their
  // corresponding bin centers
  float distCtrD = d0 - binCtrD;
  float distCtrZ = z0 - binCtrZ;

  DensityMap trackDensityMap = createTrackGrid(dOffset, cov, distCtrD, zBin, distCtrZ);

  TrackGridVector trackGrid;
  int cnt4 = 0;
  for (const auto& densityEntry : trackDensityMap) {
    trackGrid[cnt4] = densityEntry.second;
    cnt4 ++;
  }

  for (const auto& densityEntry : trackDensityMap) {
    int zBin2 = densityEntry.first;
    float trackDensity = densityEntry.second;
    // Check if z bin is already part of the main grid
    if (mainDensityMap.count(zBin2) == 1) {
      mainDensityMap[zBin2] += trackDensity;
    }
    else {
      mainDensityMap[zBin2] = trackDensity;
    }
  }

  mainGridZValues.clear();
  mainGridDensity.clear();
  for (const auto& densityEntry : mainDensityMap) {
    mainGridZValues.push_back(densityEntry.first);
    mainGridDensity.push_back(densityEntry.second);
  }

  /*
  std::cout << "\n new track:\n";
  std::cout << "\n central z bin: " << zBin << "\n";
  int cnt2 = 0;
  for (auto i = mainDensityMap.begin(); i!= mainDensityMap.end(); i++) {
    std::cout << "\n zind: " << i->first << " density: " << i->second << " trackGridDensity: " << mainGridDensity[cnt2] << "\n";
    cnt2 ++;
  }
*/
  return {zBin, trackGrid};
}

template <int trkGridSize>
typename Acts::AdaptiveGridTrackDensity<trkGridSize>::DensityMap
Acts::AdaptiveGridTrackDensity<trkGridSize>::addTrack(
    const Acts::BoundTrackParameters& trk, DensityMap& mainDensityMap) const {
  SymMatrix2 cov = trk.covariance().value().block<2, 2>(0, 0);
  float d0 = trk.parameters()[0];
  float z0 = trk.parameters()[1];

  // Calculate offset in d direction to central bin at z-axis
  int dOffset = static_cast<int>(std::floor(d0 / m_cfg.binSize - 0.5) + 1);
  // Check if current track affects grid density
  // in central bins at z-axis
  if (std::abs(dOffset) > (trkGridSize - 1) / 2.) {
    DensityMap emptyTrackDensityMap;
    return emptyTrackDensityMap;
  }
  // Calculate bin in z
  int zBin = int(z0 / m_cfg.binSize);

  // Calculate the positions of the bin centers
  float binCtrD = dOffset * m_cfg.binSize;

  int sign = (z0 > 0) ? +1 : -1;
  float binCtrZ = (zBin + sign * 0.5f) * m_cfg.binSize;

  // Calculate the distance between IP values and their
  // corresponding bin centers
  float distCtrD = d0 - binCtrD;
  float distCtrZ = z0 - binCtrZ;

  DensityMap trackDensityMap = createTrackGrid(dOffset, cov, distCtrD, zBin, distCtrZ);

  for (const auto& densityEntry : trackDensityMap) {
    int zBin2 = densityEntry.first;
    float trackDensity = densityEntry.second;
    // Check if z bin is already part of the main grid
    if (mainDensityMap.count(zBin2) == 1) {
      mainDensityMap[zBin2] += trackDensity;
    }
    else {
      mainDensityMap[zBin2] = trackDensity;
    }
  }

  return trackDensityMap;
}

template <int trkGridSize>
void Acts::AdaptiveGridTrackDensity<trkGridSize>::removeTrackGridFromMainGrid(
    int zBin, const TrackGridVector& trkGrid,
    std::vector<float>& mainGridDensity,
    const std::vector<int>& mainGridZValues) const {
  // Find position of current z bin in mainGridZValues
  auto findIter =
      std::find(mainGridZValues.begin(), mainGridZValues.end(), zBin);
  // Calculate corresponding index in mainGridDensity
  int densityIdx = std::distance(mainGridZValues.begin(), findIter);

  // Go over trkGrid and remove it from mainDensityGrid
  int startEnd = int((trkGridSize - 1) / 2);
  for (int i = 0; i < trkGridSize; i++) {
    mainGridDensity[int(densityIdx + (i - startEnd))] -= trkGrid[i];
  }
}

template <int trkGridSize>
void Acts::AdaptiveGridTrackDensity<trkGridSize>::subtractTrack(const DensityMap& trackDensityMap, DensityMap& mainDensityMap) const {
  for (auto it = trackDensityMap.begin(); it != trackDensityMap.end(); it ++) {
    mainDensityMap.at(it->first) -= it->second;
  }
}

template <int trkGridSize>
typename Acts::AdaptiveGridTrackDensity<trkGridSize>::DensityMap
Acts::AdaptiveGridTrackDensity<trkGridSize>::createTrackGrid(
    int offset, const Acts::SymMatrix2& cov, float distCtrD,
    int centralZBin, float distCtrZ) const {
  DensityMap trackDensityMap;

  int halfTrkGridSize = (trkGridSize - 1) / 2;
  // TODO: I think we can just pass d0
  float i = halfTrkGridSize + offset;
  float d = (i - static_cast<float>(trkGridSize) / 2 + 0.5f) * m_cfg.binSize;

  int firstZBin = centralZBin - halfTrkGridSize;
  // Loop over columns
  for (int j = 0; j < trkGridSize; j++) {
    float z = (j - static_cast<float>(trkGridSize) / 2 + 0.5f) * m_cfg.binSize;
    trackDensityMap[firstZBin + j] = normal2D(d + distCtrD, z + distCtrZ, cov);
  }
  return trackDensityMap;
}

template <int trkGridSize>
Acts::Result<float>
Acts::AdaptiveGridTrackDensity<trkGridSize>::estimateSeedWidth(
    const std::vector<float>& mainGridDensity,
    const std::vector<int>& mainGridZValues, float maxZ) const {
  
  DensityMap mainDensityMap;
  for (unsigned int i = 0; i < mainGridDensity.size(); i++) {
    mainDensityMap[mainGridZValues[i]] = mainGridDensity[i];
  }

  if (mainDensityMap.empty()) {
    return VertexingError::EmptyInput;
  }

  // Get z bin of max density and max density value
  int sign = (maxZ > 0) ? +1 : -1;
  int zMaxBin = int(maxZ / m_cfg.binSize - sign * 0.5f);
  const float maxValue = mainDensityMap[zMaxBin];

  int rhmBin = zMaxBin;
  float gridValue = maxValue;
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continuous z values
    if (mainDensityMap.count(rhmBin+1) == 0) {
      // TODO: I think in this case the linear approximation does not work, we would need to use 0 instead of mainDensityMap.at(rhmBin) in the computation of deltaZ1
      break;
    }
    rhmBin += 1;
    gridValue = mainDensityMap.at(rhmBin);
  }

  // Use linear approximation to find better z value for FWHM between bins
  float deltaZ1 =
      (maxValue / 2 - mainDensityMap.at(rhmBin - 1)) *
      (m_cfg.binSize / (mainDensityMap.at(rhmBin - 1) - mainDensityMap.at(rhmBin)));
  // TODO: I think we should have here:
  // float deltaZ1 =
  //       (maxValue / 2 - mainDensityMap.at(rhmBin)) *
  //       (m_cfg.binSize / (mainDensityMap.at(rhmBin - 1) - mainDensityMap.at(rhmBin)));
  int lhmBin = zMaxBin;
  gridValue = maxValue;
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continuous z values
    if (mainDensityMap.count(lhmBin-1) == 0) {
      // TODO: I think in this case the linear approximation does not work
      break;
    }
    lhmBin -= 1;
    gridValue = mainDensityMap.at(lhmBin);
  }

  // Use linear approximation to find better z value for FWHM between bins
  float deltaZ2 =
      (maxValue / 2 - mainDensityMap.at(lhmBin + 1)) *
      (m_cfg.binSize / (mainDensityMap.at(rhmBin + 1) - mainDensityMap.at(rhmBin)));
  // TODO: I think we should have here:
  // float deltaZ2 =
  //    (maxValue / 2 - mainDensityMap.at(lhmBin)) *
  //    (m_cfg.binSize / (mainDensityMap.at(lhmBin + 1) - mainDensityMap.at(lhmBin)));

  float fwhm = (rhmBin - lhmBin) * m_cfg.binSize - deltaZ1 - deltaZ2;

  // FWHM = 2.355 * sigma
  float width = fwhm / 2.355f;

  return std::isnormal(width) ? width : 0.0f;
}

template <int trkGridSize>
Acts::Result<float>
Acts::AdaptiveGridTrackDensity<trkGridSize>::estimateSeedWidth(const DensityMap& densityMap, float maxZ) const {

  if (densityMap.empty()) {
    return VertexingError::EmptyInput;
  }

  // Get z bin of max density and max density value
  int sign = (maxZ > 0) ? +1 : -1;
  int zMaxBin = int(maxZ / m_cfg.binSize - sign * 0.5f);
  const float maxValue = densityMap.at(zMaxBin);

  int rhmBin = zMaxBin;
  float gridValue = maxValue;
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continuous z values
    if (densityMap.count(rhmBin+1) == 0) {
      // TODO: I think in this case the linear approximation does not work, we would need to use 0 instead of densityMap.at(rhmBin) in the computation of deltaZ1
      break;
    }
    rhmBin += 1;
    gridValue = densityMap.at(rhmBin);
  }

  // Use linear approximation to find better z value for FWHM between bins
  float deltaZ1 =
      (maxValue / 2 - densityMap.at(rhmBin - 1)) *
      (m_cfg.binSize / (densityMap.at(rhmBin - 1) - densityMap.at(rhmBin)));
  // TODO: I think we should have here:
  // float deltaZ1 =
  //       (maxValue / 2 - densityMap.at(rhmBin)) *
  //       (m_cfg.binSize / (densityMap.at(rhmBin - 1) - densityMap.at(rhmBin)));
  int lhmBin = zMaxBin;
  gridValue = maxValue;
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continuous z values
    if (densityMap.count(lhmBin-1) == 0) {
      // TODO: I think in this case the linear approximation does not work
      break;
    }
    lhmBin -= 1;
    gridValue = densityMap.at(lhmBin);
  }

  // Use linear approximation to find better z value for FWHM between bins
  float deltaZ2 =
      (maxValue / 2 - densityMap.at(lhmBin + 1)) *
      (m_cfg.binSize / (densityMap.at(rhmBin + 1) - densityMap.at(rhmBin)));
  // TODO: I think we should have here:
  // float deltaZ2 =
  //    (maxValue / 2 - densityMap.at(lhmBin)) *
  //    (m_cfg.binSize / (densityMap.at(lhmBin + 1) - densityMap.at(lhmBin)));

  float fwhm = (rhmBin - lhmBin) * m_cfg.binSize - deltaZ1 - deltaZ2;

  // FWHM = 2.355 * sigma
  float width = fwhm / 2.355f;

  return std::isnormal(width) ? width : 0.0f;
}

template <int trkGridSize>
float Acts::AdaptiveGridTrackDensity<trkGridSize>::normal2D(
    float d, float z, const Acts::SymMatrix2& cov) const {
  float det = cov.determinant();
  float coef = 1 / (2 * M_PI * std::sqrt(det));
  float expo =
      -1 / (2 * det) *
      (cov(1, 1) * d * d - d * z * (cov(0, 1) + cov(1, 0)) + cov(0, 0) * z * z);
  return coef * std::exp(expo);
}

template <int trkGridSize>
int Acts::AdaptiveGridTrackDensity<trkGridSize>::highestDensitySumBin(
    DensityMap& densityMap) const {
  // Checks the first (up to) 3 density maxima, if they are close, checks which
  // one has the highest surrounding density sum (the two neighboring bins)

  // The global maximum
  auto firstMax = highestDensityEntry(densityMap);
  int zFirstMax = firstMax->first;
  float valueFirstMax = firstMax->second;
  float firstSum = getDensitySum(densityMap, zFirstMax);
  float densityDeviation = valueFirstMax * m_cfg.maxRelativeDensityDev;

  // Get the second highest maximum
  densityMap[zFirstMax] = 0;
  auto secondMax = highestDensityEntry(densityMap);
  int zSecondMax = secondMax->first;
  float valueSecondMax = secondMax->second;
  float secondSum = 0;
  if (valueFirstMax - valueSecondMax < densityDeviation) {
    secondSum = getDensitySum(densityMap, zSecondMax);
  }

  // Get the third highest maximum
  densityMap[zSecondMax] = 0;
  auto thirdMax = highestDensityEntry(densityMap);
  int zThirdMax = thirdMax->first;
  float valueThirdMax = thirdMax->second;
  float thirdSum = 0;
  if (valueFirstMax - valueThirdMax < densityDeviation) {
    thirdSum = getDensitySum(densityMap, zThirdMax);
  }

  // Revert back to original values
  densityMap[zFirstMax] = valueFirstMax;
  densityMap[zSecondMax] = valueSecondMax;

  // Return the z-bin position of the highest density sum
  if (secondSum > firstSum && secondSum > thirdSum) {
    return zSecondMax;
  }
  if (thirdSum > secondSum && thirdSum > firstSum) {
    return zThirdMax;
  }
  return zFirstMax;
}

template <int trkGridSize>
float Acts::AdaptiveGridTrackDensity<trkGridSize>::getDensitySum(
     const DensityMap& densityMap, int zBin) const {
  float sum = densityMap.at(zBin);
  // Check if neighboring bins are part of the densityMap and add them (if they are not part of the map, we assume them to be 0)
  // Note that each key in a map is unique; the .count() function therefore returns either 0 or 1
  if (densityMap.count(zBin - 1) == 1) {
    sum += densityMap.at(zBin - 1);
  }
  if (densityMap.count(zBin + 1) == 1) {
    sum += densityMap.at(zBin + 1);
  }
  return sum;
}
