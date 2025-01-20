// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/AccumulatedMaterialSlab.hpp"

#include "Acts/Material/Material.hpp"
#include "Acts/Material/detail/AverageMaterials.hpp"

void Acts::AccumulatedMaterialSlab::accumulate(MaterialSlab slab,
                                               double pathCorrection) {
  // scale the recorded material to the equivalence contribution along the
  // surface normal
  slab.scaleThickness(static_cast<float>(1. / pathCorrection));
  m_trackAverage = detail::combineSlabs(m_trackAverage, slab);
}

void Acts::AccumulatedMaterialSlab::trackVariance(MaterialSlab slabReference,
                                                  bool useEmptyTrack) {
  // Only use real tracks or if empty tracks are allowed.
  if (useEmptyTrack || (0 < m_trackAverage.thickness())) {
    float variance = ((1 / m_trackAverage.material().X0()) -
                      (1 / slabReference.material().X0())) *
                     ((1 / m_trackAverage.material().X0()) -
                      (1 / slabReference.material().X0()));
    if (m_totalCount == 0u) {
      m_totalVariance = variance;
    } else {
      double weightTotal = m_totalCount / (m_totalCount + 1.0);
      double weightTrack = 1 / (m_totalCount + 1.0);
      m_totalVariance = weightTotal * m_totalVariance + weightTrack * variance;
    }
  }
}

void Acts::AccumulatedMaterialSlab::trackAverage(bool useEmptyTrack) {
  // average only real tracks or if empty tracks are allowed.
  if (useEmptyTrack || (0 < m_trackAverage.thickness())) {
    if (m_totalCount == 0u) {
      m_totalAverage = m_trackAverage;
    } else {
      double weightTotal = m_totalCount / (m_totalCount + 1.0);
      double weightTrack = 1 / (m_totalCount + 1.0);
      // average such that each track contributes equally.
      MaterialSlab fromTotal(m_totalAverage.material(),
                             weightTotal * m_totalAverage.thickness());
      MaterialSlab fromTrack(m_trackAverage.material(),
                             weightTrack * m_trackAverage.thickness());
      m_totalAverage = detail::combineSlabs(fromTotal, fromTrack);
    }
    m_totalCount += 1;
  }
  // reset track average
  m_trackAverage = MaterialSlab();
}

std::pair<Acts::MaterialSlab, unsigned int>
Acts::AccumulatedMaterialSlab::totalAverage() const {
  return {m_totalAverage, m_totalCount};
}

std::pair<float, unsigned int> Acts::AccumulatedMaterialSlab::totalVariance()
    const {
  return {m_totalVariance, m_totalCount};
}
