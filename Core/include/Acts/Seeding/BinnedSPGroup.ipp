// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename external_spacepoint_t>
template <typename spacepoint_iterator_t>
Acts::BinnedSPGroup<external_spacepoint_t>::BinnedSPGroup(
    spacepoint_iterator_t spBegin, spacepoint_iterator_t spEnd,
    std::function<std::pair<Acts::Vector3, Acts::Vector2>(
        const external_spacepoint_t&, float, float, float)>
        globTool,
    std::shared_ptr<Acts::BinFinder<external_spacepoint_t>> botBinFinder,
    std::shared_ptr<Acts::BinFinder<external_spacepoint_t>> tBinFinder,
    std::unique_ptr<SpacePointGrid<external_spacepoint_t>> grid,
    const SeedfinderConfig<external_spacepoint_t>& _config) {
  auto config = _config.toInternalUnits();
  static_assert(
      std::is_same<
          typename std::iterator_traits<spacepoint_iterator_t>::value_type,
          const external_spacepoint_t*>::value,
      "Iterator does not contain type this class was templated with");

  // get region of interest (or full detector if configured accordingly)
  float phiMin = config.phiMin;
  float phiMax = config.phiMax;
  float zMin = config.zMin;
  float zMax = config.zMax;

  // sort by radius
  // add magnitude of beamPos to rMax to avoid excluding measurements
  // create number of bins equal to number of millimeters rMax
  // (worst case minR: configured minR + 1mm)
  size_t numRBins = (config.rMax + config.beamPos.norm());
  std::vector<std::vector<
      std::unique_ptr<const InternalSpacePoint<external_spacepoint_t>>>>
      rBins(numRBins);
  for (spacepoint_iterator_t it = spBegin; it != spEnd; it++) {
    if (*it == nullptr) {
      continue;
    }
    const external_spacepoint_t& sp = **it;
    const auto& [spPosition, variance] =
        globTool(sp, config.zAlign, config.rAlign, config.sigmaError);

    float spX = spPosition[0];
    float spY = spPosition[1];
    float spZ = spPosition[2];

    if (spZ > zMax || spZ < zMin) {
      continue;
    }
    float spPhi = std::atan2(spY, spX);
    if (spPhi > phiMax || spPhi < phiMin) {
      continue;
    }

    auto isp =
        std::make_unique<const InternalSpacePoint<external_spacepoint_t>>(
            sp, spPosition, config.beamPos, variance);
    // calculate r-Bin index and protect against overflow (underflow not
    // possible)
    size_t rIndex = isp->radius();
    // if index out of bounds, the SP is outside the region of interest
    if (rIndex >= numRBins) {
      continue;
    }
    rBins[rIndex].push_back(std::move(isp));
  }
  // fill rbins into grid such that each grid bin is sorted in r
  // space points with delta r < rbin size can be out of order
  for (auto& rbin : rBins) {
    for (auto& isp : rbin) {
      Acts::Vector2 spLocation(isp->phi(), isp->z());
      std::vector<
          std::unique_ptr<const InternalSpacePoint<external_spacepoint_t>>>&
          bin = grid->atPosition(spLocation);
      bin.push_back(std::move(isp));
    }
  }
  m_binnedSP = std::move(grid);
  m_bottomBinFinder = botBinFinder;
  m_topBinFinder = tBinFinder;
}
