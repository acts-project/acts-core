// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Alignment/AlignmentAlgorithm.hpp"

#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

ActsExamples::AlignmentAlgorithm::AlignmentAlgorithm(Config cfg,
                                                     Acts::Logging::Level level)
    : ActsExamples::BareAlgorithm("AlignmentAlgorithm", level),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing input measurement collection");
  }
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing input source links collection");
  }
  if (m_cfg.inputProtoTracks.empty()) {
    throw std::invalid_argument("Missing input proto tracks collection");
  }
  if (m_cfg.inputInitialTrackParameters.empty()) {
    throw std::invalid_argument(
        "Missing input initial track parameters collection");
  }
  if (m_cfg.outputAlignmentParameters.empty()) {
    throw std::invalid_argument(
        "Missing output alignment parameters collection");
  }
}

ActsExamples::ProcessCode ActsExamples::AlignmentAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input data
  const auto& measurements =
      ctx.eventStore.get<MeasurementContainer>(m_cfg.inputMeasurements);
  const auto& sourceLinks =
      ctx.eventStore.get<IndexSourceLinkContainer>(m_cfg.inputSourceLinks);
  const auto& protoTracks =
      ctx.eventStore.get<ProtoTrackContainer>(m_cfg.inputProtoTracks);
  const auto& initialParameters = ctx.eventStore.get<TrackParametersContainer>(
      m_cfg.inputInitialTrackParameters);

  // Consistency cross checks
  if (protoTracks.size() != initialParameters.size()) {
    ACTS_FATAL("Inconsistent number of proto tracks and parameters");
    return ProcessCode::ABORT;
  }

  size_t numTracksUsed = protoTracks.size();
  if (m_cfg.maxNumTracks > 0 and
      m_cfg.maxNumTracks < static_cast<int>(protoTracks.size())) {
    numTracksUsed = m_cfg.maxNumTracks;
  }

  // Prepare the input track collection
  std::vector<std::vector<IndexSourceLink>> sourceLinkTrackContainer;
  sourceLinkTrackContainer.reserve(numTracksUsed);
  std::vector<IndexSourceLink> trackSourceLinks;
  for (std::size_t itrack = 0; itrack < numTracksUsed; ++itrack) {
    // The list of hits and the initial start parameters
    const auto& protoTrack = protoTracks[itrack];

    // Clear & reserve the right size
    trackSourceLinks.clear();
    trackSourceLinks.reserve(protoTrack.size());

    // Fill the source links via their indices from the container
    for (auto hitIndex : protoTrack) {
      auto sourceLink = sourceLinks.nth(hitIndex);
      if (sourceLink == sourceLinks.end()) {
        ACTS_FATAL("Proto track " << itrack << " contains invalid hit index"
                                  << hitIndex);
        return ProcessCode::ABORT;
      }
      trackSourceLinks.push_back(*sourceLink);
    }
    sourceLinkTrackContainer.push_back(trackSourceLinks);
  }

  // Prepare the output for alignment parameters
  std::unordered_map<Acts::DetectorElementBase*, Acts::Transform3>
      alignedParameters;

  // Construct a perigee surface as the target surface for the fitter
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{0., 0., 0.});

  // Set the KalmanFitter options
  TrackFitterOptions kfOptions(
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
      MeasurementCalibrator(measurements), Acts::VoidOutlierFinder(),
      Acts::LoggerWrapper{logger()}, Acts::PropagatorPlainOptions(),
      &(*pSurface));

  // Set the alignment options
  ActsAlignment::AlignmentOptions<TrackFitterOptions> alignOptions(
      kfOptions, m_cfg.alignedTransformUpdater, Acts::LoggerWrapper{logger()},
      m_cfg.alignedDetElements, m_cfg.chi2ONdfCutOff, m_cfg.deltaChi2ONdfCutOff,
      m_cfg.maxNumIterations);

  ACTS_DEBUG("Invoke track-based alignment with " << numTracksUsed
                                                  << " input tracks");
  auto result =
      (*m_cfg.align)(sourceLinkTrackContainer, initialParameters, alignOptions);
  if (result.ok()) {
    const auto& alignOutput = result.value();
    alignedParameters = alignOutput.alignedParameters;
    ACTS_VERBOSE(
        "Alignment finished with deltaChi2 = " << result.value().deltaChi2);
  } else {
    ACTS_WARNING("Alignment failed with " << result.error());
  }

  // add alignment parameters to event store
  ctx.eventStore.add(m_cfg.outputAlignmentParameters,
                     std::move(alignedParameters));
  return ActsExamples::ProcessCode::SUCCESS;
}
