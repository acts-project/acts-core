// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFitting/RefittingAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/TrackFitting/RefittingCalibrator.hpp"
#include "ActsExamples/TrackFitting/TrackFitterFunction.hpp"

#include <algorithm>
#include <functional>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <system_error>
#include <utility>
#include <vector>

ActsExamples::RefittingAlgorithm::RefittingAlgorithm(Config config,
                                                     Acts::Logging::Level level)
    : ActsExamples::IAlgorithm("TrackFittingAlgorithm", level),
      m_cfg(std::move(config)) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing input tracks collection");
  }
  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Missing output tracks collection");
  }

  m_inputTracks.initialize(m_cfg.inputTracks);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ActsExamples::ProcessCode ActsExamples::RefittingAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  const auto& inputTracks = m_inputTracks(ctx);

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);

  // Perform the fit for each input track
  std::vector<Acts::SourceLink> trackSourceLinks;
  std::vector<const Acts::Surface*> surfSequence;
  RefittingCalibrator calibrator;

  auto itrack = 0ul;
  for (const auto& track : inputTracks) {
    // Check if you are not in picking mode
    if (m_cfg.pickTrack > -1 &&
        static_cast<std::size_t>(m_cfg.pickTrack) != itrack++) {
      continue;
    }

    if (!track.hasReferenceSurface()) {
      ACTS_VERBOSE("Skip track " << itrack << ": missing ref surface");
      continue;
    }

    TrackFitterFunction::GeneralFitterOptions options{
        ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
        &track.referenceSurface(),
        Acts::PropagatorPlainOptions(ctx.geoContext, ctx.magFieldContext)};
    options.doRefit = true;

    Acts::BoundTrackParameters initialParams(
        track.referenceSurface().getSharedPtr(), track.parameters(),
        track.covariance(), track.particleHypothesis());

    if (initialParams.covariance()) {
      for (auto i = 0ul; i < m_cfg.initialVarInflation.size(); ++i) {
        (*initialParams.covariance())(i, i) *= m_cfg.initialVarInflation.at(i);
      }
    }

    trackSourceLinks.clear();
    surfSequence.clear();

    for (auto state : track.trackStatesReversed()) {
      surfSequence.push_back(&state.referenceSurface());

      if (!state.hasCalibrated()) {
        continue;
      }

      auto sl = RefittingCalibrator::RefittingSourceLink{state};
      trackSourceLinks.push_back(Acts::SourceLink{sl});
    }

    if (surfSequence.empty()) {
      ACTS_WARNING("Empty track " << itrack << " found.");
      continue;
    }

    std::ranges::reverse(surfSequence);

    ACTS_VERBOSE("Initial parameters: "
                 << initialParams.fourPosition(ctx.geoContext).transpose()
                 << " -> " << initialParams.direction().transpose());

    ACTS_DEBUG("Invoke direct fitter for track " << itrack);
    auto result = (*m_cfg.fit)(trackSourceLinks, initialParams, options,
                               calibrator, surfSequence, tracks);

    if (result.ok()) {
      // Get the fit output object
      const auto& refittedTrack = result.value();
      if (refittedTrack.hasReferenceSurface()) {
        ACTS_VERBOSE("Refitted parameters for track " << itrack);
        ACTS_VERBOSE("  " << track.parameters().transpose());
      } else {
        ACTS_DEBUG("No refitted parameters for track " << itrack);
      }
    } else {
      ACTS_WARNING("Fit failed for track "
                   << itrack << " with error: " << result.error() << ", "
                   << result.error().message());
    }
  }

  std::stringstream ss;
  trackStateContainer->statistics().toStream(ss);
  ACTS_DEBUG(ss.str());

  ConstTrackContainer constTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackContainer)),
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*trackStateContainer))};

  m_outputTracks(ctx, std::move(constTracks));
  return ActsExamples::ProcessCode::SUCCESS;
}
