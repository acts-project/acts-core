// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <cmath>
#include <cstddef>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <system_error>
#include <utility>
#include <vector>

ActsExamples::TrackParamsEstimationAlgorithm::TrackParamsEstimationAlgorithm(
    ActsExamples::TrackParamsEstimationAlgorithm::Config cfg,
    Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("TrackParamsEstimationAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds input collection");
  }
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing track parameters output collection");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  if (not m_cfg.magneticField) {
    throw std::invalid_argument("Missing magnetic field");
  }

  m_inputSeeds.initialize(m_cfg.inputSeeds);
  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);

  // Set up the track parameters covariance (the same for all tracks)
  m_covariance(Acts::eBoundLoc0, Acts::eBoundLoc0) =
      m_cfg.initialVarInflation[Acts::eBoundLoc0] * cfg.sigmaLoc0 *
      m_cfg.sigmaLoc0;
  m_covariance(Acts::eBoundLoc1, Acts::eBoundLoc1) =
      m_cfg.initialVarInflation[Acts::eBoundLoc1] * cfg.sigmaLoc1 *
      m_cfg.sigmaLoc1;
  m_covariance(Acts::eBoundPhi, Acts::eBoundPhi) =
      m_cfg.initialVarInflation[Acts::eBoundPhi] * cfg.sigmaPhi *
      m_cfg.sigmaPhi;
  m_covariance(Acts::eBoundTheta, Acts::eBoundTheta) =
      m_cfg.initialVarInflation[Acts::eBoundTheta] * cfg.sigmaTheta *
      m_cfg.sigmaTheta;
  m_covariance(Acts::eBoundQOverP, Acts::eBoundQOverP) =
      m_cfg.initialVarInflation[Acts::eBoundQOverP] * cfg.sigmaQOverP *
      m_cfg.sigmaQOverP;
  m_covariance(Acts::eBoundTime, Acts::eBoundTime) =
      m_cfg.initialVarInflation[Acts::eBoundTime] * m_cfg.sigmaT0 *
      m_cfg.sigmaT0;
}

ActsExamples::ProcessCode ActsExamples::TrackParamsEstimationAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  auto const& seeds = m_inputSeeds(ctx);
  ACTS_VERBOSE("Read " << seeds.size() << " seeds");

  TrackParametersContainer trackParameters;
  trackParameters.reserve(seeds.size());

  auto bCache = m_cfg.magneticField->makeCache(ctx.magFieldContext);

  IndexSourceLink::SurfaceAccessor surfaceAccessor{*m_cfg.trackingGeometry};

  // Loop over all found seeds to estimate track parameters
  for (size_t iseed = 0; iseed < seeds.size(); ++iseed) {
    const auto& seed = seeds[iseed];
    // Get the bottom space point and its reference surface
    const auto bottomSP = seed.sp().front();
    if (bottomSP->sourceLinks().empty()) {
      ACTS_WARNING("Missing source link in the space point")
      continue;
    }
    const auto& sourceLink = bottomSP->sourceLinks()[0];
    const Acts::Surface* surface = surfaceAccessor(sourceLink);

    if (surface == nullptr) {
      ACTS_WARNING(
          "Surface from source link is not found in the tracking geometry");
      continue;
    }

    // Get the magnetic field at the bottom space point
    auto fieldRes = m_cfg.magneticField->getField(
        {bottomSP->x(), bottomSP->y(), bottomSP->z()}, bCache);
    if (!fieldRes.ok()) {
      ACTS_ERROR("Field lookup error: " << fieldRes.error());
      return ProcessCode::ABORT;
    }
    Acts::Vector3 field = *fieldRes;

    // Estimate the track parameters from seed
    auto optParams = Acts::estimateTrackParamsFromSeed(
        ctx.geoContext, seed.sp().begin(), seed.sp().end(), *surface, field,
        m_cfg.bFieldMin, logger());
    if (not optParams.has_value()) {
      ACTS_WARNING("Estimation of track parameters for seed " << iseed
                                                              << " failed.");
      continue;
    } else {
      const auto& params = optParams.value();
      trackParameters.emplace_back(surface->getSharedPtr(), params,
                                   m_covariance,
                                   Acts::ParticleHypothesis::pion());
    }
  }

  ACTS_VERBOSE("Estimated " << trackParameters.size() << " track parameters");

  m_outputTrackParameters(ctx, std::move(trackParameters));
  return ProcessCode::SUCCESS;
}
