// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"

namespace ActsExamples {

///@brief Propagator wrapper
///
/// This class encapsulates a propagator instance and exposes
/// propagation functions for neutral and charged particles.
/// @note This interface exists to decouple the concrete propagators from
///       the @c PropagationAlgorithm
class PropagatorInterface {
 public:
  virtual ~PropagatorInterface() = default;

  ///@brief  Execute a propagation for charged particle parameters
  ///
  ///@param context The algorithm context
  ///@param cfg  The propagation algorithm configuration
  ///@param logger A logger wrapper instance
  ///@param startParameters The start parameters
  ///@return PropagationOutput
  virtual Acts::Result<PropagationOutput> execute(
      const AlgorithmContext& context, const PropagationAlgorithm::Config& cfg,
      const Acts::Logger& logger,
      const Acts::BoundTrackParameters& startParameters) const = 0;
};

///@brief Concrete instance of a propagator
///
/// This class implements a wrapped version of a concrete propagator.
/// It can be instantiated to comply with @c PropagatorInterface.
/// Setting up a @c PropagationAlgorithm with this looks like:
/// ```cpp
/// config.propagatorImpl =
///     std::make_shared<ActsExamples::ConcretePropagator<Propagator>>(
///         std::move(propagator));
/// ```
///@tparam propagator_t The concrete propagator to instantiate with
template <typename propagator_t>
class ConcretePropagator : public PropagatorInterface {
 public:
  explicit ConcretePropagator(propagator_t propagator)
      : m_propagator{std::move(propagator)} {}

  Acts::Result<PropagationOutput> execute(
      const AlgorithmContext& context, const PropagationAlgorithm::Config& cfg,
      const Acts::Logger& logger,
      const Acts::BoundTrackParameters& startParameters) const override {
    ACTS_DEBUG("Test propagation/extrapolation starts");

    PropagationSummary summary(startParameters);
    RecordedMaterial recordedMaterial;

    // The step length logger for testing & end of world aborter
    using MaterialInteractor = Acts::MaterialInteractor;
    using SteppingLogger = Acts::detail::SteppingLogger;
    using EndOfWorld = Acts::EndOfWorldReached;

    // Actor list
    using ActorList =
        Acts::ActorList<SteppingLogger, MaterialInteractor, EndOfWorld>;
    using PropagatorOptions =
        typename propagator_t::template Options<ActorList>;

    PropagatorOptions options(context.geoContext, context.magFieldContext);
    // Activate loop protection at some pt value
    options.loopProtection =
        startParameters.transverseMomentum() < cfg.ptLoopers;

    // Switch the material interaction on/off & eventually into logging mode
    auto& mInteractor = options.actorList.template get<MaterialInteractor>();
    mInteractor.multipleScattering = cfg.multipleScattering;
    mInteractor.energyLoss = cfg.energyLoss;
    mInteractor.recordInteractions = cfg.recordMaterialInteractions;

    // Switch the logger to sterile, e.g. for timing checks
    auto& sLogger = options.actorList.template get<SteppingLogger>();
    sLogger.sterile = cfg.sterileLogger;
    // Set a maximum step size
    options.stepping.maxStepSize = cfg.maxStepSize;

    auto state = m_propagator.makeState(startParameters, options);

    auto resultInit = m_propagator.initialize(state);
    if (!resultInit.ok()) {
      return resultInit.error();
    }

    // Propagate using the propagator
    auto resultTmp = m_propagator.propagate(state);
    if (!resultTmp.ok()) {
      return resultTmp.error();
    }

    // Collect internal stepping information
    summary.nStepTrials = state.stepping.nStepTrials;

    auto result =
        m_propagator.makeResult(std::move(state), resultTmp, options, true);
    if (!result.ok()) {
      return result.error();
    }
    auto& resultValue = result.value();

    // Collect general summary information
    summary.nSteps = resultValue.steps;
    summary.pathLength = resultValue.pathLength;

    // Collect the steps
    auto& steppingResults =
        resultValue.template get<SteppingLogger::result_type>();
    summary.steps = std::move(steppingResults.steps);

    summary.statistics = resultValue.statistics;

    // Also set the material recording result - if configured
    if (cfg.recordMaterialInteractions) {
      auto materialResult =
          resultValue.template get<MaterialInteractor::result_type>();
      recordedMaterial = std::move(materialResult);
    }

    return std::pair{std::move(summary), std::move(recordedMaterial)};
  }

 private:
  propagator_t m_propagator;
};

}  // namespace ActsExamples
