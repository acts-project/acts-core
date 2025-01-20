// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/EventData/PropagationSummary.hpp"
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"
#include "ActsExamples/Propagation/PropagatorInterface.hpp"
#include "ActsExamples/Traccc/DetrayStore.hpp"

#include <detray/navigation/navigator.hpp>
#include <detray/propagator/actor_chain.hpp>
#include <detray/propagator/propagator.hpp>
#include <detray/test/utils/inspectors.hpp>
#include <detray/test/validation/material_validation_utils.hpp>

namespace ActsExamples {

/// Define the algebra type
using DetrayAlgebraType = typename Acts::DetrayHostDetector::algebra_type;

/// Type that holds the intersection information
using DetrayIntersection =
    detray::intersection2D<typename Acts::DetrayHostDetector::surface_type,
                           DetrayAlgebraType>;

/// Inspector that records all encountered surfaces
using DetrayObjectTracer =
    detray::navigation::object_tracer<DetrayIntersection, detray::dvector,
                                      detray::navigation::status::e_on_module,
                                      detray::navigation::status::e_on_portal>;

template <typename stepper_t, typename detray_store_t>
class DetrayPropagator : public PropagatorInterface {
 public:
  /// Configuration struct
  struct Config {
    /// The detray store
    std::shared_ptr<const detray_store_t> detrayStore = nullptr;
    /// Switch to sterile
    bool sterile = false;
  };

  /// Create a DetrayPropagator
  ///
  /// @param cfg configuration struct
  /// @param logger The logger instance
  DetrayPropagator(const Config& cfg,
                   std::unique_ptr<const Acts::Logger> logger =
                       Acts::getDefaultLogger("DetrayPropagator",
                                              Acts::Logging::INFO))
      : PropagatorInterface(), m_cfg(cfg), m_logger(std::move(logger)) {}

  ///@brief  Execute a propagation for charged particle parameters
  ///
  ///@param context The algorithm context
  ///@param cfg  The propagation algorithm configuration
  ///@param logger A logger wrapper instance
  ///@param startParameters The start parameters
  ///@return PropagationOutput
  Acts::Result<PropagationOutput> execute(
      const AlgorithmContext& context,
      [[maybe_unused]] const PropagationAlgorithm::Config& cfg,
      const Acts::Logger& logger,
      const Acts::BoundTrackParameters& startParameters) const final {
    // Get the geometry context form the algorithm context
    const auto& geoContext = context.geoContext;
    // Get the track information
    const Acts::Vector3 position = startParameters.position(geoContext);
    const Acts::Vector3 direction = startParameters.momentum().normalized();

    ACTS_VERBOSE("Starting propagation at " << position.transpose()
                                            << " with direction "
                                            << direction.transpose());

    // Now follow that ray with the same track and check, if we find
    // the same volumes and distances along the way
    detray::free_track_parameters<DetrayAlgebraType> track(
        {position.x(), position.y(), position.z()}, 0.f,
        {direction.x(), direction.y(), direction.z()},
        startParameters.charge());

    // Return material
    RecordedMaterial recordedMaterial;
    PropagationSummary summary(startParameters);

    if (!m_cfg.sterile) {
      /// Aggregation of multiple inspectors
      using DetrayInspector = detray::aggregate_inspector<DetrayObjectTracer>;

      // Navigation with inspection
      using DetrayNavigator =
          detray::navigator<Acts::DetrayHostDetector,
                            detray::navigation::default_cache_size,
                            DetrayInspector>;

      // Propagator with actor chain
      using MaterialTracer =
          detray::material_validator::material_tracer<Acts::ActsScalar,
                                                      detray::dvector>;
      using ActorChain = detray::actor_chain<detray::dtuple, MaterialTracer>;

      using Propagator = detray::propagator<
          stepper_t, DetrayNavigator, ActorChain>;

      typename Propagator::state pState(track,
                                             m_cfg.detrayStore->detector);

      Propagator propagator;

      MaterialTracer::state mTracerState{*m_cfg.detrayStore->memoryResource};
      auto actorStates = detray::tie(mTracerState);

      // Run the actual propagation
      propagator.propagate(pState, actorStates);

      // Retrieve navigation information
      auto& inspector = pState._navigation.inspector();
      auto& objectTracer = inspector.template get<DetrayObjectTracer>();

      // Translate the objects into the steps
      for (const auto& object : objectTracer.object_trace) {
        // Get the position of the object
        const auto& dposition = object.pos;
        const auto& sfDesription = object.intersection.sf_desc;
        const auto sf =
            detray::tracking_surface{m_cfg.detrayStore->detector, sfDesription};
        Acts::GeometryIdentifier geoID{sf.source()};
        // Create a step from the object
        Acts::detail::Step step;
        step.position = Acts::Vector3(dposition[0], dposition[1], dposition[2]);
        step.geoID = geoID;
        step.navDir = object.intersection.direction
                          ? Acts::Direction::Forward()
                          : Acts::Direction::Backward();
        summary.steps.emplace_back(step);
      }

      // Convert the material
      const auto& mRecord = mTracerState.get_material_record();
      recordedMaterial.materialInX0 = mRecord.sX0;
      recordedMaterial.materialInL0 = mRecord.sL0;

    } else {
      // Navigation with inspection
      using DetrayNavigator =
          detray::navigator<Acts::DetrayHostDetector,
                            detray::navigation::default_cache_size>;

      // Propagator with empty actor chain (for the moment)
      using Propagator =
          detray::propagator<stepper_t, DetrayNavigator, detray::actor_chain<>>;

      typename Propagator::state pState(track,
                                             m_cfg.detrayStore->detector);

      Propagator propagator;
      // Run the actual propagation
      propagator.propagate(pState);
    }

    return std::pair{std::move(summary), std::move(recordedMaterial)};
  }

 private:
  /// The detray detector store and memory resource
  Config m_cfg;

  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger = nullptr;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
