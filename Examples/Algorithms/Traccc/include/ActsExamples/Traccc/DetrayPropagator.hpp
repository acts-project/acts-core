// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"
#include "ActsExamples/Propagation/PropagatorInterface.hpp"
#include "ActsExamples/Traccc/DetrayStore.hpp"
#include "ActsExamples/EventData/PropagationSummary.hpp"

#include <detray/navigation/navigator.hpp>
#include <detray/utils/inspectors.hpp>

namespace ActsExamples {

  /// Type that holds the intersection information
  using DetrayIntersection =
      detray::intersection2D<typename Acts::DetrayDetector::surface_type,
                             detray::cmath<detray::scalar>>;

  /// Inspector that records all encountered surfaces
  using DetrayObjectTracer = detray::navigation::object_tracer<
      DetrayIntersection, detray::dvector, detray::navigation::status::e_on_module,
      detray::navigation::status::e_on_portal>;

  /// Inspector that prints the navigator state from within the navigator's
  /// method calls (cannot be done with an actor)
  using DetrayPrintInspector = detray::navigation::print_inspector;

template <typename propagator_t, typename detray_store_t>
class DetrayPropagator : public PropagatorInterface {
 public:

  /// Create a DetrayPropagator
  ///
  /// @param propagator The actual detray propagator to wrap
  /// @param detrayStore The detray store to access the detector
  /// @param logger The logger instance
  DetrayPropagator(propagator_t&& propagator,
                   std::shared_ptr<const detray_store_t> detrayStore,
                   std::unique_ptr<const Acts::Logger> logger =
                       Acts::getDefaultLogger("DetrayPropagator",
                                              Acts::Logging::INFO))
      : PropagatorInterface(),
        m_propagator(std::move(propagator)),
        m_detrayStore(std::move(detrayStore)),
        m_logger(std::move(logger)) {}

  ///@brief  Execute a propagation for charged particle parameters
  ///
  ///@param context The algorithm context
  ///@param cfg  The propagation algorithm configuration
  ///@param logger A logger wrapper instance
  ///@param startParameters The start parameters
  ///@return PropagationOutput
  Acts::Result<PropagationOutput> execute(
      const AlgorithmContext& context, const PropagationAlgorithm::Config& cfg,
      const Acts::Logger& logger,
      const Acts::BoundTrackParameters& startParameters) const final {
    const auto& geoContext = context.geoContext;
    // Get the detector
    const Acts::Vector3 position = startParameters.position(geoContext);
    const Acts::Vector3 direction = startParameters.momentum().normalized();

    ACTS_VERBOSE("Starting propagation at " << position.transpose()
                                         << " with direction "
                                         << direction.transpose());

    // Now follow that ray with the same track and check, if we find
    // the same volumes and distances along the way
    detray::free_track_parameters<detray::cmath<detray::scalar>> track(
        {position.x(), position.y(), position.z()}, 0.f,
        {direction.x(), direction.y(), direction.z()}, -1.f);

    typename decltype(m_propagator)::state propagation(track,
                                                       m_detrayStore->detector);

    // Run the actual propagation
    m_propagator.propagate(propagation);

    // Retrieve navigation information
    auto& inspector = propagation._navigation.inspector();
    auto& objectTracer = inspector.template get<DetrayObjectTracer>();
    auto& debugPrinter = inspector.template get<DetrayPrintInspector>();

    PropagationSummary summary(startParameters);
    for (const auto& object : objectTracer.object_trace) {
      //const auto& intersection = objectTracer[intr_idx].intersection;
       const auto& dposition = object.pos;
       Acts::detail::Step step;
       step.position = Acts::Vector3(dposition[0], dposition[1], dposition[2]);
       summary.steps.emplace_back(step);
    }


    RecordedMaterial recordedMaterial;

    return std::pair{std::move(summary), std::move(recordedMaterial)};
  }

 private:
  /// The propagator @todo fix when propagate() method is const in detray
  mutable propagator_t m_propagator;

  /// The detray detector store and memory ressource
  std::shared_ptr<const detray_store_t> m_detrayStore = nullptr;

  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger = nullptr;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples