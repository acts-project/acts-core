// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/detail/CuboidalDetectorHelper.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Logger.hpp"

using namespace Acts::UnitLiterals;

namespace Acts {
class Surface;
}  // namespace Acts

Acts::GeometryContext geoContext;
Acts::MagneticFieldContext mfContext;

/// A simple action to record the navigation state
struct StateRecorder {
  using result_type = std::vector<Acts::Experimental::DetectorNavigator::State>;

  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  void act(propagator_state_t& state, const stepper_t& /*stepper*/,
           const navigator_t& /*navigator*/, result_type& result,
           const Acts::Logger& /*logger*/) const {
    result.push_back(state.navigation);
  }
};

BOOST_AUTO_TEST_SUITE(DetectorNavigator)

// Initialization tests
BOOST_AUTO_TEST_CASE(DetectorNavigatorTestsInitialization) {
  auto bounds = std::make_unique<Acts::CuboidVolumeBounds>(3, 3, 3);
  auto volume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), geoContext,
      "volume", Acts::Transform3::Identity(), std::move(bounds), {}, {},
      Acts::Experimental::tryNoVolumes(),
      Acts::Experimental::tryAllPortalsAndSurfaces());
  volume->assignGeometryId(1);
  auto detector = Acts::Experimental::Detector::makeShared(
      "detector", {volume}, Acts::Experimental::tryRootVolumes());

  using Stepper = Acts::StraightLineStepper;
  using Navigator = Acts::Experimental::DetectorNavigator;
  using Propagator = Acts::Propagator<Stepper, Navigator>;
  using ActorList = Acts::ActorList<>;
  using PropagatorOptions = Propagator::Options<ActorList>;

  PropagatorOptions options(geoContext, mfContext);

  Stepper stepper;

  Acts::Vector4 pos(-2, 0, 0, 0);
  Acts::CurvilinearTrackParameters start(pos, 0_degree, 90_degree, 1_e / 1_GeV,
                                         std::nullopt,
                                         Acts::ParticleHypothesis::electron());

  //
  // (1) Test for inactivity
  //
  // Run without anything present
  {
    Navigator::Config navCfg;
    navCfg.resolveSensitive = false;
    navCfg.resolveMaterial = false;
    navCfg.resolvePassive = false;

    Navigator navigator(navCfg);

    Propagator propagator(stepper, navigator);

    auto state = propagator.makeState(start, options);

    BOOST_CHECK_THROW(propagator.initialize(state), std::invalid_argument);
  }

  // Run with geometry but without resolving
  {
    Acts::Experimental::DetectorNavigator::Config navCfg;
    navCfg.resolveSensitive = false;
    navCfg.resolveMaterial = false;
    navCfg.resolvePassive = false;
    navCfg.detector = detector.get();

    Acts::Experimental::DetectorNavigator navigator(navCfg);

    Acts::Propagator<Acts::StraightLineStepper,
                     Acts::Experimental::DetectorNavigator>
        propagator(stepper, navigator);

    auto state = propagator.makeState(start, options);

    navigator.initialize(state.navigation, stepper.position(state.stepping),
                         stepper.direction(state.stepping),
                         state.options.direction);

    navigator.nextTarget(state.navigation, stepper.position(state.stepping),
                         stepper.direction(state.stepping));
    auto preStepState = state.navigation;
    BOOST_CHECK_EQUAL(preStepState.currentSurface, nullptr);
    BOOST_CHECK_EQUAL(preStepState.currentPortal, nullptr);
  }

  //
  // (2) Initialization tests
  //
  // Run from endOfWorld
  {
    Acts::Vector4 posEoW(-20, 0, 0, 0);
    Acts::CurvilinearTrackParameters startEoW(
        posEoW, 0_degree, 90_degree, 1_e / 1_GeV, std::nullopt,
        Acts::ParticleHypothesis::electron());

    Acts::Experimental::DetectorNavigator::Config navCfg;
    navCfg.detector = detector.get();

    Acts::Experimental::DetectorNavigator navigator(navCfg);

    Acts::Propagator<Acts::StraightLineStepper,
                     Acts::Experimental::DetectorNavigator>
        propagator(stepper, navigator);

    auto state = propagator.makeState(startEoW, options);

    BOOST_CHECK_THROW(propagator.initialize(state), std::invalid_argument);
  }

  // Initialize properly
  {
    Acts::Experimental::DetectorNavigator::Config navCfg;
    navCfg.detector = detector.get();

    Acts::Experimental::DetectorNavigator navigator(navCfg);

    Acts::Propagator<Acts::StraightLineStepper,
                     Acts::Experimental::DetectorNavigator>
        propagator(stepper, navigator);

    auto state = propagator.makeState(start, options);

    navigator.initialize(state.navigation, stepper.position(state.stepping),
                         stepper.direction(state.stepping),
                         state.options.direction);

    auto initState = state.navigation;
    BOOST_CHECK_EQUAL(initState.currentDetector, detector.get());
    BOOST_CHECK_EQUAL(
        initState.currentVolume,
        detector->findDetectorVolume(geoContext, start.position()));
    BOOST_CHECK_EQUAL(initState.currentSurface, nullptr);
    BOOST_CHECK_EQUAL(initState.currentPortal, nullptr);
    BOOST_CHECK_EQUAL(initState.surfaceCandidates.size(), 1u);
  }
}

// Stadard forward and backward propagation
// through cubic volumes with planar surfaces
// and no surprises
BOOST_AUTO_TEST_CASE(DetectorNavigatorTestsForwardBackward) {
  // Construct a cubic detector with 3 volumes
  Acts::RotationMatrix3 rotation;
  double angle = 90_degree;
  Acts::Vector3 xPos(cos(angle), 0., sin(angle));
  Acts::Vector3 yPos(0., 1., 0.);
  Acts::Vector3 zPos(-sin(angle), 0., cos(angle));
  rotation.col(0) = xPos;
  rotation.col(1) = yPos;
  rotation.col(2) = zPos;

  auto bounds1 = std::make_unique<Acts::CuboidVolumeBounds>(3, 3, 3);
  auto transform1 = Acts::Transform3::Identity();
  auto surface1 = Acts::Surface::makeShared<Acts::PlaneSurface>(
      transform1 * Acts::Transform3(rotation),
      std::make_shared<Acts::RectangleBounds>(2, 2));
  auto volume1 = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), geoContext,
      "volume1", transform1, std::move(bounds1), {surface1}, {},
      Acts::Experimental::tryNoVolumes(),
      Acts::Experimental::tryAllPortalsAndSurfaces());

  auto bounds2 = std::make_unique<Acts::CuboidVolumeBounds>(3, 3, 3);
  auto transform2 =
      Acts::Transform3::Identity() * Acts::Translation3(Acts::Vector3(6, 0, 0));
  auto surface2 = Acts::Surface::makeShared<Acts::PlaneSurface>(
      transform2 * Acts::Transform3(rotation),
      std::make_shared<Acts::RectangleBounds>(2, 2));
  auto volume2 = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), geoContext,
      "volume2", transform2, std::move(bounds2), {surface2}, {},
      Acts::Experimental::tryNoVolumes(),
      Acts::Experimental::tryAllPortalsAndSurfaces());

  auto bounds3 = std::make_unique<Acts::CuboidVolumeBounds>(3, 3, 3);
  auto transform3 = Acts::Transform3::Identity() *
                    Acts::Translation3(Acts::Vector3(12, 0, 0));
  auto surface3 = Acts::Surface::makeShared<Acts::PlaneSurface>(
      transform3 * Acts::Transform3(rotation),
      std::make_shared<Acts::RectangleBounds>(2, 2));
  auto volume3 = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), geoContext,
      "volume3", transform3, std::move(bounds3), {surface3}, {},
      Acts::Experimental::tryNoVolumes(),
      Acts::Experimental::tryAllPortalsAndSurfaces());

  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
      detectorVolumes = {volume1, volume2, volume3};

  auto portalContainer =
      Acts::Experimental::detail::CuboidalDetectorHelper::connect(
          geoContext, detectorVolumes, Acts::AxisDirection::AxisX, {},
          Acts::Logging::VERBOSE);

  // Make sure that the geometry ids are
  // independent of the potential Id generation
  // changes
  int id = 1;

  // Volume ids: 1-3
  for (auto& volume : detectorVolumes) {
    volume->assignGeometryId(id);
    id++;
  }
  // Intervolume portal ids: 6,7,10,11
  for (auto& volume : detectorVolumes) {
    for (auto& port : volume->portalPtrs()) {
      if (port->surface().geometryId() == 0) {
        port->surface().assignGeometryId(id);
        id++;
      }
    }
  }
  // Surface ids: 12-14
  for (auto& surf : {surface1, surface2, surface3}) {
    surf->assignGeometryId(id);
    id++;
  }

  auto detector = Acts::Experimental::Detector::makeShared(
      "cubicDetector", detectorVolumes, Acts::Experimental::tryRootVolumes());

  using Stepper = Acts::StraightLineStepper;
  using Navigator = Acts::Experimental::DetectorNavigator;
  using Propagator = Acts::Propagator<Stepper, Navigator>;
  using ActorList = Acts::ActorList<StateRecorder, Acts::EndOfWorldReached>;
  using PropagatorOptions = Propagator::Options<ActorList>;

  Navigator::Config navCfg;
  navCfg.detector = detector.get();

  Stepper stepper;

  Navigator navigator(navCfg,
                      Acts::getDefaultLogger("DetectorNavigator",
                                             Acts::Logging::Level::VERBOSE));

  PropagatorOptions options(geoContext, mfContext);
  options.direction = Acts::Direction::Forward();

  Propagator propagator(
      stepper, navigator,
      Acts::getDefaultLogger("Propagator", Acts::Logging::Level::VERBOSE));

  // Forward and backward propagation
  // should be consistent between each other
  Acts::Vector4 posFwd(-2, 0, 0, 0);
  Acts::CurvilinearTrackParameters startFwd(
      posFwd, 0_degree, 90_degree, 1_e / 1_GeV, std::nullopt,
      Acts::ParticleHypothesis::electron());

  auto resultFwd = propagator.propagate(startFwd, options).value();
  auto statesFwd = resultFwd.get<StateRecorder::result_type>();

  options.direction = Acts::Direction::Backward();

  Acts::Vector4 posBwd(14, 0, 0, 0);
  Acts::CurvilinearTrackParameters startBwd(
      posBwd, 0_degree, 90_degree, 1_e / 1_GeV, std::nullopt,
      Acts::ParticleHypothesis::electron());

  auto resultBwd = propagator.propagate(startBwd, options).value();
  auto statesBwd = resultBwd.get<StateRecorder::result_type>();

  // 7 steps to reach the end of world
  // + 1 recording in the post-step
  // + 1 recording before the stepping loop
  BOOST_CHECK_EQUAL(statesFwd.size(), 8u);
  BOOST_CHECK_EQUAL(statesFwd.size(), statesBwd.size());
  BOOST_CHECK_EQUAL(statesFwd[0].surfaceCandidates.size(), 2u);
  BOOST_CHECK_EQUAL(statesBwd[0].surfaceCandidates.size(), 2u);

  // Action list call before the first step
  // Starting in the volume1
  BOOST_CHECK_EQUAL(statesFwd[0].currentVolume->geometryId(), 1);
  BOOST_CHECK_EQUAL(statesFwd[0].currentSurface, nullptr);
  BOOST_CHECK_EQUAL(statesFwd[0].currentPortal, nullptr);

  // Step to the surface inside volume1
  BOOST_CHECK_EQUAL(statesFwd[1].currentVolume->geometryId(), 1);
  BOOST_CHECK_EQUAL(statesFwd[1].currentSurface->geometryId(), 12);
  BOOST_CHECK_EQUAL(statesFwd[1].currentPortal, nullptr);

  // Step to the volume1|volume2 boundary (portal has witched volume id)
  BOOST_CHECK_EQUAL(statesFwd[2].currentVolume->geometryId(), 2);
  BOOST_CHECK_EQUAL(statesFwd[2].currentSurface,
                    &(statesFwd[2].currentPortal->surface()));
  BOOST_CHECK_EQUAL(statesFwd[2].currentPortal->surface().geometryId(), 7);

  // Step to the surface inside volume2
  BOOST_CHECK_EQUAL(statesFwd[3].currentVolume->geometryId(), 2);
  BOOST_CHECK_EQUAL(statesFwd[3].currentSurface->geometryId(), 13);
  BOOST_CHECK_EQUAL(statesFwd[3].currentPortal, nullptr);

  // Step to the volume2|volume3 boundary - volume has switched
  BOOST_CHECK_EQUAL(statesFwd[4].currentVolume->geometryId(), 3);
  BOOST_CHECK_EQUAL(statesFwd[4].currentSurface,
                    &(statesFwd[4].currentPortal->surface()));
  BOOST_CHECK_EQUAL(statesFwd[4].currentPortal->surface().geometryId(), 10);

  // Step to the surface inside volume3
  BOOST_CHECK_EQUAL(statesFwd[5].currentVolume->geometryId(), 3);
  BOOST_CHECK_EQUAL(statesFwd[5].currentSurface->geometryId(), 14);
  BOOST_CHECK_EQUAL(statesFwd[5].currentPortal, nullptr);

  // Step to the volume3|endOfWorld boundary
  BOOST_CHECK_EQUAL(statesFwd[6].currentVolume, nullptr);
  BOOST_CHECK_EQUAL(statesFwd[6].currentSurface,
                    &(statesFwd[6].currentPortal->surface()));
  BOOST_CHECK_EQUAL(statesFwd[6].currentPortal->surface().geometryId(), 11);

  // Step to the end of world
  BOOST_CHECK(navigator.endOfWorldReached(statesFwd[7]));
  BOOST_CHECK(navigator.endOfWorldReached(statesBwd[7]));

  // Action list call before the first step
  // Starting in the volume3
  BOOST_CHECK_EQUAL(statesBwd[6].currentVolume, nullptr);
  BOOST_CHECK_EQUAL(statesBwd[6].currentSurface,
                    &(statesBwd[6].currentPortal->surface()));
  BOOST_CHECK_EQUAL(statesBwd[6].currentPortal->surface().geometryId(), 6);

  // Step to the surface inside volume1
  BOOST_CHECK_EQUAL(statesBwd[5].currentVolume->geometryId(), 1);
  BOOST_CHECK_EQUAL(statesBwd[5].currentSurface->geometryId(), 12);
  BOOST_CHECK_EQUAL(statesBwd[5].currentPortal, nullptr);

  // Step to the volume1|volume2 boundary / preStep not yet set
  BOOST_CHECK_EQUAL(statesBwd[4].currentVolume->geometryId(), 1);
  BOOST_CHECK_EQUAL(statesBwd[4].currentSurface,
                    &(statesBwd[4].currentPortal->surface()));
  BOOST_CHECK_EQUAL(statesBwd[4].currentPortal->surface().geometryId(), 7);

  // Step to the surface inside volume2
  BOOST_CHECK_EQUAL(statesBwd[3].currentVolume->geometryId(), 2);
  BOOST_CHECK_EQUAL(statesBwd[3].currentSurface->geometryId(), 13);
  BOOST_CHECK_EQUAL(statesBwd[3].currentPortal, nullptr);

  // Step to the volume2|volume3 boundary / pre-step not yet set
  BOOST_CHECK_EQUAL(statesBwd[2].currentVolume->geometryId(), 2);
  BOOST_CHECK_EQUAL(statesBwd[2].currentSurface,
                    &(statesBwd[2].currentPortal->surface()));
  BOOST_CHECK_EQUAL(statesBwd[2].currentPortal->surface().geometryId(), 10);
  BOOST_CHECK_EQUAL(statesBwd[2].surfaceCandidates.size(), 2u);

  // Step to the surface inside volume3
  BOOST_CHECK_EQUAL(statesBwd[1].currentVolume->geometryId(), 3);
  BOOST_CHECK_EQUAL(statesBwd[1].currentSurface->geometryId(), 14);
  BOOST_CHECK_EQUAL(statesBwd[1].currentPortal, nullptr);

  // Step to the volume3|endOfWorld boundary
  BOOST_CHECK_EQUAL(statesBwd[0].currentVolume->geometryId(), 3);
  BOOST_CHECK_EQUAL(statesBwd[0].currentSurface, nullptr);
  BOOST_CHECK_EQUAL(statesBwd[0].currentPortal, nullptr);
}

// Check how the navigator handles the case
// when the same surface may be reached
// in multiple points
BOOST_AUTO_TEST_CASE(DetectorNavigatorTestsAmbiguity) {
  // Construct a cubic detector with a cylindrical surface
  auto bounds = std::make_unique<Acts::CuboidVolumeBounds>(10, 10, 10);
  auto surface = Acts::Surface::makeShared<Acts::CylinderSurface>(
      Acts::Transform3::Identity(),
      std::make_shared<Acts::CylinderBounds>(4, 9));
  auto volume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), geoContext,
      "volume", Acts::Transform3::Identity(), std::move(bounds), {surface}, {},
      Acts::Experimental::tryNoVolumes(),
      Acts::Experimental::tryAllPortalsAndSurfaces());

  volume->assignGeometryId(1);
  surface->assignGeometryId(2);
  int id = 3;
  for (auto& port : volume->portalPtrs()) {
    port->surface().assignGeometryId(id);
    id++;
  }

  auto detector = Acts::Experimental::Detector::makeShared(
      "detector", {volume}, Acts::Experimental::tryRootVolumes());

  using Stepper = Acts::StraightLineStepper;
  using Navigator = Acts::Experimental::DetectorNavigator;
  using Propagator = Acts::Propagator<Stepper, Navigator>;
  using ActorList = Acts::ActorList<StateRecorder, Acts::EndOfWorldReached>;
  using PropagatorOptions = Propagator::Options<ActorList>;

  Navigator::Config navCfg;
  navCfg.detector = detector.get();

  Stepper stepper;

  Navigator navigator(navCfg,
                      Acts::getDefaultLogger("DetectorNavigator",
                                             Acts::Logging::Level::VERBOSE));

  PropagatorOptions options(geoContext, mfContext);
  options.direction = Acts::Direction::Forward();

  Propagator propagator(
      stepper, navigator,
      Acts::getDefaultLogger("Propagator", Acts::Logging::Level::VERBOSE));

  // Depending on the direction, the same surface
  // may be reached in different points
  Acts::Vector4 pos(0, 0, 0, 0);
  Acts::CurvilinearTrackParameters start(pos, 0_degree, 90_degree, 1_e / 1_GeV,
                                         std::nullopt,
                                         Acts::ParticleHypothesis::electron());

  // Has to properly handle propagation in the
  // forward and backward direction
  auto resultFwd = propagator.propagate(start, options).value();
  auto statesFwd = resultFwd.get<StateRecorder::result_type>();

  options.direction = Acts::Direction::Backward();

  auto resultBwd = propagator.propagate(start, options).value();
  auto statesBwd = resultBwd.get<StateRecorder::result_type>();

  // 3 steps to reach the end of world
  // + 1 recording in the post-step
  // + 1 recording before the stepping loop
  BOOST_CHECK_EQUAL(statesFwd.size(), 4u);
  BOOST_CHECK_EQUAL(statesFwd.size(), statesBwd.size());
  BOOST_CHECK_EQUAL(statesFwd[0].surfaceCandidates.size(), 2u);
  BOOST_CHECK_EQUAL(statesBwd[0].surfaceCandidates.size(), 2u);

  // Action list call before the first step
  // Starting in the volume
  BOOST_CHECK_EQUAL(statesFwd[0].currentVolume->geometryId(), 1);
  BOOST_CHECK_EQUAL(statesFwd[0].currentSurface, nullptr);
  BOOST_CHECK_EQUAL(statesFwd[0].currentPortal, nullptr);

  // Step to the cylindrical surface
  BOOST_CHECK_EQUAL(statesFwd[1].currentVolume->geometryId(), 1);
  BOOST_CHECK_EQUAL(statesFwd[1].currentSurface->geometryId(), 2);
  BOOST_CHECK_EQUAL(statesFwd[1].currentPortal, nullptr);
  BOOST_CHECK_EQUAL(statesFwd[1].surfaceCandidates.size(), 2u);
  CHECK_CLOSE_REL(statesFwd[1].position.x(), 4, 1e-6);

  // Step to the volume|endOfWorld boundary
  BOOST_CHECK_EQUAL(statesFwd[2].currentVolume, nullptr);
  BOOST_CHECK_EQUAL(statesFwd[2].currentSurface,
                    &(statesFwd[2].currentPortal->surface()));
  BOOST_CHECK_EQUAL(statesFwd[2].currentPortal->surface().geometryId(), 6);

  // Step to the end of world
  BOOST_CHECK(navigator.endOfWorldReached(statesFwd[3]));
  BOOST_CHECK(navigator.endOfWorldReached(statesBwd[3]));

  // Step to the endOfWorld|volume boundary
  BOOST_CHECK_EQUAL(statesBwd[2].currentVolume, nullptr);
  BOOST_CHECK_EQUAL(statesBwd[2].currentSurface,
                    &(statesBwd[2].currentPortal->surface()));
  BOOST_CHECK_EQUAL(statesBwd[2].currentPortal->surface().geometryId(), 5);

  // Step to the cylindrical surface
  BOOST_CHECK_EQUAL(statesBwd[1].currentVolume->geometryId(), 1);
  BOOST_CHECK_EQUAL(statesBwd[1].currentSurface->geometryId(), 2);
  BOOST_CHECK_EQUAL(statesBwd[1].currentPortal, nullptr);
  CHECK_CLOSE_REL(statesBwd[1].position.x(), -4, 1e-6);

  // Action list call before the first step
  // Starting in the volume
  BOOST_CHECK_EQUAL(statesBwd[0].currentVolume->geometryId(), 1);
  BOOST_CHECK_EQUAL(statesBwd[0].currentSurface, nullptr);
  BOOST_CHECK_EQUAL(statesBwd[0].currentPortal, nullptr);
}

// Check how the navigator handles the case
// when the same surface has multiple valid
// intersections
BOOST_AUTO_TEST_CASE(DetectorNavigatorTestsMultipleIntersection) {
  // Construct a cubic detector with a cylindrical surface
  auto bounds = std::make_unique<Acts::CuboidVolumeBounds>(10, 10, 10);
  auto surface = Acts::Surface::makeShared<Acts::CylinderSurface>(
      Acts::Transform3::Identity(),
      std::make_shared<Acts::CylinderBounds>(4, 9));
  auto volume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), geoContext,
      "volume", Acts::Transform3::Identity(), std::move(bounds), {surface}, {},
      Acts::Experimental::tryNoVolumes(),
      Acts::Experimental::tryAllPortalsAndSurfaces());

  volume->assignGeometryId(1);
  surface->assignGeometryId(2);
  int id = 3;
  for (auto& port : volume->portalPtrs()) {
    port->surface().assignGeometryId(id);
    id++;
  }

  auto detector = Acts::Experimental::Detector::makeShared(
      "detector", {volume}, Acts::Experimental::tryRootVolumes());

  using Stepper = Acts::StraightLineStepper;
  using Navigator = Acts::Experimental::DetectorNavigator;
  using Propagator = Acts::Propagator<Stepper, Navigator>;
  using ActorList = Acts::ActorList<StateRecorder, Acts::EndOfWorldReached>;
  using PropagatorOptions = Propagator::Options<ActorList>;

  Navigator::Config navCfg;
  navCfg.detector = detector.get();

  Stepper stepper;

  Navigator navigator(navCfg,
                      Acts::getDefaultLogger("DetectorNavigator",
                                             Acts::Logging::Level::VERBOSE));

  PropagatorOptions options(geoContext, mfContext);
  options.direction = Acts::Direction::Forward();

  Propagator propagator(
      stepper, navigator,
      Acts::getDefaultLogger("Propagator", Acts::Logging::Level::VERBOSE));

  // Forward and backward propagation
  // should be consistent between each other
  // and the cylindrical surface should be
  // reached in two points during navigation
  Acts::Vector4 posFwd(-5, 0, 0, 0);
  Acts::CurvilinearTrackParameters startFwd(
      posFwd, 0_degree, 90_degree, 1_e / 1_GeV, std::nullopt,
      Acts::ParticleHypothesis::electron());

  auto resultFwd = propagator.propagate(startFwd, options).value();
  auto statesFwd = resultFwd.get<StateRecorder::result_type>();

  options.direction = Acts::Direction::Backward();
  Acts::Vector4 posBwd(5, 0, 0, 0);
  Acts::CurvilinearTrackParameters startBwd(
      posBwd, 0_degree, 90_degree, 1_e / 1_GeV, std::nullopt,
      Acts::ParticleHypothesis::electron());

  auto resultBwd = propagator.propagate(startBwd, options).value();
  auto statesBwd = resultBwd.get<StateRecorder::result_type>();

  // 4 steps to reach the end of world
  // + 1 recording in the post-step
  // + 1 recording before the stepping loop
  BOOST_CHECK_EQUAL(statesFwd.size(), 5u);
  BOOST_CHECK_EQUAL(statesFwd.size(), statesBwd.size());
  BOOST_CHECK_EQUAL(statesFwd[0].surfaceCandidates.size(), 3u);
  BOOST_CHECK_EQUAL(statesBwd[0].surfaceCandidates.size(), 3u);

  // Action list call before the first step
  // Starting in the volume
  BOOST_CHECK_EQUAL(statesFwd[0].currentVolume->geometryId(), 1);
  BOOST_CHECK_EQUAL(statesFwd[0].currentSurface, nullptr);
  BOOST_CHECK_EQUAL(statesFwd[0].currentPortal, nullptr);

  // First intersection of the cylindrical surface
  BOOST_CHECK_EQUAL(statesFwd[1].currentVolume->geometryId(), 1);
  BOOST_CHECK_EQUAL(statesFwd[1].currentSurface->geometryId(), 2);
  BOOST_CHECK_EQUAL(statesFwd[1].currentPortal, nullptr);
  CHECK_CLOSE_REL(statesFwd[1].position.x(), -4, 1e-6);

  // Second intersection of the cylindrical surface
  BOOST_CHECK_EQUAL(statesFwd[2].currentVolume->geometryId(), 1);
  BOOST_CHECK_EQUAL(statesFwd[2].currentSurface->geometryId(), 2);
  BOOST_CHECK_EQUAL(statesFwd[2].currentPortal, nullptr);
  CHECK_CLOSE_REL(statesFwd[2].position.x(), 4, 1e-6);

  // Step to the volume|endOfWorld boundary
  BOOST_CHECK_EQUAL(statesFwd[3].currentVolume, nullptr);
  BOOST_CHECK_EQUAL(statesFwd[3].currentSurface,
                    &(statesFwd[3].currentPortal->surface()));
  BOOST_CHECK_EQUAL(statesFwd[3].currentPortal->surface().geometryId(), 6);

  // Step to the end of world
  BOOST_CHECK(navigator.endOfWorldReached(statesFwd[4]));
  BOOST_CHECK(navigator.endOfWorldReached(statesBwd[4]));

  // Step to the endOfWorld|volume boundary
  BOOST_CHECK_EQUAL(statesBwd[3].currentVolume, nullptr);
  BOOST_CHECK_EQUAL(statesBwd[3].currentSurface,
                    &(statesBwd[3].currentPortal->surface()));
  BOOST_CHECK_EQUAL(statesBwd[3].currentPortal->surface().geometryId(), 5);

  // Second intersection of the cylindrical surface
  BOOST_CHECK_EQUAL(statesBwd[2].currentVolume->geometryId(), 1);
  BOOST_CHECK_EQUAL(statesBwd[2].currentSurface->geometryId(), 2);
  BOOST_CHECK_EQUAL(statesBwd[2].currentPortal, nullptr);
  CHECK_CLOSE_REL(statesBwd[2].position.x(), -4, 1e-6);

  // First intersection of the cylindrical surface
  BOOST_CHECK_EQUAL(statesBwd[1].currentVolume->geometryId(), 1);
  BOOST_CHECK_EQUAL(statesBwd[1].currentSurface->geometryId(), 2);
  BOOST_CHECK_EQUAL(statesBwd[1].currentPortal, nullptr);
  CHECK_CLOSE_REL(statesBwd[1].position.x(), 4, 1e-6);

  // Action list call before the first step
  // Starting in the volume
  BOOST_CHECK_EQUAL(statesBwd[0].currentVolume->geometryId(), 1);
  BOOST_CHECK_EQUAL(statesBwd[0].currentSurface, nullptr);
  BOOST_CHECK_EQUAL(statesBwd[0].currentPortal, nullptr);
}

BOOST_AUTO_TEST_SUITE_END()
