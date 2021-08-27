// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/TrialAndErrorNavigator.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "ActsExamples/GenericDetector/GenericDetector.hpp"

#include <fstream>

#include "PropagationDatasets.hpp"

namespace ds = ActsTests::PropagationDatasets;

using Stepper = Acts::EigenStepper<>;

const Acts::GeometryContext gctx;
const Acts::MagneticFieldContext mctx;
const auto logger = Acts::getDefaultLogger("Logger", Acts::Logging::INFO);

enum class NavigatorType { standard, trialAndError };

const std::string datestring = []() {
  char text[100];
  time_t now = std::time(nullptr);
  std::tm* t = std::localtime(&now);
  std::strftime(text, sizeof(text), "-%Y%m%d-%H%M%S", t);
  return std::string(text);
}();

auto export_track_to_csv(const std::vector<Acts::detail::Step>& steps,
                         NavigatorType type) {
  static int count_standard = 0;
  static int count_trialerror = 0;

  auto file = [&]() {
    if (type == NavigatorType::standard) {
      const auto filename = std::to_string(count_standard++) +
                            std::string("-standard") + datestring + ".csv";
      return std::ofstream(filename, std::ios_base::trunc);
    } else {
      const auto filename = std::to_string(count_trialerror++) +
                            std::string("-trialerror") + datestring + ".csv";
      return std::ofstream(filename, std::ios_base::trunc);
    }
  }();

  file << "surface,x,y,z\n";

  for (const auto& step : steps) {
    file << step.surface->geometryId().value() << "," << step.position[0] << ","
         << step.position[1] << "," << step.position[2] << "\n";
  }
}

void print_difference(const std::vector<Acts::detail::Step>& standard,
                      const std::vector<Acts::detail::Step>& trialAndError) {
  std::vector<Acts::GeometryIdentifier> standardIds;
  for (const auto& step : standard) {
    standardIds.push_back(step.surface->geometryId());
  }

  std::vector<Acts::GeometryIdentifier> trialErrorFound(
      standardIds.size(), Acts::GeometryIdentifier{});
  std::vector<Acts::GeometryIdentifier> trialErrorNotFound;

  for (const auto& step : trialAndError) {
    const auto found =
        std::find(standardIds.begin(), standardIds.end(), step.surface->geometryId());

    if (found == standardIds.end()) {
      trialErrorNotFound.push_back(step.surface->geometryId());
    } else {
      const auto i = std::distance(standardIds.begin(), found);
      trialErrorFound[i] = step.surface->geometryId();
    }
  }

  std::cout << "Standard\tTrial&Error\n";
  for (auto i = 0ul; i < standardIds.size(); ++i) {
    std::cout << standardIds[i] << "\t" << trialErrorFound[i] << "\n";
  }

  std::cout << "\nNot found: ";
  for (const auto& geoid : trialErrorNotFound) {
    std::cout << geoid << ",  ";
  }
  std::cout << std::endl;
}

auto getGenericTrackingGeometry() {
  static std::shared_ptr<const Acts::TrackingGeometry> geometry;
  static std::vector<const Acts::Surface*> surfaces;

  if (!geometry) {
    static GenericDetector detector;
    geometry = std::shared_ptr(detector.finalize({}, {}).first);

    geometry->visitSurfaces([&](auto surface) { surfaces.push_back(surface); });
  }

  return std::tie(geometry, surfaces);
}

void compareNavigation(
    const std::shared_ptr<const Acts::TrackingGeometry>& trackingGeometry,
    const std::shared_ptr<Acts::MagneticFieldProvider>& magneticField,
    const Acts::BoundTrackParameters& startPars,
    const std::vector<const Acts::Surface*>& allSurfaces) {
  // Do normal navigation
  const auto normalStepLog = [&]() {
    Acts::Navigator::Config cfg;
    cfg.trackingGeometry = trackingGeometry;
    Acts::Navigator navigator(cfg);
    using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
    Propagator prop(Stepper(magneticField), navigator);

    using Actors = Acts::ActionList<Acts::detail::SteppingLogger>;
    using Aborters = Acts::AbortList<Acts::EndOfWorldReached>;

    Acts::PropagatorOptions<Actors, Aborters> propOptions(
        gctx, mctx, Acts::LoggerWrapper(*logger));

    auto res = prop.propagate(startPars, propOptions);

    auto stepLog =
        (*res).get<Acts::detail::SteppingLogger::result_type>().steps;

    stepLog.erase(
        std::remove_if(stepLog.begin(), stepLog.end(),
                       [](const auto& s) { return s.surface == nullptr; }),
        stepLog.end());

    return stepLog;
  }();

  // Trial-and-error navigation
  const auto trialAndErrorStepLog = [&]() {
    using Navigator = Acts::TrialAndErrorNavigator;
    Navigator navigator;
    navigator.guider.possibleSurfaces = allSurfaces;

    using Propagator = Acts::Propagator<Stepper, Navigator>;
    Propagator prop(Stepper(magneticField), navigator);

    using Actors = Acts::ActionList<Acts::TrialAndErrorNavigatorInitializer,
                                    Acts::detail::SteppingLogger>;
    using Aborters = Acts::AbortList<>;

    Acts::PropagatorOptions<Actors, Aborters> propOptions(
        gctx, mctx, Acts::LoggerWrapper(*logger));

    propOptions.actionList.get<Acts::TrialAndErrorNavigatorInitializer>()
        .surfaces = &navigator.guider.possibleSurfaces;

    auto res = prop.propagate(startPars, propOptions);

    auto stepLog =
        (*res).get<Acts::detail::SteppingLogger::result_type>().steps;

    stepLog.erase(
        std::remove_if(stepLog.begin(), stepLog.end(),
                       [](const auto& s) { return s.surface == nullptr; }),
        stepLog.end());

    return stepLog;
  }();

  if (trialAndErrorStepLog.empty() && normalStepLog.empty()) {
    return;
  }

  BOOST_REQUIRE(!trialAndErrorStepLog.empty());
  BOOST_REQUIRE(!normalStepLog.empty());

  auto normal_it = normalStepLog.cbegin();
  for (const auto& tae_step : trialAndErrorStepLog) {
    normal_it =
        std::find_if(normal_it, normalStepLog.cend(), [&](const auto& s) {
          return s.surface->geometryId() == tae_step.surface->geometryId();
        });

    if (normal_it == normalStepLog.end()) {
      export_track_to_csv(normalStepLog, NavigatorType::standard);
      export_track_to_csv(trialAndErrorStepLog, NavigatorType::trialAndError);
      print_difference(normalStepLog, trialAndErrorStepLog);
    }

    BOOST_REQUIRE(normal_it != normalStepLog.end());
  }
}

BOOST_DATA_TEST_CASE(CompareNavigationGenericDetector,
                     ds::phi* ds::thetaCentral* ds::absMomentum*
                         ds::chargeNonZero* ds::magneticField,
                     phi, theta, p, q, bz) {
  const auto& [trackingGeometry, surfaces] = getGenericTrackingGeometry();

  using namespace Acts::UnitLiterals;

  Acts::BoundVector bound;
  bound[Acts::eBoundLoc0] = 0.0;
  bound[Acts::eBoundLoc1] = 0.0;
  bound[Acts::eBoundPhi] = phi;
  bound[Acts::eBoundTheta] = theta;
  bound[Acts::eBoundQOverP] = q / p;
  bound[Acts::eBoundTime] = 0.0;

  Acts::BoundTrackParameters pars(
      trackingGeometry->getBeamline()->getSharedPtr(), bound);

  auto magneticField =
      std::make_shared<Acts::ConstantBField>(Acts::Vector3{0.0, 0.0, bz});

  compareNavigation(trackingGeometry, magneticField, pars, surfaces);
}
