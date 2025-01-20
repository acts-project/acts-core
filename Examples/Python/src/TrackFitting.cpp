// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/TrackFitting/BetheHeitlerApprox.hpp"
#include "Acts/TrackFitting/GsfOptions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/EventData/ScalingCalibrator.hpp"
#include "ActsExamples/TrackFitting/RefittingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFitterFunction.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"

#include <cstddef>
#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;
using namespace py::literals;

namespace Acts::Python {

void addTrackFitting(Context& ctx) {
  auto mex = ctx.get("examples");

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::TrackFittingAlgorithm, mex, "TrackFittingAlgorithm",
      inputMeasurements, inputProtoTracks, inputInitialTrackParameters,
      inputClusters, outputTracks, fit, pickTrack, calibrator);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::RefittingAlgorithm, mex,
                                "RefittingAlgorithm", inputTracks, outputTracks,
                                fit, pickTrack, initialVarInflation);

  {
    py::class_<TrackFitterFunction, std::shared_ptr<TrackFitterFunction>>(
        mex, "TrackFitterFunction");

    mex.def(
        "makeKalmanFitterFunction",
        [](std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
           std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
           bool multipleScattering, bool energyLoss,
           double reverseFilteringMomThreshold,
           double reverseFilteringCovarianceScaling,
           Acts::FreeToBoundCorrection freeToBoundCorrection,
           Logging::Level level) {
          return ActsExamples::makeKalmanFitterFunction(
              trackingGeometry, magneticField, multipleScattering, energyLoss,
              reverseFilteringMomThreshold, reverseFilteringCovarianceScaling,
              freeToBoundCorrection, *Acts::getDefaultLogger("Kalman", level));
        },
        "trackingGeometry"_a, "magneticField"_a, "multipleScattering"_a,
        "energyLoss"_a, "reverseFilteringMomThreshold"_a,
        "reverseFilteringCovarianceScaling"_a, "freeToBoundCorrection"_a,
        "level"_a);

    py::class_<MeasurementCalibrator, std::shared_ptr<MeasurementCalibrator>>(
        mex, "MeasurementCalibrator");

    mex.def("makePassThroughCalibrator",
            []() -> std::shared_ptr<MeasurementCalibrator> {
              return std::make_shared<PassThroughCalibrator>();
            });

    mex.def(
        "makeScalingCalibrator",
        [](const char* path) -> std::shared_ptr<MeasurementCalibrator> {
          return std::make_shared<ActsExamples::ScalingCalibrator>(path);
        },
        py::arg("path"));

    py::enum_<Acts::ComponentMergeMethod>(mex, "ComponentMergeMethod")
        .value("mean", Acts::ComponentMergeMethod::eMean)
        .value("maxWeight", Acts::ComponentMergeMethod::eMaxWeight);

    py::enum_<ActsExamples::MixtureReductionAlgorithm>(
        mex, "MixtureReductionAlgorithm")
        .value("weightCut", MixtureReductionAlgorithm::weightCut)
        .value("KLDistance", MixtureReductionAlgorithm::KLDistance);

    py::class_<ActsExamples::BetheHeitlerApprox>(mex, "AtlasBetheHeitlerApprox")
        .def_static(
            "loadFromFiles", &ActsExamples::BetheHeitlerApprox::loadFromFiles,
            "lowParametersPath"_a, "highParametersPath"_a, "lowLimit"_a = 0.1,
            "highLimit"_a = 0.2, "clampToRange"_a = false)
        .def_static(
            "makeDefault",
            [](bool clampToRange) {
              return Acts::makeDefaultBetheHeitlerApprox(clampToRange);
            },
            "clampToRange"_a = false);

    mex.def(
        "makeGsfFitterFunction",
        [](std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
           std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
           BetheHeitlerApprox betheHeitlerApprox, std::size_t maxComponents,
           double weightCutoff, Acts::ComponentMergeMethod componentMergeMethod,
           ActsExamples::MixtureReductionAlgorithm mixtureReductionAlgorithm,
           double reverseFilteringCovarianceScaling, Logging::Level level) {
          return ActsExamples::makeGsfFitterFunction(
              trackingGeometry, magneticField, betheHeitlerApprox,
              maxComponents, weightCutoff, componentMergeMethod,
              mixtureReductionAlgorithm, reverseFilteringCovarianceScaling,
              *Acts::getDefaultLogger("GSFFunc", level));
        },
        "trackingGeometry"_a, "magneticField"_a, "betheHeitlerApprox"_a,
        "maxComponents"_a, "weightCutoff"_a, "componentMergeMethod"_a,
        "mixtureReductionAlgorithm"_a, "reverseFilteringCovarianceScaling"_a,
        "level"_a);

    mex.def(
        "makeGlobalChiSquareFitterFunction",
        [](std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
           std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
           bool multipleScattering, bool energyLoss,
           Acts::FreeToBoundCorrection freeToBoundCorrection,
           std::size_t nUpdateMax, double relChi2changeCutOff,
           Logging::Level level) {
          return ActsExamples::makeGlobalChiSquareFitterFunction(
              trackingGeometry, magneticField, multipleScattering, energyLoss,
              freeToBoundCorrection, nUpdateMax, relChi2changeCutOff,
              *Acts::getDefaultLogger("Gx2f", level));
        },
        py::arg("trackingGeometry"), py::arg("magneticField"),
        py::arg("multipleScattering"), py::arg("energyLoss"),
        py::arg("freeToBoundCorrection"), py::arg("nUpdateMax"),
        py::arg("relChi2changeCutOff"), py::arg("level"));
  }

  {
    py::class_<FreeToBoundCorrection>(mex, "FreeToBoundCorrection")
        .def(py::init<>())
        .def(py::init<bool>(), py::arg("apply") = false)
        .def(py::init<bool, double, double>(), py::arg("apply") = false,
             py::arg("alpha") = 0.1, py::arg("beta") = 2);
  }
}

}  // namespace Acts::Python
