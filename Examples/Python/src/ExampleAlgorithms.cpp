// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "ActsExamples/Fatras/FatrasAlgorithm.hpp"
#include "ActsExamples/Io/Json/JsonGeometryList.hpp"
#include "ActsExamples/Printers/HitsPrinter.hpp"
#include "ActsExamples/Printers/ParticlesPrinter.hpp"
#include "ActsExamples/Printers/TrackParametersPrinter.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addExampleAlgorithms(Context& ctx) {
  auto mex = ctx.get("examples");

  mex.def("readJsonGeometryList", ActsExamples::readJsonGeometryList);

  {
    using Config = ActsExamples::FatrasAlgorithm::Config;

    auto alg =
        py::class_<ActsExamples::FatrasAlgorithm, ActsExamples::BareAlgorithm,
                   std::shared_ptr<ActsExamples::FatrasAlgorithm>>(
            mex, "FatrasAlgorithm")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config",
                                   &ActsExamples::FatrasAlgorithm::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(outputParticlesInitial);
    ACTS_PYTHON_MEMBER(outputParticlesFinal);
    ACTS_PYTHON_MEMBER(outputSimHits);
    ACTS_PYTHON_MEMBER(imputParametrisationNuclearInteraction);
    ACTS_PYTHON_MEMBER(randomNumbers);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(magneticField);
    ACTS_PYTHON_MEMBER(pMin);
    ACTS_PYTHON_MEMBER(emScattering);
    ACTS_PYTHON_MEMBER(emEnergyLossIonisation);
    ACTS_PYTHON_MEMBER(emEnergyLossRadiation);
    ACTS_PYTHON_MEMBER(emPhotonConversion);
    ACTS_PYTHON_MEMBER(generateHitsOnSensitive);
    ACTS_PYTHON_MEMBER(generateHitsOnMaterial);
    ACTS_PYTHON_MEMBER(generateHitsOnPassive);
    ACTS_PYTHON_MEMBER(averageHitsPerParticle);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Alg = ActsExamples::ParticlesPrinter;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            mex, "ParticlesPrinter")
            .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    py::class_<Alg::Config>(alg, "Config")
        .def(py::init<>())
        .def_readwrite("inputParticles", &Alg::Config::inputParticles);
  }

  {
    using Alg = ActsExamples::HitsPrinter;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            mex, "HitsPrinter")
            .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    auto c = py::class_<Alg::Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Alg::Config);
    ACTS_PYTHON_MEMBER(inputClusters);
    ACTS_PYTHON_MEMBER(inputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(inputHitIds);
    ACTS_PYTHON_MEMBER(selectIndexStart);
    ACTS_PYTHON_MEMBER(selectIndexLength);
    ACTS_PYTHON_MEMBER(selectVolume);
    ACTS_PYTHON_MEMBER(selectLayer);
    ACTS_PYTHON_MEMBER(selectModule);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Alg = ActsExamples::TrackParametersPrinter;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            mex, "TrackParametersPrinter")
            .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    py::class_<Alg::Config>(alg, "Config")
        .def(py::init<>())
        .def_readwrite("inputTrackParameters",
                       &Alg::Config::inputTrackParameters);
  }
}

}  // namespace Acts::Python