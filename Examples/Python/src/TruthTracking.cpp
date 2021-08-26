// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/TruthTracking/ParticleSelector.hpp"
#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"
#include "ActsExamples/TruthTracking/TrackSelector.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"
#include "ActsExamples/TruthTracking/TruthVertexFinder.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addTruthTracking(Context& ctx) {
  auto mex = ctx.get("examples");

  {
    using Alg = ActsExamples::TruthTrackFinder;
    using Config = Alg::Config;

    auto alg = py::class_<Alg, BareAlgorithm, std::shared_ptr<Alg>>(
                   mex, "TruthTrackFinder")
                   .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(inputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(outputProtoTracks);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Alg = ActsExamples::TruthSeedSelector;
    using Config = Alg::Config;

    auto alg = py::class_<Alg, BareAlgorithm, std::shared_ptr<Alg>>(
                   mex, "TruthSeedSelector")
                   .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(inputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(outputParticles);
    ACTS_PYTHON_MEMBER(rhoMin);
    ACTS_PYTHON_MEMBER(rhoMax);
    ACTS_PYTHON_MEMBER(zMin);
    ACTS_PYTHON_MEMBER(zMax);
    ACTS_PYTHON_MEMBER(phiMin);
    ACTS_PYTHON_MEMBER(phiMax);
    ACTS_PYTHON_MEMBER(etaMin);
    ACTS_PYTHON_MEMBER(etaMax);
    ACTS_PYTHON_MEMBER(absEtaMin);
    ACTS_PYTHON_MEMBER(absEtaMax);
    ACTS_PYTHON_MEMBER(ptMin);
    ACTS_PYTHON_MEMBER(ptMax);
    ACTS_PYTHON_MEMBER(keepNeutral);
    ACTS_PYTHON_MEMBER(nHitsMin);
    ACTS_PYTHON_MEMBER(nHitsMax);
    ACTS_PYTHON_STRUCT_END();

    pythonRangeProperty(c, "rho", &Config::rhoMin, &Config::rhoMax);
    pythonRangeProperty(c, "z", &Config::zMin, &Config::zMax);
    pythonRangeProperty(c, "phi", &Config::phiMin, &Config::phiMax);
    pythonRangeProperty(c, "eta", &Config::etaMin, &Config::etaMax);
    pythonRangeProperty(c, "absEta", &Config::absEtaMin, &Config::absEtaMax);
    pythonRangeProperty(c, "pt", &Config::ptMin, &Config::ptMax);
    pythonRangeProperty(c, "nHits", &Config::nHitsMin, &Config::nHitsMax);
  }

  {
    using Alg = ActsExamples::ParticleSmearing;
    using Config = Alg::Config;

    auto alg = py::class_<Alg, BareAlgorithm, std::shared_ptr<Alg>>(
                   mex, "ParticleSmearing")
                   .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(outputTrackParameters);
    ACTS_PYTHON_MEMBER(sigmaD0);
    ACTS_PYTHON_MEMBER(sigmaD0PtA);
    ACTS_PYTHON_MEMBER(sigmaD0PtB);
    ACTS_PYTHON_MEMBER(sigmaZ0);
    ACTS_PYTHON_MEMBER(sigmaZ0PtA);
    ACTS_PYTHON_MEMBER(sigmaZ0PtB);
    ACTS_PYTHON_MEMBER(sigmaT0);
    ACTS_PYTHON_MEMBER(sigmaPhi);
    ACTS_PYTHON_MEMBER(sigmaTheta);
    ACTS_PYTHON_MEMBER(sigmaPRel);
    ACTS_PYTHON_MEMBER(initialVarInflation);
    ACTS_PYTHON_MEMBER(randomNumbers);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Alg = ActsExamples::ParticleSelector;
    using Config = Alg::Config;

    auto alg = py::class_<Alg, BareAlgorithm, std::shared_ptr<Alg>>(
                   mex, "ParticleSelector")
                   .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(outputParticles);
    ACTS_PYTHON_MEMBER(rhoMin);
    ACTS_PYTHON_MEMBER(rhoMax);
    ACTS_PYTHON_MEMBER(absZMin);
    ACTS_PYTHON_MEMBER(absZMax);
    ACTS_PYTHON_MEMBER(timeMin);
    ACTS_PYTHON_MEMBER(timeMax);
    ACTS_PYTHON_MEMBER(phiMin);
    ACTS_PYTHON_MEMBER(phiMax);
    ACTS_PYTHON_MEMBER(etaMin);
    ACTS_PYTHON_MEMBER(etaMax);
    ACTS_PYTHON_MEMBER(absEtaMin);
    ACTS_PYTHON_MEMBER(absEtaMax);
    ACTS_PYTHON_MEMBER(ptMin);
    ACTS_PYTHON_MEMBER(ptMax);
    ACTS_PYTHON_MEMBER(removeCharged);
    ACTS_PYTHON_MEMBER(removeNeutral);
    ACTS_PYTHON_STRUCT_END();

    pythonRangeProperty(c, "rho", &Config::rhoMin, &Config::rhoMax);
    pythonRangeProperty(c, "absZ", &Config::absZMin, &Config::absZMax);
    pythonRangeProperty(c, "time", &Config::timeMin, &Config::timeMax);
    pythonRangeProperty(c, "phi", &Config::phiMin, &Config::phiMax);
    pythonRangeProperty(c, "eta", &Config::etaMin, &Config::etaMax);
    pythonRangeProperty(c, "absEta", &Config::absEtaMin, &Config::absEtaMax);
    pythonRangeProperty(c, "pt", &Config::ptMin, &Config::ptMax);
  }

  {
    using Alg = ActsExamples::TrackSelector;
    using Config = Alg::Config;

    auto alg = py::class_<Alg, BareAlgorithm, std::shared_ptr<Alg>>(
                   mex, "TrackSelector")
                   .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputTrackParameters);
    ACTS_PYTHON_MEMBER(outputTrackParameters);
    ACTS_PYTHON_MEMBER(outputTrackIndices);
    ACTS_PYTHON_MEMBER(loc0Min);
    ACTS_PYTHON_MEMBER(loc0Max);
    ACTS_PYTHON_MEMBER(loc1Min);
    ACTS_PYTHON_MEMBER(loc1Max);
    ACTS_PYTHON_MEMBER(timeMin);
    ACTS_PYTHON_MEMBER(timeMax);
    ACTS_PYTHON_MEMBER(phiMin);
    ACTS_PYTHON_MEMBER(phiMax);
    ACTS_PYTHON_MEMBER(etaMin);
    ACTS_PYTHON_MEMBER(etaMax);
    ACTS_PYTHON_MEMBER(absEtaMin);
    ACTS_PYTHON_MEMBER(absEtaMax);
    ACTS_PYTHON_MEMBER(ptMin);
    ACTS_PYTHON_MEMBER(ptMax);
    ACTS_PYTHON_MEMBER(removeCharged);
    ACTS_PYTHON_MEMBER(removeNeutral);
    ACTS_PYTHON_STRUCT_END();

    pythonRangeProperty(c, "loc0", &Config::loc0Min, &Config::loc0Max);
    pythonRangeProperty(c, "loc1", &Config::loc1Min, &Config::loc1Max);
    pythonRangeProperty(c, "time", &Config::timeMin, &Config::timeMax);
    pythonRangeProperty(c, "phi", &Config::phiMin, &Config::phiMax);
    pythonRangeProperty(c, "eta", &Config::etaMin, &Config::etaMax);
    pythonRangeProperty(c, "absEta", &Config::absEtaMin, &Config::absEtaMax);
    pythonRangeProperty(c, "pt", &Config::ptMin, &Config::ptMax);
  }

  {
    using Alg = ActsExamples::TruthVertexFinder;
    using Config = Alg::Config;

    auto alg = py::class_<Alg, BareAlgorithm, std::shared_ptr<Alg>>(
                   mex, "TruthVertexFinder")
                   .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(outputProtoVertices);
    ACTS_PYTHON_MEMBER(excludeSecondaries);
    ACTS_PYTHON_MEMBER(separateSecondaries);
    ACTS_PYTHON_STRUCT_END();
  }
}

}  // namespace Acts::Python