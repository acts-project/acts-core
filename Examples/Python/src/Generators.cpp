// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Generators/EventGenerator.hpp"
#include "ActsExamples/Generators/MultiplicityGenerators.hpp"
#include "ActsExamples/Generators/ParametricParticleGenerator.hpp"
#include "ActsExamples/Generators/VertexGenerators.hpp"
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"

#include <cmath>
#include <memory>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace {
double thetaToEta(double theta) {
  return -1 * std::log(std::tan(theta / 2.));
}
double etaToTheta(double eta) {
  return 2 * std::atan(std::exp(-eta));
}
}  // namespace

namespace Acts::Python {

void addGenerators(Context& ctx) {
  auto mex = ctx.get("examples");
  {
    using Config = ActsExamples::EventGenerator::Config;
    auto gen = py::class_<ActsExamples::EventGenerator, ActsExamples::IReader,
                          std::shared_ptr<ActsExamples::EventGenerator>>(
                   mex, "EventGenerator")
                   .def(py::init<const Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly(
                       "config", &ActsExamples::EventGenerator::config);

    py::class_<ActsExamples::EventGenerator::VertexGenerator,
               std::shared_ptr<ActsExamples::EventGenerator::VertexGenerator>>(
        gen, "VertexGenerator");
    py::class_<
        ActsExamples::EventGenerator::ParticlesGenerator,
        std::shared_ptr<ActsExamples::EventGenerator::ParticlesGenerator>>(
        gen, "ParticlesGenerator");
    py::class_<
        ActsExamples::EventGenerator::MultiplicityGenerator,
        std::shared_ptr<ActsExamples::EventGenerator::MultiplicityGenerator>>(
        gen, "MultiplicityGenerator");

    using EventGenerator = ActsExamples::EventGenerator;
    using Generator = EventGenerator::Generator;
    py::class_<Generator>(gen, "Generator")
        .def(py::init<>())
        .def(py::init<std::shared_ptr<EventGenerator::MultiplicityGenerator>,
                      std::shared_ptr<EventGenerator::VertexGenerator>,
                      std::shared_ptr<EventGenerator::ParticlesGenerator>>(),
             py::arg("multiplicity"), py::arg("vertex"), py::arg("particles"))
        .def_readwrite("multiplicity", &Generator::multiplicity)
        .def_readwrite("vertex", &Generator::vertex)
        .def_readwrite("particles", &Generator::particles);

    py::class_<Config>(gen, "Config")
        .def(py::init<>())
        .def_readwrite("outputParticles", &Config::outputParticles)
        .def_readwrite("generators", &Config::generators)
        .def_readwrite("randomNumbers", &Config::randomNumbers);
  }

  py::class_<ActsExamples::GaussianVertexGenerator,
             ActsExamples::EventGenerator::VertexGenerator,
             std::shared_ptr<ActsExamples::GaussianVertexGenerator>>(
      mex, "GaussianVertexGenerator")
      .def(py::init<>())
      .def(py::init([](Acts::Vector4 stddev, Acts::Vector4 mean) {
             ActsExamples::GaussianVertexGenerator g;
             g.stddev = stddev;
             g.mean = mean;
             return g;
           }),
           py::arg("stddev"), py::arg("mean"))
      .def_readwrite("stddev", &ActsExamples::GaussianVertexGenerator::stddev)
      .def_readwrite("mean", &ActsExamples::GaussianVertexGenerator::mean);

  py::class_<ActsExamples::FixedVertexGenerator,
             ActsExamples::EventGenerator::VertexGenerator,
             std::shared_ptr<ActsExamples::FixedVertexGenerator>>(
      mex, "FixedVertexGenerator")
      .def(py::init<>())
      .def(py::init([](Acts::Vector4 v) {
             ActsExamples::FixedVertexGenerator g;
             g.fixed = v;
             return g;
           }),
           py::arg("fixed"))
      .def_readwrite("fixed", &ActsExamples::FixedVertexGenerator::fixed);

  py::class_<ActsExamples::SimParticle>(mex, "SimParticle");
  py::class_<ActsExamples::SimParticleContainer>(mex, "SimParticleContainer");

  {
    using Config = ActsExamples::ParametricParticleGenerator::Config;
    auto gen =
        py::class_<ActsExamples::ParametricParticleGenerator,
                   ActsExamples::EventGenerator::ParticlesGenerator,
                   std::shared_ptr<ActsExamples::ParametricParticleGenerator>>(
            mex, "ParametricParticleGenerator")
            .def(py::init<const Config&>());

    py::class_<Config>(gen, "Config")
        .def(py::init<>())
        .def_readwrite("phiMin", &Config::phiMin)
        .def_readwrite("phiMax", &Config::phiMax)
        .def_readwrite("thetaMin", &Config::thetaMin)
        .def_readwrite("thetaMax", &Config::thetaMax)
        .def_readwrite("pMin", &Config::pMin)
        .def_readwrite("pMax", &Config::pMax)
        .def_readwrite("pTransverse", &Config::pTransverse)
        .def_readwrite("pdg", &Config::pdg)
        .def_readwrite("randomizeCharge", &Config::randomizeCharge)
        .def_readwrite("numParticles", &Config::numParticles)
        .def_property(
            "p",
            [](Config& cfg) {
              return std::pair{cfg.pMin, cfg.pMax};
            },
            [](Config& cfg, std::pair<double, double> value) {
              cfg.pMin = value.first;
              cfg.pMax = value.second;
            })
        .def_property(
            "phi",
            [](Config& cfg) {
              return std::pair{cfg.phiMin, cfg.phiMax};
            },
            [](Config& cfg, std::pair<double, double> value) {
              cfg.phiMin = value.first;
              cfg.phiMax = value.second;
            })
        .def_property(
            "theta",
            [](Config& cfg) {
              return std::pair{cfg.thetaMin, cfg.thetaMax};
            },
            [](Config& cfg, std::pair<double, double> value) {
              cfg.thetaMin = value.first;
              cfg.thetaMax = value.second;
            })
        .def_property(
            "eta",
            [](Config& cfg) {
              return std::pair{thetaToEta(cfg.thetaMin),
                               thetaToEta(cfg.thetaMax)};
            },
            [](Config& cfg, std::pair<double, double> value) {
              cfg.thetaMin = etaToTheta(value.first);
              cfg.thetaMax = etaToTheta(value.second);
            });
  }

  py::class_<ActsExamples::FixedMultiplicityGenerator,
             ActsExamples::EventGenerator::MultiplicityGenerator,
             std::shared_ptr<ActsExamples::FixedMultiplicityGenerator>>(
      mex, "FixedMultiplicityGenerator")
      .def(py::init<>())
      .def(py::init([](size_t n) {
             ActsExamples::FixedMultiplicityGenerator g;
             g.n = n;
             return g;
           }),
           py::arg("n"))
      .def_readwrite("n", &ActsExamples::FixedMultiplicityGenerator::n);

  py::class_<ActsExamples::PoissonMultiplicityGenerator,
             ActsExamples::EventGenerator::MultiplicityGenerator,
             std::shared_ptr<ActsExamples::PoissonMultiplicityGenerator>>(
      mex, "PoissonMultiplicityGenerator")
      .def(py::init<>())
      .def(py::init([](double mean) {
             ActsExamples::PoissonMultiplicityGenerator g;
             g.mean = mean;
             return g;
           }),
           py::arg("mean"))
      .def_readwrite("mean", &ActsExamples::PoissonMultiplicityGenerator::mean);
}
}  // namespace Acts::Python
