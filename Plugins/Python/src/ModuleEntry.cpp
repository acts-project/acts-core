// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/ActsVersion.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;

namespace {
class PyIAlgorithm : public IAlgorithm {
 public:
  using IAlgorithm::IAlgorithm;

  std::string name() const override {
    py::gil_scoped_acquire acquire{};
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu-zero-variadic-macro-arguments"
#endif
    PYBIND11_OVERRIDE_PURE(std::string, IAlgorithm, name);
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
  }

  ProcessCode execute(const AlgorithmContext& ctx) const override {
    py::gil_scoped_acquire acquire{};
    PYBIND11_OVERRIDE_PURE(ProcessCode, IAlgorithm, execute, ctx);
  }
};

class PyBareAlgorithm : public BareAlgorithm {
 public:
  using BareAlgorithm::BareAlgorithm;

  ProcessCode execute(const AlgorithmContext& ctx) const override {
    py::gil_scoped_acquire acquire{};
    try {
      PYBIND11_OVERRIDE_PURE(ProcessCode, BareAlgorithm, execute, ctx);
    } catch (py::error_already_set& e) {
      throw;  // Error from python, handle in python.
    } catch (std::runtime_error& e) {
      throw py::type_error("Python algorithm did not conform to interface");
    }
  }
};
}  // namespace

namespace Acts::Python {
void addUnits(Context& ctx);
void addLogging(Context& ctx);
void addPdgParticle(Context& ctx);
void addAlgebra(Context& ctx);

void addPropagation(Context& ctx);

void addGeometry(Context& ctx);

void addMagneticField(Context& ctx);

void addMaterial(Context& ctx);
void addOutput(Context& ctx);
void addDetector(Context& ctx);
void addExampleAlgorithms(Context& ctx);
void addInput(Context& ctx);
void addGenerators(Context& ctx);
void addTruthTracking(Context& ctx);
void addTrackFitting(Context& ctx);
void addTrackFinding(Context& ctx);
void addVertexing(Context& ctx);

// Plugins
void addDigitization(Context& ctx);
void addPythia8(Context& ctx);
void addJson(Context& ctx);
void addHepMC3(Context& ctx);

}  // namespace Acts::Python

using namespace Acts::Python;

PYBIND11_MODULE(ActsPythonBindings, m) {
  Acts::Python::Context ctx;
  ctx.modules["main"] = &m;
  auto mex = m.def_submodule("_examples");
  ctx.modules["examples"] = &mex;
  auto prop = m.def_submodule("_propagator");
  ctx.modules["propagation"] = &prop;
  m.doc() = "Acts";

  m.attr("__version__") =
      std::tuple{Acts::VersionMajor, Acts::VersionMinor, Acts::VersionPatch};

  // {
  //   py::class_<Acts::Logger, std::unique_ptr<Acts::Logger>>(m, "Logger");
  //   m.def("getDefaultLogger", [](const std::string& name, Logging::Level lvl)
  //   {return Acts::getDefaultLogger);
  // }

  py::class_<ActsExamples::IWriter, std::shared_ptr<ActsExamples::IWriter>>(
      mex, "IWriter");

  py::class_<ActsExamples::IReader, std::shared_ptr<ActsExamples::IReader>>(
      mex, "IReader");

  py::enum_<ProcessCode>(mex, "ProcessCode")
      .value("SUCCESS", ProcessCode::SUCCESS)
      .value("ABORT", ProcessCode::ABORT)
      .value("END", ProcessCode::END);

  py::class_<WhiteBoard>(mex, "WhiteBoard")
      .def(py::init([](Acts::Logging::Level level, const std::string& name) {
             return std::make_unique<WhiteBoard>(
                 Acts::getDefaultLogger(name, level));
           }),
           py::arg("level"), py::arg("name") = "WhiteBoard")
      .def("exists", &WhiteBoard::exists);

  py::class_<Acts::GeometryContext>(m, "GeometryContext");
  //   py::class_<Acts::MagneticFieldContext>(m, "MagneticFieldContext");

  py::class_<AlgorithmContext>(mex, "AlgorithmContext")
      .def(py::init<size_t, size_t, WhiteBoard&>())
      .def_readonly("algorithmNumber", &AlgorithmContext::algorithmNumber)
      .def_readonly("eventNumber", &AlgorithmContext::eventNumber)
      .def_property_readonly("eventStore",
                             [](const AlgorithmContext& self) -> WhiteBoard& {
                               return self.eventStore;
                             })
      .def_readonly("magFieldContext", &AlgorithmContext::magFieldContext)
      .def_readonly("geoContext", &AlgorithmContext::geoContext)
      .def_readonly("calibContext", &AlgorithmContext::calibContext);

  auto iAlgorithm =
      py::class_<ActsExamples::IAlgorithm, PyIAlgorithm,
                 std::shared_ptr<ActsExamples::IAlgorithm>>(mex, "IAlgorithm")
          .def(py::init_alias<>())
          .def("execute", &IAlgorithm::execute)
          .def("name", &IAlgorithm::name);

  auto bareAlgorithm =
      py::class_<ActsExamples::BareAlgorithm,
                 std::shared_ptr<ActsExamples::BareAlgorithm>, IAlgorithm,
                 PyBareAlgorithm>(mex, "BareAlgorithm")
          .def(py::init_alias<const std::string&, Acts::Logging::Level>(),
               py::arg("name"), py::arg("level"))
          .def("execute", &BareAlgorithm::execute);

  py::class_<Acts::GeometryIdentifier>(m, "GeometryIdentifier")
      .def(py::init<>());

  using ActsExamples::Sequencer;
  using Config = Sequencer::Config;
  auto sequencer =
      py::class_<Sequencer>(mex, "Sequencer")
          .def(py::init([](Config cfg) {
            cfg.iterationCallback = []() {
              py::gil_scoped_acquire gil;
              if (PyErr_CheckSignals() != 0) {
                throw py::error_already_set{};
              }
            };
            return std::make_unique<Sequencer>(cfg);
          }))
          .def("run",
               [](Sequencer& self) {
                 py::gil_scoped_release gil;
                 int res = self.run();
                 if (res != EXIT_SUCCESS) {
                   throw std::runtime_error{"Sequencer terminated abnormally"};
                 }
               })
          .def("addContextDecorator", &Sequencer::addContextDecorator)
          .def("addAlgorithm", &Sequencer::addAlgorithm, py::keep_alive<1, 2>())
          .def("addReader", &Sequencer::addReader)
          .def("addWriter", &Sequencer::addWriter)
          .def_property_readonly("config", &Sequencer::config);

  py::class_<Config>(sequencer, "Config")
      .def(py::init<>())
      .def_readwrite("skip", &Config::skip)
      .def_readwrite("events", &Config::events)
      .def_readwrite("logLevel", &Config::logLevel)
      .def_readwrite("numThreads", &Config::numThreads)
      .def_readwrite("outputDir", &Config::outputDir);

  using ActsExamples::RandomNumbers;
  auto randomNumbers =
      py::class_<RandomNumbers, std::shared_ptr<RandomNumbers>>(mex,
                                                                "RandomNumbers")
          .def(py::init<const RandomNumbers::Config&>());

  py::class_<ActsExamples::RandomEngine>(mex, "RandomEngine").def(py::init<>());

  py::class_<RandomNumbers::Config>(randomNumbers, "Config")
      .def(py::init<>())
      .def_readwrite("seed", &RandomNumbers::Config::seed);

  addUnits(ctx);
  addLogging(ctx);
  addPdgParticle(ctx);
  addAlgebra(ctx);

  addPropagation(ctx);
  addGeometry(ctx);
  addMagneticField(ctx);
  addMaterial(ctx);
  addOutput(ctx);
  addDetector(ctx);
  addExampleAlgorithms(ctx);
  addInput(ctx);
  addGenerators(ctx);
  addTruthTracking(ctx);
  addTrackFitting(ctx);
  addTrackFinding(ctx);
  addVertexing(ctx);

  addDigitization(ctx);
  addPythia8(ctx);
  addJson(ctx);
  addHepMC3(ctx);
}