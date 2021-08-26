// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/HepMC/HepMCProcessExtractor.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Reader.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Writer.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;

namespace Acts::Python {
void addHepMC3(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  auto hepmc3 = mex.def_submodule("_hepmc3");

  {
    using Alg = ActsExamples::HepMCProcessExtractor;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            hepmc3, "HepMCProcessExtractor")
            .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    auto c = py::class_<Alg::Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Alg::Config);
    ACTS_PYTHON_MEMBER(inputEvents);
    ACTS_PYTHON_MEMBER(outputSimulationProcesses);
    ACTS_PYTHON_MEMBER(extractionProcess);
    ACTS_PYTHON_MEMBER(absPdgMin);
    ACTS_PYTHON_MEMBER(absPdgMax);
    ACTS_PYTHON_MEMBER(pMin);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::HepMC3AsciiWriter;

    auto alg =
        py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
            hepmc3, "HepMC3AsciiWriter")
            .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_MEMBER(outputStem);
    ACTS_PYTHON_MEMBER(inputEvents);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Reader = ActsExamples::HepMC3AsciiReader;

    auto alg =
        py::class_<Reader, ActsExamples::IReader, std::shared_ptr<Reader>>(
            hepmc3, "HepMC3AsciiReader")
            .def(py::init<const Reader::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    auto c = py::class_<Reader::Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Reader::Config);
    ACTS_PYTHON_MEMBER(inputDir);
    ACTS_PYTHON_MEMBER(inputStem);
    ACTS_PYTHON_MEMBER(outputEvents);
    ACTS_PYTHON_STRUCT_END();
  }
}
}  // namespace Acts::Python