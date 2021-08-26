// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Vertexing/AdaptiveMultiVertexFinderAlgorithm.hpp"
#include "ActsExamples/Vertexing/IterativeVertexFinderAlgorithm.hpp"
#include "ActsExamples/Vertexing/TutorialVertexFinderAlgorithm.hpp"
#include "ActsExamples/Vertexing/VertexFitterAlgorithm.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addVertexing(Context& ctx) {
  auto mex = ctx.get("examples");

  {
    using Alg = ActsExamples::AdaptiveMultiVertexFinderAlgorithm;
    using Config = Alg::Config;

    auto alg = py::class_<Alg, BareAlgorithm, std::shared_ptr<Alg>>(
                   mex, "AdaptiveMultiVertexFinderAlgorithm")
                   .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputTrackParameters);
    ACTS_PYTHON_MEMBER(outputProtoVertices);
    ACTS_PYTHON_MEMBER(outputVertices);
    ACTS_PYTHON_MEMBER(outputTime);
    ACTS_PYTHON_MEMBER(bField);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Alg = ActsExamples::IterativeVertexFinderAlgorithm;
    using Config = Alg::Config;

    auto alg = py::class_<Alg, BareAlgorithm, std::shared_ptr<Alg>>(
                   mex, "IterativeVertexFinderAlgorithm")
                   .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputTrackParameters);
    ACTS_PYTHON_MEMBER(outputProtoVertices);
    ACTS_PYTHON_MEMBER(outputVertices);
    ACTS_PYTHON_MEMBER(outputTime);
    ACTS_PYTHON_MEMBER(bField);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Alg = ActsExamples::TutorialVertexFinderAlgorithm;
    using Config = Alg::Config;

    auto alg = py::class_<Alg, BareAlgorithm, std::shared_ptr<Alg>>(
                   mex, "TutorialVertexFinderAlgorithm")
                   .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputTrackParameters);
    ACTS_PYTHON_MEMBER(outputProtoVertices);
    ACTS_PYTHON_MEMBER(bField);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Alg = ActsExamples::VertexFitterAlgorithm;
    using Config = Alg::Config;

    auto alg = py::class_<Alg, BareAlgorithm, std::shared_ptr<Alg>>(
                   mex, "VertexFitterAlgorithm")
                   .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputTrackParameters);
    ACTS_PYTHON_MEMBER(inputProtoVertices);
    ACTS_PYTHON_MEMBER(outputVertices);
    ACTS_PYTHON_MEMBER(bField);
    ACTS_PYTHON_MEMBER(doConstrainedFit);
    ACTS_PYTHON_MEMBER(constraintPos);
    ACTS_PYTHON_MEMBER(constraintCov);
    ACTS_PYTHON_STRUCT_END();
  }
}

}  // namespace Acts::Python