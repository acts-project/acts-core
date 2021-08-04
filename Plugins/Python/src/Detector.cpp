// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/ContextualDetector/AlignedDetector.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/GenericDetector/GenericDetector.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace ActsExamples;

namespace Acts::Python {
void addDetector(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");
  {
    py::class_<IContextDecorator, std::shared_ptr<IContextDecorator>>(
        mex, "IContextDecorator")
        .def("decorate", &IContextDecorator::decorate)
        .def("name", &IContextDecorator::name);
  }

  {
    using Config = GenericDetector::Config;

    auto gd = py::class_<GenericDetector, std::shared_ptr<GenericDetector>>(
                  mex, "GenericDetector")
                  .def(py::init<>())
                  .def("finalize",
                       py::overload_cast<
                           const Config&,
                           std::shared_ptr<const Acts::IMaterialDecorator>>(
                           &GenericDetector::finalize));

    py::class_<Config>(gd, "Config")
        .def(py::init<>())
        .def_readwrite("buildLevel", &Config::buildLevel)
        .def_readwrite("surfaceLogLevel", &Config::surfaceLogLevel)
        .def_readwrite("layerLogLevel", &Config::layerLogLevel)
        .def_readwrite("volumeLogLevel", &Config::volumeLogLevel)
        .def_readwrite("buildProto", &Config::buildProto);
  }

  {
    using Config = AlignedDetector::Config;

    auto d = py::class_<AlignedDetector, std::shared_ptr<AlignedDetector>>(
                 mex, "AlignedDetector")
                 .def(py::init<>())
                 .def("finalize",
                      py::overload_cast<
                          const Config&,
                          std::shared_ptr<const Acts::IMaterialDecorator>>(
                          &AlignedDetector::finalize));

    auto c = py::class_<Config, GenericDetector::Config>(d, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(seed);
    ACTS_PYTHON_MEMBER(iovSize);
    ACTS_PYTHON_MEMBER(flushSize);
    ACTS_PYTHON_MEMBER(sigmaInPlane);
    ACTS_PYTHON_MEMBER(sigmaOutPlane);
    ACTS_PYTHON_MEMBER(sigmaInRot);
    ACTS_PYTHON_MEMBER(sigmaOutRot);
    ACTS_PYTHON_MEMBER(firstIovNominal);
    ACTS_PYTHON_STRUCT_END();
  }
}
}  // namespace Acts::Python