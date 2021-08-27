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
#include "ActsExamples/TGeoDetector/TGeoDetector.hpp"

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

  {
    using Config = TGeoDetector::Config;

    auto d = py::class_<TGeoDetector, std::shared_ptr<TGeoDetector>>(
                 mex, "TGeoDetector")
                 .def(py::init<>())
                 .def("finalize",
                      py::overload_cast<
                          const Config&,
                          std::shared_ptr<const Acts::IMaterialDecorator>>(
                          &TGeoDetector::finalize));

    py::class_<Options::Interval>(mex, "Interval")
        .def(py::init<>())
        .def(py::init<std::optional<double>, std::optional<double>>())
        .def_readwrite("lower", &Options::Interval::lower)
        .def_readwrite("upper", &Options::Interval::upper);

    auto c = py::class_<Config>(d, "Config").def(py::init<>());

    py::enum_<Config::SubVolume>(c, "SubVolume")
        .value("Negative", Config::SubVolume::Negative)
        .value("Central", Config::SubVolume::Central)
        .value("Positive", Config::SubVolume::Positive);

    auto volume = py::class_<Config::Volume>(c, "Volume").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(volume, Config::Volume);
    ACTS_PYTHON_MEMBER(name);
    ACTS_PYTHON_MEMBER(binToleranceR);
    ACTS_PYTHON_MEMBER(binTolerancePhi);
    ACTS_PYTHON_MEMBER(binToleranceZ);
    ACTS_PYTHON_MEMBER(cylinderDiscSplit);
    ACTS_PYTHON_MEMBER(cylinderNZSegments);
    ACTS_PYTHON_MEMBER(cylinderNPhiSegments);
    ACTS_PYTHON_MEMBER(discNRSegments);
    ACTS_PYTHON_MEMBER(discNPhiSegments);

    ACTS_PYTHON_MEMBER(layers);
    ACTS_PYTHON_MEMBER(subVolumeName);
    ACTS_PYTHON_MEMBER(sensitiveNames);
    ACTS_PYTHON_MEMBER(sensitiveAxes);
    ACTS_PYTHON_MEMBER(rRange);
    ACTS_PYTHON_MEMBER(zRange);
    ACTS_PYTHON_MEMBER(splitTolR);
    ACTS_PYTHON_MEMBER(splitTolZ);
    ACTS_PYTHON_STRUCT_END();

    auto regTriplet = [&c](const std::string& name, auto v) {
      using type = decltype(v);
      py::class_<Config::LayerTriplet<type>>(c, name.c_str())
          .def(py::init<>())
          .def(py::init<type>())
          .def(py::init<type, type, type>())
          .def_readwrite("negative", &Config::LayerTriplet<type>::negative)
          .def_readwrite("central", &Config::LayerTriplet<type>::central)
          .def_readwrite("positive", &Config::LayerTriplet<type>::positive)
          .def("at", py::overload_cast<Config::SubVolume>(
                         &Config::LayerTriplet<type>::at));
    };

    regTriplet("LayerTripletBool", true);
    regTriplet("LayerTripletString", std::string{""});
    regTriplet("LayerTripletVectorString", std::vector<std::string>{});
    regTriplet("LayerTripletInterval", Options::Interval{});
    regTriplet("LayerTripletDouble", double{5.5});

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(surfaceLogLevel);
    ACTS_PYTHON_MEMBER(layerLogLevel);
    ACTS_PYTHON_MEMBER(volumeLogLevel);
    ACTS_PYTHON_MEMBER(fileName);
    ACTS_PYTHON_MEMBER(buildBeamPipe);
    ACTS_PYTHON_MEMBER(beamPipeRadius);
    ACTS_PYTHON_MEMBER(beamPipeHalflengthZ);
    ACTS_PYTHON_MEMBER(beamPipeLayerThickness);
    ACTS_PYTHON_MEMBER(unitScalor);
    ACTS_PYTHON_MEMBER(volumes);
    ACTS_PYTHON_STRUCT_END();

    patchKwargsConstructor(c);
  }
}
}  // namespace Acts::Python