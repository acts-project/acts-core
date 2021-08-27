// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Plugins/Json/JsonMaterialDecorator.hpp"
#include "Acts/Plugins/Json/MaterialMapJsonConverter.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Io/Json/JsonMaterialWriter.hpp"
#include "ActsExamples/Io/Json/JsonSurfacesWriter.hpp"

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;

namespace Acts::Python {
void addJson(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  {
    py::class_<JsonMaterialDecorator, Acts::IMaterialDecorator,
               std::shared_ptr<JsonMaterialDecorator>>(m,
                                                       "JsonMaterialDecorator")
        .def(py::init<const MaterialMapJsonConverter::Config&,
                      const std::string&, Acts::Logging::Level, bool, bool>(),
             py::arg("rConfig"), py::arg("jFileName"), py::arg("level"),
             py::arg("clearSurfaceMaterial") = true,
             py::arg("clearVolumeMaterial") = true);
  }

  {
    auto cls =
        py::class_<MaterialMapJsonConverter>(m, "MaterialMapJsonConverter")
            .def(py::init<const MaterialMapJsonConverter::Config&,
                          Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    auto c = py::class_<MaterialMapJsonConverter::Config>(cls, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, MaterialMapJsonConverter::Config);
    ACTS_PYTHON_MEMBER(context);
    ACTS_PYTHON_MEMBER(processSensitives);
    ACTS_PYTHON_MEMBER(processApproaches);
    ACTS_PYTHON_MEMBER(processRepresenting);
    ACTS_PYTHON_MEMBER(processBoundaries);
    ACTS_PYTHON_MEMBER(processVolumes);
    ACTS_PYTHON_MEMBER(processDenseVolumes);
    ACTS_PYTHON_MEMBER(processNonMaterial);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    py::enum_<JsonFormat>(mex, "JsonFormat")
        .value("NoOutput", JsonFormat::NoOutput)
        .value("Json", JsonFormat::Json)
        .value("Cbor", JsonFormat::Cbor)
        .value("All", JsonFormat::All);
  }

  {
    auto cls =
        py::class_<JsonMaterialWriter, IMaterialWriter,
                   std::shared_ptr<JsonMaterialWriter>>(mex,
                                                        "JsonMaterialWriter")
            .def(py::init<const JsonMaterialWriter::Config&,
                          Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def("writeMaterial", &JsonMaterialWriter::writeMaterial)
            .def("write", &JsonMaterialWriter::write)
            .def_property_readonly("config", &JsonMaterialWriter::config);

    auto c =
        py::class_<JsonMaterialWriter::Config>(cls, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, JsonMaterialWriter::Config);
    ACTS_PYTHON_MEMBER(converterCfg);
    ACTS_PYTHON_MEMBER(fileName);
    ACTS_PYTHON_MEMBER(writeFormat);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    auto cls =
        py::class_<JsonSurfacesWriter, IWriter,
                   std::shared_ptr<JsonSurfacesWriter>>(mex,
                                                        "JsonSurfacesWriter")
            .def(py::init<const JsonSurfacesWriter::Config&,
                          Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def("write", &JsonSurfacesWriter::write)
            .def_property_readonly("config", &JsonSurfacesWriter::config);

    auto c =
        py::class_<JsonSurfacesWriter::Config>(cls, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, JsonSurfacesWriter::Config);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_MEMBER(outputPrecision);
    ACTS_PYTHON_MEMBER(writeLayer);
    ACTS_PYTHON_MEMBER(writeApproach);
    ACTS_PYTHON_MEMBER(writeSensitive);
    ACTS_PYTHON_MEMBER(writeBoundary);
    ACTS_PYTHON_MEMBER(writePerEvent);
    ACTS_PYTHON_MEMBER(writeOnlyNames);
    ACTS_PYTHON_STRUCT_END();
  }
}
}  // namespace Acts::Python