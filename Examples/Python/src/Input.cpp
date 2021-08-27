// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvPlanarClusterReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackReader.hpp"
#include "ActsExamples/Io/Root/RootParticleReader.hpp"
#include "ActsExamples/Io/Root/RootTrajectorySummaryReader.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;

namespace Acts::Python {
void addInput(Context& ctx) {
  auto mex = ctx.get("examples");

  // ROOT READERS

  {
    using Reader = ActsExamples::RootParticleReader;
    using Config = Reader::Config;
    auto reader =
        py::class_<Reader, ActsExamples::IReader, std::shared_ptr<Reader>>(
            mex, "RootParticleReader")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Reader::config);

    auto c = py::class_<Config>(reader, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(particleCollection);
    ACTS_PYTHON_MEMBER(vertexPrimaryCollection);
    ACTS_PYTHON_MEMBER(vertexSecondaryCollection);
    ACTS_PYTHON_MEMBER(treeName);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(orderedEvents);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Reader = ActsExamples::RootMaterialTrackReader;
    using Config = Reader::Config;
    auto reader =
        py::class_<Reader, ActsExamples::IReader, std::shared_ptr<Reader>>(
            mex, "RootMaterialTrackReader")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Reader::config);

    auto c = py::class_<Config>(reader, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(collection);
    ACTS_PYTHON_MEMBER(treeName);
    ACTS_PYTHON_MEMBER(fileList);
    ACTS_PYTHON_MEMBER(orderedEvents);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Reader = ActsExamples::RootTrajectorySummaryReader;
    using Config = Reader::Config;
    auto reader =
        py::class_<Reader, ActsExamples::IReader, std::shared_ptr<Reader>>(
            mex, "RootTrajectorySummaryReader")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Reader::config);

    auto c = py::class_<Config>(reader, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(outputTracks);
    ACTS_PYTHON_MEMBER(outputParticles);
    ACTS_PYTHON_MEMBER(treeName);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(orderedEvents);
    ACTS_PYTHON_STRUCT_END();
  }

  // CSV READERS

  {
    using Reader = ActsExamples::CsvParticleReader;
    using Config = Reader::Config;
    auto reader =
        py::class_<Reader, ActsExamples::IReader, std::shared_ptr<Reader>>(
            mex, "CsvParticleReader")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Reader::config);

    auto c = py::class_<Config>(reader, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputDir);
    ACTS_PYTHON_MEMBER(inputStem);
    ACTS_PYTHON_MEMBER(outputParticles);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Reader = ActsExamples::CsvMeasurementReader;
    using Config = Reader::Config;
    auto reader =
        py::class_<Reader, ActsExamples::IReader, std::shared_ptr<Reader>>(
            mex, "CsvMeasurementReader")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Reader::config);

    auto c = py::class_<Config>(reader, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputDir);
    ACTS_PYTHON_MEMBER(outputMeasurements);
    ACTS_PYTHON_MEMBER(outputMeasurementSimHitsMap);
    ACTS_PYTHON_MEMBER(outputSourceLinks);
    ACTS_PYTHON_MEMBER(outputClusters);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Reader = ActsExamples::CsvPlanarClusterReader;
    using Config = Reader::Config;
    auto reader =
        py::class_<Reader, ActsExamples::IReader, std::shared_ptr<Reader>>(
            mex, "CsvPlanarClusterReader")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Reader::config);

    auto c = py::class_<Config>(reader, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputDir);
    ACTS_PYTHON_MEMBER(outputClusters);
    ACTS_PYTHON_MEMBER(outputHitIds);
    ACTS_PYTHON_MEMBER(outputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(outputSimHits);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Reader = ActsExamples::CsvSimHitReader;
    using Config = Reader::Config;
    auto reader =
        py::class_<Reader, ActsExamples::IReader, std::shared_ptr<Reader>>(
            mex, "CsvSimHitReader")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Reader::config);

    auto c = py::class_<Config>(reader, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputDir);
    ACTS_PYTHON_MEMBER(inputStem);
    ACTS_PYTHON_MEMBER(outputSimHits);
    ACTS_PYTHON_STRUCT_END();
  }
}
}  // namespace Acts::Python