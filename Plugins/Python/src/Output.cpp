// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementWriter.hpp"
#include "ActsExamples/Io/Csv/CsvMultiTrajectoryWriter.hpp"
#include "ActsExamples/Io/Csv/CsvParticleWriter.hpp"
#include "ActsExamples/Io/Csv/CsvPlanarClusterWriter.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitWriter.hpp"
#include "ActsExamples/Io/Csv/CsvTrackingGeometryWriter.hpp"
#include "ActsExamples/Io/NuclearInteractions/RootNuclearInteractionParametersWriter.hpp"
#include "ActsExamples/Io/Performance/CKFPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/SeedingPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFitterPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootBFieldWriter.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackWriter.hpp"
#include "ActsExamples/Io/Root/RootMaterialWriter.hpp"
#include "ActsExamples/Io/Root/RootMeasurementWriter.hpp"
#include "ActsExamples/Io/Root/RootParticleWriter.hpp"
#include "ActsExamples/Io/Root/RootPlanarClusterWriter.hpp"
#include "ActsExamples/Io/Root/RootPropagationStepsWriter.hpp"
#include "ActsExamples/Io/Root/RootSimHitWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackParameterWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryStatesWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectorySummaryWriter.hpp"
#include "ActsExamples/Io/Root/RootVertexPerformanceWriter.hpp"
#include "ActsExamples/Plugins/Obj/ObjPropagationStepsWriter.hpp"
#include "ActsExamples/Plugins/Obj/ObjTrackingGeometryWriter.hpp"
#include "ActsExamples/Validation/ResPlotTool.hpp"

#include <memory>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;

namespace Acts::Python {
void addOutput(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");
  {
    using Writer = ActsExamples::ObjPropagationStepsWriter<Acts::detail::Step>;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 mex, "ObjPropagationStepsWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("cfg"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(collection);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_MEMBER(outputScalor);
    ACTS_PYTHON_MEMBER(outputPrecision);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    auto c = py::class_<ViewConfig>(m, "ViewConfig").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, ViewConfig);
    ACTS_PYTHON_MEMBER(visible);
    ACTS_PYTHON_MEMBER(color);
    ACTS_PYTHON_MEMBER(offset);
    ACTS_PYTHON_MEMBER(lineThickness);
    ACTS_PYTHON_MEMBER(surfaceThickness);
    ACTS_PYTHON_MEMBER(nSegments);
    ACTS_PYTHON_MEMBER(triangulate);
    ACTS_PYTHON_MEMBER(outputName);
    ACTS_PYTHON_STRUCT_END();

    patchKwargsConstructor(c);
  }

  {
    using Writer = ActsExamples::ObjTrackingGeometryWriter;
    auto w = py::class_<Writer, std::shared_ptr<Writer>>(
                 mex, "ObjTrackingGeometryWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def("write", py::overload_cast<const AlgorithmContext&,
                                                 const Acts::TrackingGeometry&>(
                                   &Writer::write));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(outputScalor);
    ACTS_PYTHON_MEMBER(outputPrecision);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_MEMBER(containerView);
    ACTS_PYTHON_MEMBER(volumeView);
    ACTS_PYTHON_MEMBER(sensitiveView);
    ACTS_PYTHON_MEMBER(passiveView);
    ACTS_PYTHON_MEMBER(gridView);
    ACTS_PYTHON_STRUCT_END();
  }

  // ROOT WRITERS

  {
    using Writer = ActsExamples::RootPropagationStepsWriter;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 mex, "RootPropagationStepsWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("cfg"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(collection);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_MEMBER(treeName);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::RootParticleWriter;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 mex, "RootParticleWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("cfg"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_MEMBER(treeName);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::TrackFinderPerformanceWriter;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 mex, "TrackFinderPerformanceWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputProtoTracks);
    ACTS_PYTHON_MEMBER(inputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_MEMBER(treeNameTracks);
    ACTS_PYTHON_MEMBER(treeNameParticles);
    ACTS_PYTHON_STRUCT_END();
  }

  {

  }

  {
    using Writer = ActsExamples::TrackFitterPerformanceWriter;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 mex, "TrackFitterPerformanceWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputTrajectories);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(inputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(resPlotToolConfig);
    ACTS_PYTHON_MEMBER(effPlotToolConfig);
    ACTS_PYTHON_MEMBER(trackSummaryPlotToolConfig);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::SeedingPerformanceWriter;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 mex, "SeedingPerformanceWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("cfg"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputProtoTracks);
    ACTS_PYTHON_MEMBER(inputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_MEMBER(effPlotToolConfig);
    ACTS_PYTHON_MEMBER(duplicationPlotToolConfig);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::RootTrackParameterWriter;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 mex, "RootTrackParameterWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("cfg"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputTrackParameters);
    ACTS_PYTHON_MEMBER(inputProtoTracks);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(inputSimHits);
    ACTS_PYTHON_MEMBER(inputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(inputMeasurementSimHitsMap);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(treeName);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::RootMaterialTrackWriter;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 mex, "RootMaterialTrackWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(collection);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_MEMBER(treeName);
    ACTS_PYTHON_MEMBER(recalculateTotals);
    ACTS_PYTHON_MEMBER(prePostStep);
    ACTS_PYTHON_MEMBER(storeSurface);
    ACTS_PYTHON_MEMBER(storeVolume);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::RootBFieldWriter;
    auto w =
        py::class_<Writer>(mex, "RootBFieldWriter")
            .def_static(
                "run",
                [](const Writer::Config& config, Acts::Logging::Level level) {
                  Writer::run(config, Acts::getDefaultLogger("RootBFieldWriter",
                                                             level));
                },
                py::arg("config"), py::arg("level"));

    py::enum_<Writer::GridType>(w, "GridType")
        .value("rz", Writer::GridType::rz)
        .value("xyz", Writer::GridType::xyz);

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(treeName);
    ACTS_PYTHON_MEMBER(fileName);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_MEMBER(bField);
    ACTS_PYTHON_MEMBER(gridType);
    ACTS_PYTHON_MEMBER(rBounds);
    ACTS_PYTHON_MEMBER(zBounds);
    ACTS_PYTHON_MEMBER(rBins);
    ACTS_PYTHON_MEMBER(zBins);
    ACTS_PYTHON_MEMBER(phiBins);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::RootMeasurementWriter;
    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "RootMeasurementWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    c.def("addBoundIndicesFromDigiConfig",
          [](Writer::Config& self, const DigitizationConfig& digiCfg) {
            self.boundIndices =
                Acts::GeometryHierarchyMap<std::vector<Acts::BoundIndices>>(
                    digiCfg.getBoundIndices());
          });

    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputMeasurements);
    ACTS_PYTHON_MEMBER(inputClusters);
    ACTS_PYTHON_MEMBER(inputSimHits);
    ACTS_PYTHON_MEMBER(inputMeasurementSimHitsMap);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_MEMBER(boundIndices);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_STRUCT_END();
  }

  py::class_<IMaterialWriter, std::shared_ptr<IMaterialWriter>>(
      mex, "IMaterialWriter");

  {
    using Writer = ActsExamples::RootMaterialWriter;
    auto w = py::class_<Writer, IMaterialWriter, std::shared_ptr<Writer>>(
                 mex, "RootMaterialWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def("write", py::overload_cast<const Acts::TrackingGeometry&>(
                                   &Writer::write));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(processSensitives);
    ACTS_PYTHON_MEMBER(processApproaches);
    ACTS_PYTHON_MEMBER(processRepresenting);
    ACTS_PYTHON_MEMBER(processBoundaries);
    ACTS_PYTHON_MEMBER(processVolumes);
    ACTS_PYTHON_MEMBER(folderSurfaceNameBase);
    ACTS_PYTHON_MEMBER(folderVolumeNameBase);
    ACTS_PYTHON_MEMBER(voltag);
    ACTS_PYTHON_MEMBER(boutag);
    ACTS_PYTHON_MEMBER(laytag);
    ACTS_PYTHON_MEMBER(apptag);
    ACTS_PYTHON_MEMBER(sentag);
    ACTS_PYTHON_MEMBER(ntag);
    ACTS_PYTHON_MEMBER(vtag);
    ACTS_PYTHON_MEMBER(otag);
    ACTS_PYTHON_MEMBER(mintag);
    ACTS_PYTHON_MEMBER(maxtag);
    ACTS_PYTHON_MEMBER(ttag);
    ACTS_PYTHON_MEMBER(x0tag);
    ACTS_PYTHON_MEMBER(l0tag);
    ACTS_PYTHON_MEMBER(atag);
    ACTS_PYTHON_MEMBER(ztag);
    ACTS_PYTHON_MEMBER(rhotag);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::RootPlanarClusterWriter;
    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "RootPlanarClusterWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputClusters);
    ACTS_PYTHON_MEMBER(inputSimHits);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_MEMBER(treeName);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::RootSimHitWriter;
    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "RootSimHitWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputSimHits);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_MEMBER(treeName);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::RootTrajectoryStatesWriter;
    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "RootTrajectoryStatesWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputTrajectories);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(inputSimHits);
    ACTS_PYTHON_MEMBER(inputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(inputMeasurementSimHitsMap);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(treeName);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::RootTrajectorySummaryWriter;
    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "RootTrajectorySummaryWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputTrajectories);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(inputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(treeName);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::RootVertexPerformanceWriter;
    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "RootVertexPerformanceWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputAllTruthParticles);
    ACTS_PYTHON_MEMBER(inputSelectedTruthParticles);
    ACTS_PYTHON_MEMBER(inputAssociatedTruthParticles);
    ACTS_PYTHON_MEMBER(inputFittedTracks);
    ACTS_PYTHON_MEMBER(inputVertices);
    ACTS_PYTHON_MEMBER(inputTime);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(treeName);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_MEMBER(minTrackVtxMatchFraction);
    ACTS_PYTHON_STRUCT_END();
  }

  // CSV WRITERS

  {
    using Writer = ActsExamples::CsvParticleWriter;
    auto w = py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
                 mex, "CsvParticleWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_MEMBER(outputStem);
    ACTS_PYTHON_MEMBER(outputPrecision);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::CsvMeasurementWriter;
    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "CsvMeasurementWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputMeasurements);
    ACTS_PYTHON_MEMBER(inputClusters);
    ACTS_PYTHON_MEMBER(inputSimHits);
    ACTS_PYTHON_MEMBER(inputMeasurementSimHitsMap);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_MEMBER(outputPrecision);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::CsvPlanarClusterWriter;
    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "CsvPlanarClusterWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputClusters);
    ACTS_PYTHON_MEMBER(inputSimHits);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_MEMBER(outputPrecision);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::CsvSimHitWriter;
    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "CsvSimHitWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputSimHits);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_MEMBER(outputStem);
    ACTS_PYTHON_MEMBER(outputPrecision);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::CsvMultiTrajectoryWriter;
    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "CsvMultiTrajectoryWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputTrajectories);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_MEMBER(inputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(outputPrecision);
    ACTS_PYTHON_MEMBER(nMeasurementsMin);
    ACTS_PYTHON_MEMBER(truthMatchProbMin);
    ACTS_PYTHON_MEMBER(ptMin);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::CsvTrackingGeometryWriter;
    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "CsvTrackingGeometryWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def("write", &Writer::write);

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_MEMBER(outputPrecision);
    ACTS_PYTHON_MEMBER(writeSensitive);
    ACTS_PYTHON_MEMBER(writeBoundary);
    ACTS_PYTHON_MEMBER(writeSurfaceGrid);
    ACTS_PYTHON_MEMBER(writeLayerVolume);
    ACTS_PYTHON_MEMBER(writePerEvent);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::CKFPerformanceWriter;
    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "CKFPerformanceWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def("write", &Writer::write);

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputTrajectories);
    ACTS_PYTHON_MEMBER(inputParticles);
    ACTS_PYTHON_MEMBER(inputMeasurementParticlesMap);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_MEMBER(effPlotToolConfig);
    ACTS_PYTHON_MEMBER(fakeRatePlotToolConfig);
    ACTS_PYTHON_MEMBER(duplicationPlotToolConfig);
    ACTS_PYTHON_MEMBER(trackSummaryPlotToolConfig);
    ACTS_PYTHON_MEMBER(truthMatchProbMin);
    ACTS_PYTHON_MEMBER(nMeasurementsMin);
    ACTS_PYTHON_MEMBER(ptMin);
    ACTS_PYTHON_MEMBER(duplicatedPredictor);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Writer = ActsExamples::RootNuclearInteractionParametersWriter;
    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "RootNuclearInteractionParametersWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def("write", &Writer::write);

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(inputSimulationProcesses);
    ACTS_PYTHON_MEMBER(filePath);
    ACTS_PYTHON_MEMBER(fileMode);
    ACTS_PYTHON_MEMBER(interactionProbabilityBins);
    ACTS_PYTHON_MEMBER(momentumBins);
    ACTS_PYTHON_MEMBER(invariantMassBins);
    ACTS_PYTHON_MEMBER(multiplicityMax);
    ACTS_PYTHON_MEMBER(writeOptionalHistograms);
    ACTS_PYTHON_MEMBER(nSimulatedEvents);
    ACTS_PYTHON_STRUCT_END();
  }
}
}  // namespace Acts::Python