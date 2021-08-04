// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addTrackFinding(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  {
    using Config = ActsExamples::SpacePointMaker::Config;
    auto alg =
        py::class_<ActsExamples::SpacePointMaker, ActsExamples::BareAlgorithm,
                   std::shared_ptr<ActsExamples::SpacePointMaker>>(
            mex, "SpacePointMaker")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config",
                                   &ActsExamples::SpacePointMaker::config);

    auto c = py::class_<ActsExamples::SpacePointMaker::Config>(alg, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputSourceLinks);
    ACTS_PYTHON_MEMBER(inputMeasurements);
    ACTS_PYTHON_MEMBER(outputSpacePoints);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(geometrySelection);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Config = ActsExamples::SeedingAlgorithm::Config;

    auto alg =
        py::class_<ActsExamples::SeedingAlgorithm, ActsExamples::BareAlgorithm,
                   std::shared_ptr<ActsExamples::SeedingAlgorithm>>(
            mex, "SeedingAlgorithm")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config",
                                   &ActsExamples::SeedingAlgorithm::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputSpacePoints);
    ACTS_PYTHON_MEMBER(outputSeeds);
    ACTS_PYTHON_MEMBER(outputProtoTracks);
    ACTS_PYTHON_MEMBER(rMax);
    ACTS_PYTHON_MEMBER(deltaRMin);
    ACTS_PYTHON_MEMBER(deltaRMax);
    ACTS_PYTHON_MEMBER(collisionRegionMin);
    ACTS_PYTHON_MEMBER(collisionRegionMax);
    ACTS_PYTHON_MEMBER(zMin);
    ACTS_PYTHON_MEMBER(zMax);
    ACTS_PYTHON_MEMBER(maxSeedsPerSpM);
    ACTS_PYTHON_MEMBER(cotThetaMax);
    ACTS_PYTHON_MEMBER(sigmaScattering);
    ACTS_PYTHON_MEMBER(radLengthPerSeed);
    ACTS_PYTHON_MEMBER(minPt);
    ACTS_PYTHON_MEMBER(bFieldInZ);
    ACTS_PYTHON_MEMBER(beamPosX);
    ACTS_PYTHON_MEMBER(beamPosY);
    ACTS_PYTHON_MEMBER(impactMax);
    ACTS_PYTHON_STRUCT_END();

    c.def_property(
        "deltaR",
        [](Config& cfg) {
          return std::pair{cfg.deltaRMin, cfg.deltaRMax};
        },
        [](Config& cfg, std::pair<double, double> values) {
          cfg.deltaRMin = values.first;
          cfg.deltaRMax = values.second;
        });
    c.def_property(
        "z",
        [](Config& cfg) {
          return std::pair{cfg.zMin, cfg.zMax};
        },
        [](Config& cfg, std::pair<double, double> values) {
          cfg.zMin = values.first;
          cfg.zMax = values.second;
        });

    c.def_property(
        "collisionRegion",
        [](Config& cfg) {
          return std::pair{cfg.collisionRegionMin, cfg.collisionRegionMax};
        },
        [](Config& cfg, std::pair<double, double> values) {
          cfg.collisionRegionMin = values.first;
          cfg.collisionRegionMax = values.second;
        });
  }

  {
    using Alg = ActsExamples::TrackParamsEstimationAlgorithm;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            mex, "TrackParamsEstimationAlgorithm")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputSeeds);
    ACTS_PYTHON_MEMBER(inputSpacePoints);
    ACTS_PYTHON_MEMBER(inputProtoTracks);
    ACTS_PYTHON_MEMBER(inputSourceLinks);
    ACTS_PYTHON_MEMBER(outputTrackParameters);
    ACTS_PYTHON_MEMBER(outputProtoTracks);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(magneticField);
    ACTS_PYTHON_MEMBER(deltaRMin);
    ACTS_PYTHON_MEMBER(deltaRMax);
    ACTS_PYTHON_MEMBER(bFieldMin);
    ACTS_PYTHON_MEMBER(sigmaLoc0);
    ACTS_PYTHON_MEMBER(sigmaLoc1);
    ACTS_PYTHON_MEMBER(sigmaPhi);
    ACTS_PYTHON_MEMBER(sigmaTheta);
    ACTS_PYTHON_MEMBER(sigmaQOverP);
    ACTS_PYTHON_MEMBER(sigmaT0);
    ACTS_PYTHON_MEMBER(initialVarInflation);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Alg = ActsExamples::TrackFindingAlgorithm;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            mex, "TrackFindingAlgorithm")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Alg::config)
            .def_static("makeTrackFinderFunction",
                        &Alg::makeTrackFinderFunction);

    py::class_<Alg::TrackFinderFunction,
               std::shared_ptr<Alg::TrackFinderFunction>>(
        alg, "TrackFinderFunction");

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputMeasurements);
    ACTS_PYTHON_MEMBER(inputSourceLinks);
    ACTS_PYTHON_MEMBER(inputInitialTrackParameters);
    ACTS_PYTHON_MEMBER(outputTrajectories);
    ACTS_PYTHON_MEMBER(findTracks);
    ACTS_PYTHON_MEMBER(measurementSelectorCfg);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    auto constructor =
        [](std::vector<std::pair<GeometryIdentifier, std::pair<double, size_t>>>
               input) {
          std::vector<std::pair<GeometryIdentifier, MeasurementSelectorCuts>>
              converted;
          converted.reserve(input.size());
          for (const auto& [id, cuts] : input) {
            const auto [chi2, num] = cuts;
            converted.emplace_back(id, MeasurementSelectorCuts{chi2, num});
          }
          return std::make_unique<MeasurementSelector::Config>(converted);
        };

    py::class_<MeasurementSelectorCuts>(m, "MeasurementSelectorCuts")
        .def(py::init<>())
        .def(py::init<double, size_t>())
        .def_readwrite("chi2CutOff", &MeasurementSelectorCuts::chi2CutOff)
        .def_readwrite("numMeasurementsCutOff",
                       &MeasurementSelectorCuts::numMeasurementsCutOff);

    auto ms = py::class_<MeasurementSelector>(m, "MeasurementSelector");
    auto c =
        py::class_<MeasurementSelector::Config>(ms, "Config")
            .def(py::init<std::vector<
                     std::pair<GeometryIdentifier, MeasurementSelectorCuts>>>())
            .def(py::init(constructor));
  }
}

}  // namespace Acts::Python