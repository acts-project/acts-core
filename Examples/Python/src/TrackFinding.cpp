// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"
#include "Acts/Seeding/SeedFilterConfig.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SeedFinderGbtsConfig.hpp"
#include "Acts/Seeding/SeedFinderOrthogonalConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SpacePointContainer.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/TrackFinding/GbtsSeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/HoughTransformSeeder.hpp"
#include "ActsExamples/TrackFinding/MuonHoughSeeder.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SeedingOrthogonalAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackParamsLookupEstimation.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace Acts {
class MagneticFieldProvider;
class TrackingGeometry;
}  // namespace Acts
namespace ActsExamples {
class IAlgorithm;
class SimSpacePoint;
}  // namespace ActsExamples

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addTrackFinding(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::SpacePointMaker, mex, "SpacePointMaker", inputMeasurements,
      outputSpacePoints, trackingGeometry, geometrySelection);

  {
    using Config = Acts::SeedFilterConfig;
    auto c = py::class_<Config>(m, "SeedFilterConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(deltaInvHelixDiameter);
    ACTS_PYTHON_MEMBER(impactWeightFactor);
    ACTS_PYTHON_MEMBER(zOriginWeightFactor);
    ACTS_PYTHON_MEMBER(compatSeedWeight);
    ACTS_PYTHON_MEMBER(deltaRMin);
    ACTS_PYTHON_MEMBER(maxSeedsPerSpM);
    ACTS_PYTHON_MEMBER(compatSeedLimit);
    ACTS_PYTHON_MEMBER(seedConfirmation);
    ACTS_PYTHON_MEMBER(centralSeedConfirmationRange);
    ACTS_PYTHON_MEMBER(forwardSeedConfirmationRange);
    ACTS_PYTHON_MEMBER(useDeltaRorTopRadius);
    ACTS_PYTHON_MEMBER(seedWeightIncrement);
    ACTS_PYTHON_MEMBER(numSeedIncrement);
    ACTS_PYTHON_MEMBER(maxSeedsPerSpMConf);
    ACTS_PYTHON_MEMBER(maxQualitySeedsPerSpMConf);
    ACTS_PYTHON_STRUCT_END();
    patchKwargsConstructor(c);
  }

  {
    using Config = Acts::SeedFinderConfig<typename Acts::SpacePointContainer<
        ActsExamples::SpacePointContainer<std::vector<const SimSpacePoint*>>,
        Acts::detail::RefHolder>::SpacePointProxyType>;
    auto c = py::class_<Config>(m, "SeedFinderConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(minPt);
    ACTS_PYTHON_MEMBER(cotThetaMax);
    ACTS_PYTHON_MEMBER(deltaRMin);
    ACTS_PYTHON_MEMBER(deltaRMax);
    ACTS_PYTHON_MEMBER(deltaRMinBottomSP);
    ACTS_PYTHON_MEMBER(deltaRMaxBottomSP);
    ACTS_PYTHON_MEMBER(deltaRMinTopSP);
    ACTS_PYTHON_MEMBER(deltaRMaxTopSP);
    ACTS_PYTHON_MEMBER(impactMax);
    ACTS_PYTHON_MEMBER(sigmaScattering);
    ACTS_PYTHON_MEMBER(maxPtScattering);
    ACTS_PYTHON_MEMBER(maxSeedsPerSpM);
    ACTS_PYTHON_MEMBER(collisionRegionMin);
    ACTS_PYTHON_MEMBER(collisionRegionMax);
    ACTS_PYTHON_MEMBER(phiMin);
    ACTS_PYTHON_MEMBER(phiMax);
    ACTS_PYTHON_MEMBER(zMin);
    ACTS_PYTHON_MEMBER(zMax);
    ACTS_PYTHON_MEMBER(rMax);
    ACTS_PYTHON_MEMBER(rMin);
    ACTS_PYTHON_MEMBER(radLengthPerSeed);
    ACTS_PYTHON_MEMBER(zAlign);
    ACTS_PYTHON_MEMBER(rAlign);
    ACTS_PYTHON_MEMBER(sigmaError);
    ACTS_PYTHON_MEMBER(maxBlockSize);
    ACTS_PYTHON_MEMBER(nTrplPerSpBLimit);
    ACTS_PYTHON_MEMBER(nAvgTrplPerSpBLimit);
    ACTS_PYTHON_MEMBER(impactMax);
    ACTS_PYTHON_MEMBER(deltaZMax);
    ACTS_PYTHON_MEMBER(zBinEdges);
    ACTS_PYTHON_MEMBER(interactionPointCut);
    ACTS_PYTHON_MEMBER(zBinsCustomLooping);
    ACTS_PYTHON_MEMBER(useVariableMiddleSPRange);
    ACTS_PYTHON_MEMBER(deltaRMiddleMinSPRange);
    ACTS_PYTHON_MEMBER(deltaRMiddleMaxSPRange);
    ACTS_PYTHON_MEMBER(rRangeMiddleSP);
    ACTS_PYTHON_MEMBER(rMinMiddle);
    ACTS_PYTHON_MEMBER(rMaxMiddle);
    ACTS_PYTHON_MEMBER(binSizeR);
    ACTS_PYTHON_MEMBER(seedConfirmation);
    ACTS_PYTHON_MEMBER(centralSeedConfirmationRange);
    ACTS_PYTHON_MEMBER(forwardSeedConfirmationRange);
    ACTS_PYTHON_MEMBER(useDetailedDoubleMeasurementInfo);
    ACTS_PYTHON_STRUCT_END();
    patchKwargsConstructor(c);
  }
  {
    using seedOptions = Acts::SeedFinderOptions;
    auto c = py::class_<seedOptions>(m, "SeedFinderOptions").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, seedOptions);
    ACTS_PYTHON_MEMBER(beamPos);
    ACTS_PYTHON_MEMBER(bFieldInZ);
    ACTS_PYTHON_STRUCT_END();
    patchKwargsConstructor(c);
  }
  {
    using Config =
        Acts::SeedFinderOrthogonalConfig<typename Acts::SpacePointContainer<
            ActsExamples::SpacePointContainer<
                std::vector<const SimSpacePoint*>>,
            Acts::detail::RefHolder>::SpacePointProxyType>;
    auto c =
        py::class_<Config>(m, "SeedFinderOrthogonalConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(minPt);
    ACTS_PYTHON_MEMBER(cotThetaMax);
    ACTS_PYTHON_MEMBER(deltaRMinBottomSP);
    ACTS_PYTHON_MEMBER(deltaRMaxBottomSP);
    ACTS_PYTHON_MEMBER(deltaRMinTopSP);
    ACTS_PYTHON_MEMBER(deltaRMaxTopSP);
    ACTS_PYTHON_MEMBER(impactMax);
    ACTS_PYTHON_MEMBER(deltaPhiMax);
    ACTS_PYTHON_MEMBER(deltaZMax);
    ACTS_PYTHON_MEMBER(sigmaScattering);
    ACTS_PYTHON_MEMBER(maxPtScattering);
    ACTS_PYTHON_MEMBER(maxSeedsPerSpM);
    ACTS_PYTHON_MEMBER(collisionRegionMin);
    ACTS_PYTHON_MEMBER(collisionRegionMax);
    ACTS_PYTHON_MEMBER(phiMin);
    ACTS_PYTHON_MEMBER(phiMax);
    ACTS_PYTHON_MEMBER(zMin);
    ACTS_PYTHON_MEMBER(zMax);
    ACTS_PYTHON_MEMBER(rMax);
    ACTS_PYTHON_MEMBER(rMin);
    ACTS_PYTHON_MEMBER(radLengthPerSeed);
    ACTS_PYTHON_MEMBER(deltaZMax);
    ACTS_PYTHON_MEMBER(interactionPointCut);
    ACTS_PYTHON_MEMBER(deltaPhiMax);
    ACTS_PYTHON_MEMBER(highland);
    ACTS_PYTHON_MEMBER(maxScatteringAngle2);
    ACTS_PYTHON_MEMBER(useVariableMiddleSPRange);
    ACTS_PYTHON_MEMBER(deltaRMiddleMinSPRange);
    ACTS_PYTHON_MEMBER(deltaRMiddleMaxSPRange);
    ACTS_PYTHON_MEMBER(rRangeMiddleSP);
    ACTS_PYTHON_MEMBER(rMinMiddle);
    ACTS_PYTHON_MEMBER(rMaxMiddle);
    ACTS_PYTHON_MEMBER(seedConfirmation);
    ACTS_PYTHON_MEMBER(centralSeedConfirmationRange);
    ACTS_PYTHON_MEMBER(forwardSeedConfirmationRange);
    ACTS_PYTHON_STRUCT_END();
    patchKwargsConstructor(c);
  }

  {
    using Config = Acts::SeedFinderGbtsConfig<SimSpacePoint>;
    auto c = py::class_<Config>(m, "SeedFinderGbtsConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(minPt);
    ACTS_PYTHON_MEMBER(sigmaScattering);
    ACTS_PYTHON_MEMBER(highland);
    ACTS_PYTHON_MEMBER(maxScatteringAngle2);
    ACTS_PYTHON_MEMBER(ConnectorInputFile);
    ACTS_PYTHON_MEMBER(m_phiSliceWidth);
    ACTS_PYTHON_MEMBER(m_nMaxPhiSlice);
    ACTS_PYTHON_MEMBER(m_useClusterWidth);
    ACTS_PYTHON_MEMBER(m_layerGeometry);
    ACTS_PYTHON_STRUCT_END();
    patchKwargsConstructor(c);
  }

  {
    using seedConf = Acts::SeedConfirmationRangeConfig;
    auto c = py::class_<seedConf>(m, "SeedConfirmationRangeConfig")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, seedConf);
    ACTS_PYTHON_MEMBER(zMinSeedConf);
    ACTS_PYTHON_MEMBER(zMaxSeedConf);
    ACTS_PYTHON_MEMBER(rMaxSeedConf);
    ACTS_PYTHON_MEMBER(nTopForLargeR);
    ACTS_PYTHON_MEMBER(nTopForSmallR);
    ACTS_PYTHON_MEMBER(seedConfMinBottomRadius);
    ACTS_PYTHON_MEMBER(seedConfMaxZOrigin);
    ACTS_PYTHON_MEMBER(minImpactSeedConf);
    ACTS_PYTHON_STRUCT_END();
    patchKwargsConstructor(c);
  }

  {
    using Config = Acts::CylindricalSpacePointGridConfig;
    auto c = py::class_<Config>(m, "SpacePointGridConfig").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(minPt);
    ACTS_PYTHON_MEMBER(rMax);
    ACTS_PYTHON_MEMBER(zMax);
    ACTS_PYTHON_MEMBER(zMin);
    ACTS_PYTHON_MEMBER(phiMin);
    ACTS_PYTHON_MEMBER(phiMax);
    ACTS_PYTHON_MEMBER(deltaRMax);
    ACTS_PYTHON_MEMBER(cotThetaMax);
    ACTS_PYTHON_MEMBER(phiBinDeflectionCoverage);
    ACTS_PYTHON_MEMBER(maxPhiBins);
    ACTS_PYTHON_MEMBER(impactMax);
    ACTS_PYTHON_MEMBER(zBinEdges);
    ACTS_PYTHON_STRUCT_END();
    patchKwargsConstructor(c);
  }
  {
    using Options = Acts::CylindricalSpacePointGridOptions;
    auto c = py::class_<Options>(m, "SpacePointGridOptions").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Options);
    ACTS_PYTHON_MEMBER(bFieldInZ);
    ACTS_PYTHON_STRUCT_END();
    patchKwargsConstructor(c);
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::SeedingAlgorithm, mex, "SeedingAlgorithm", inputSpacePoints,
      outputSeeds, seedFilterConfig, seedFinderConfig, seedFinderOptions,
      gridConfig, gridOptions, allowSeparateRMax, zBinNeighborsTop,
      zBinNeighborsBottom, numPhiNeighbors, useExtraCuts);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::SeedingOrthogonalAlgorithm, mex,
                                "SeedingOrthogonalAlgorithm", inputSpacePoints,
                                outputSeeds, seedFilterConfig, seedFinderConfig,
                                seedFinderOptions);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::GbtsSeedingAlgorithm, mex, "GbtsSeedingAlgorithm",
      inputSpacePoints, outputSeeds, seedFinderConfig, seedFinderOptions,
      layerMappingFile, geometrySelection, trackingGeometry, ActsGbtsMap,
      fill_module_csv, inputClusters);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::HoughTransformSeeder, mex, "HoughTransformSeeder",
      inputSpacePoints, outputProtoTracks, trackingGeometry, geometrySelection,
      inputMeasurements, subRegions, nLayers, xMin, xMax, yMin, yMax,
      houghHistSize_x, houghHistSize_y, hitExtend_x, threshold,
      localMaxWindowSize, kA);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::MuonHoughSeeder, mex,
                                "MuonHoughSeeder", inSimHits, inDriftCircles);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::TrackParamsEstimationAlgorithm, mex,
      "TrackParamsEstimationAlgorithm", inputSeeds, inputProtoTracks,
      outputTrackParameters, outputSeeds, outputProtoTracks, trackingGeometry,
      magneticField, bFieldMin, initialSigmas, initialSigmaPtRel,
      initialVarInflation, noTimeVarInflation, particleHypothesis);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::TrackParamsLookupEstimation, mex,
                                "TrackParamsLookupEstimation", refLayers, bins,
                                inputHits, inputParticles,
                                trackLookupGridWriters);

  {
    using Alg = ActsExamples::TrackFindingAlgorithm;
    using Config = Alg::Config;

    auto alg =
        py::class_<Alg, ActsExamples::IAlgorithm, std::shared_ptr<Alg>>(
            mex, "TrackFindingAlgorithm")
            .def(py::init<const Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Alg::config)
            .def_static("makeTrackFinderFunction",
                        [](std::shared_ptr<const Acts::TrackingGeometry>
                               trackingGeometry,
                           std::shared_ptr<const Acts::MagneticFieldProvider>
                               magneticField,
                           Logging::Level level) {
                          return Alg::makeTrackFinderFunction(
                              std::move(trackingGeometry),
                              std::move(magneticField),
                              *Acts::getDefaultLogger("TrackFinding", level));
                        });

    py::class_<Alg::TrackFinderFunction,
               std::shared_ptr<Alg::TrackFinderFunction>>(
        alg, "TrackFinderFunction");

    auto c = py::class_<Config>(alg, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(inputMeasurements);
    ACTS_PYTHON_MEMBER(inputInitialTrackParameters);
    ACTS_PYTHON_MEMBER(inputSeeds);
    ACTS_PYTHON_MEMBER(outputTracks);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(magneticField);
    ACTS_PYTHON_MEMBER(findTracks);
    ACTS_PYTHON_MEMBER(measurementSelectorCfg);
    ACTS_PYTHON_MEMBER(trackSelectorCfg);
    ACTS_PYTHON_MEMBER(maxSteps);
    ACTS_PYTHON_MEMBER(twoWay);
    ACTS_PYTHON_MEMBER(reverseSearch);
    ACTS_PYTHON_MEMBER(seedDeduplication);
    ACTS_PYTHON_MEMBER(stayOnSeed);
    ACTS_PYTHON_MEMBER(pixelVolumeIds);
    ACTS_PYTHON_MEMBER(stripVolumeIds);
    ACTS_PYTHON_MEMBER(maxPixelHoles);
    ACTS_PYTHON_MEMBER(maxStripHoles);
    ACTS_PYTHON_MEMBER(trimTracks);
    ACTS_PYTHON_MEMBER(constrainToVolumeIds);
    ACTS_PYTHON_MEMBER(endOfWorldVolumeIds);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    auto constructor =
        [](const std::vector<std::pair<
               GeometryIdentifier,
               std::tuple<std::vector<double>, std::vector<double>,
                          std::vector<double>, std::vector<std::size_t>>>>&
               input) {
          std::vector<std::pair<GeometryIdentifier, MeasurementSelectorCuts>>
              converted;
          converted.reserve(input.size());
          for (const auto& [id, cuts] : input) {
            const auto& [bins, chi2Measurement, chi2Outlier, num] = cuts;
            converted.emplace_back(
                id, MeasurementSelectorCuts{bins, chi2Measurement, num,
                                            chi2Outlier});
          }
          return std::make_unique<MeasurementSelector::Config>(converted);
        };

    py::class_<MeasurementSelectorCuts>(m, "MeasurementSelectorCuts")
        .def(py::init<>())
        .def(py::init<std::vector<double>, std::vector<double>,
                      std::vector<std::size_t>, std::vector<double>>())
        .def_readwrite("etaBins", &MeasurementSelectorCuts::etaBins)
        .def_readwrite("chi2CutOffMeasurement",
                       &MeasurementSelectorCuts::chi2CutOff)
        .def_readwrite("chi2CutOffOutlier",
                       &MeasurementSelectorCuts::chi2CutOffOutlier)
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
