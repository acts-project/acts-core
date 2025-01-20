// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/Framework/BufferedReader.hpp"
#include "ActsExamples/Io/Csv/CsvDriftCircleReader.hpp"
#include "ActsExamples/Io/Csv/CsvExaTrkXGraphReader.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementReader.hpp"
#include "ActsExamples/Io/Csv/CsvMuonSimHitReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Csv/CsvSpacePointReader.hpp"
#include "ActsExamples/Io/Csv/CsvTrackParameterReader.hpp"
#include "ActsExamples/Io/Root/RootAthenaDumpReader.hpp"
#include "ActsExamples/Io/Root/RootAthenaNTupleReader.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackReader.hpp"
#include "ActsExamples/Io/Root/RootParticleReader.hpp"
#include "ActsExamples/Io/Root/RootSimHitReader.hpp"
#include "ActsExamples/Io/Root/RootTrackSummaryReader.hpp"
#include "ActsExamples/Io/Root/RootVertexReader.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;

namespace Acts::Python {

void addInput(Context& ctx) {
  auto mex = ctx.get("examples");

  // Buffered reader
  ACTS_PYTHON_DECLARE_READER(ActsExamples::BufferedReader, mex,
                             "BufferedReader", upstreamReader, selectionSeed,
                             bufferSize);

  // ROOT READERS
  ACTS_PYTHON_DECLARE_READER(ActsExamples::RootParticleReader, mex,
                             "RootParticleReader", outputParticles, treeName,
                             filePath);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::RootVertexReader, mex,
                             "RootVertexReader", outputVertices, treeName,
                             filePath);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::RootMaterialTrackReader, mex,
                             "RootMaterialTrackReader", outputMaterialTracks,
                             treeName, fileList, readCachedSurfaceInformation);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::RootTrackSummaryReader, mex,
                             "RootTrackSummaryReader", outputTracks,
                             outputParticles, treeName, filePath);

  // CSV READERS
  ACTS_PYTHON_DECLARE_READER(ActsExamples::CsvParticleReader, mex,
                             "CsvParticleReader", inputDir, inputStem,
                             outputParticles);

  ACTS_PYTHON_DECLARE_READER(
      ActsExamples::CsvMeasurementReader, mex, "CsvMeasurementReader", inputDir,
      outputMeasurements, outputMeasurementSimHitsMap, outputClusters,
      outputMeasurementParticlesMap, inputSimHits);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::CsvSimHitReader, mex,
                             "CsvSimHitReader", inputDir, inputStem,
                             outputSimHits);
  ACTS_PYTHON_DECLARE_READER(ActsExamples::CsvMuonSimHitReader, mex,
                             "CsvMuonSimHitReader", inputDir, inputStem,
                             outputSimHits);
  ACTS_PYTHON_DECLARE_READER(ActsExamples::CsvDriftCircleReader, mex,
                             "CsvDriftCircleReader", inputDir, inputStem,
                             outputDriftCircles);

  ACTS_PYTHON_DECLARE_READER(
      ActsExamples::CsvSpacePointReader, mex, "CsvSpacePointReader", inputDir,
      inputStem, inputCollection, outputSpacePoints, extendCollection);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::CsvTrackParameterReader, mex,
                             "CsvTrackParameterReader", inputDir, inputStem,
                             outputTrackParameters, beamspot);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::RootAthenaNTupleReader, mex,
                             "RootAthenaNTupleReader", inputTreeName,
                             inputFilePath, outputTrackParameters,
                             outputTruthVtxParameters, outputRecoVtxParameters,
                             outputBeamspotConstraint);

  ACTS_PYTHON_DECLARE_READER(
      ActsExamples::RootAthenaDumpReader, mex, "RootAthenaDumpReader", treename,
      inputfiles, outputMeasurements, outputPixelSpacePoints,
      outputStripSpacePoints, outputSpacePoints, outputClusters,
      outputMeasurementParticlesMap, outputParticleMeasurementsMap,
      outputParticles, onlyPassedParticles, skipOverlapSPsPhi,
      skipOverlapSPsEta, geometryIdMap, trackingGeometry, absBoundaryTolerance,
      onlySpacepoints, noTruth);

#ifdef WITH_GEOMODEL_PLUGIN
  ACTS_PYTHON_DECLARE_READER(ActsExamples::RootAthenaDumpGeoIdCollector, mex,
                             "RootAthenaDumpGeoIdCollector", treename,
                             inputfile, trackingGeometry, geometryIdMap);
#endif

  ACTS_PYTHON_DECLARE_READER(ActsExamples::RootSimHitReader, mex,
                             "RootSimHitReader", treeName, filePath,
                             outputSimHits);

  ACTS_PYTHON_DECLARE_READER(ActsExamples::CsvExaTrkXGraphReader, mex,
                             "CsvExaTrkXGraphReader", inputDir, inputStem,
                             outputGraph);
}

}  // namespace Acts::Python
