// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>
#include <chrono>
#include <fstream>
#include <iterator>
#include <regex>

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/GridDensityVertexFinder.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Tests/CommonHelpers/DataDirectory.hpp"

namespace Acts {
namespace Test {

using namespace Acts::UnitLiterals;

using Covariance = BoundSymMatrix;
using Propagator = Propagator<EigenStepper<ConstantBField>>;
using Linearizer = HelicalTrackLinearizer<Propagator>;

// Create a test context
GeometryContext geoContext = GeometryContext();
MagneticFieldContext magFieldContext = MagneticFieldContext();

const auto DATAPATH = Acts::Test::getDataPath("vertexing_AMVF_data.csv");

/// @brief Helper struct to store reference vertex related information
struct VertexInfo
{
  // The position
  Vector3D position;
  // The covariance
  ActsSymMatrixD<3> covariance;
  // Number of tracks
  int nTracks;
  // Weight of first track
  double trk1Weight;
  // Vertex compatibility value of first track
  double trk1Comp;
  // Chi2 of first track
  double trk1Chi2;
};

std::tuple<Vertex<BoundParameters>, std::vector<VertexInfo>, std::vector<BoundParameters>> 
readTracksAndVertexCSV(std::string file){

  const std::regex comma(",");

  // Open source file.
  std::ifstream mesh(file);

  // Here we will store the result
  std::vector<std::vector<std::string>> point_coordinates;

  // We want to read all lines of the file
  std::string line{};
  bool isBeamSpot = false;
  bool isTrack = false;
  bool isVertex = false;

  std::shared_ptr<PerigeeSurface> perigeeSurface; 
  std::vector<BoundParameters> tracks;
  std::vector<VertexInfo> vertices;
  Vertex<BoundParameters> beamspotConstraint;

  while (mesh && getline(mesh, line)) {
      // Tokenize line and store result in vector
      std::vector<std::string> row{ std::sregex_token_iterator(line.begin(),line.end(),comma,-1), std::sregex_token_iterator() };
  
      if(row[0] == std::string("beamspot")){
        isBeamSpot = true;
        continue;
      }   
      if(row[0] == "tracks"){
        isTrack = true;
        isBeamSpot = false;
        continue;
      }   
      if(row[0] == "vertices"){
        isTrack = false;
        isBeamSpot = false;
        isVertex = true;
        continue;
      }   

      if(isBeamSpot){
        Vector3D beamspotPos;
        ActsSymMatrixD<3> beamspotCov;
        beamspotPos << std::stod(row[0])* (1_mm), std::stod(row[1])* (1_mm), std::stod(row[2])* (1_mm);
        beamspotCov << std::stod(row[3]), 0, 0, 
                       0, std::stod(row[4]), 0,
                       0, 0, std::stod(row[5]);
        beamspotConstraint.setPosition(beamspotPos);
        beamspotConstraint.setCovariance(beamspotCov);
        perigeeSurface = Surface::makeShared<PerigeeSurface>(beamspotPos);
      } 

      if(isTrack){
        BoundVector params;
        params << std::stod(row[0]), std::stod(row[1]), std::stod(row[2]),
        std::stod(row[3]), std::stod(row[4]) * 1. / (1_MeV), std::stod(row[5]);
        Covariance covMat;
        covMat << std::stod(row[6]), std::stod(row[7]),
        std::stod(row[8]), std::stod(row[9]),
        std::stod(row[10]) * 1. / (1_MeV), std::stod(row[11]), std::stod(row[12]),
        std::stod(row[13]), std::stod(row[14]), std::stod(row[15]),
        std::stod(row[16]) * 1. / (1_MeV), std::stod(row[17]), std::stod(row[18]),
        std::stod(row[19]), std::stod(row[20]), std::stod(row[21]),
        std::stod(row[22]) * 1. / (1_MeV), std::stod(row[23]), std::stod(row[24]),
        std::stod(row[25]), std::stod(row[26]), std::stod(row[27]),
        std::stod(row[28]) * 1. / (1_MeV), std::stod(row[29]),
        std::stod(row[30]) * 1. / (1_MeV),
        std::stod(row[31]) * 1. / (1_MeV),
        std::stod(row[32]) * 1. / (1_MeV),
        std::stod(row[33]) * 1. / (1_MeV),
        std::stod(row[34]) * 1. / (1_MeV * 1_MeV), std::stod(row[35]), std::stod(row[36]), std::stod(row[37]),
        std::stod(row[38]), std::stod(row[39]), std::stod(row[40]), std::stod(row[41]);

        auto boundParams =
        BoundParameters(geoContext, std::move(covMat), params, perigeeSurface);
        tracks.push_back(boundParams);
      }

      if(isVertex){
        Vector3D pos;
        pos << std::stod(row[0]) * (1_mm), std::stod(row[1])* (1_mm), std::stod(row[2])* (1_mm);
        ActsSymMatrixD<3> cov;
        cov << std::stod(row[3]), std::stod(row[4]), std::stod(row[5]), 
                       std::stod(row[6]), std::stod(row[7]), std::stod(row[8]),
                       std::stod(row[9]), std::stod(row[10]), std::stod(row[11]);
        VertexInfo vertexInfo;
        vertexInfo.position = pos;
        vertexInfo.covariance = cov;
        vertexInfo.nTracks = std::stoi(row[12]);
        vertexInfo.trk1Weight = std::stod(row[13]);
        vertexInfo.trk1Comp = std::stod(row[14]);
        vertexInfo.trk1Chi2 = std::stod(row[15]);
        vertices.push_back(vertexInfo);
      }
  }

  return std::make_tuple(beamspotConstraint, vertices, tracks);
}

/// @brief AMVF test with Gaussian seed finder
BOOST_AUTO_TEST_CASE(adaptive_multi_vertex_finder_test) {
  // Set debug mode
  bool debugMode = false;
  // Set up constant B-Field
  ConstantBField bField(Vector3D(0., 0., 2_T));

  // Set up EigenStepper
  // EigenStepper<ConstantBField> stepper(bField);
  EigenStepper<ConstantBField> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  // IP 3D Estimator
  using IPEstimator = ImpactPointEstimator<BoundParameters, Propagator>;

  IPEstimator::Config ipEstimatorCfg(bField, propagator);
  IPEstimator ipEstimator(ipEstimatorCfg);

  std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  AnnealingUtility::Config annealingConfig(temperatures);
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter<BoundParameters, Linearizer>;

  Fitter::Config fitterCfg(ipEstimator);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer for BoundParameters type test
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;

  Fitter fitter(fitterCfg);

  using SeedFinder =
      TrackDensityVertexFinder<Fitter, GaussianTrackDensity<BoundParameters>>;

  SeedFinder seedFinder;

  using Finder = AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEstimator,
                              linearizer);

  // TODO: test this as well!
  // finderConfig.useBeamSpotConstraint = false;

  Finder finder(finderConfig);
  Finder::State state;

  auto csvData = readTracksAndVertexCSV(DATAPATH);
  auto tracks = std::get<2>(csvData);

  if (debugMode) {
    std::cout << "Number of tracks in event: " << tracks.size() << std::endl;
    int maxCout = 10;
    int count = 0;
    for (const auto& trk : tracks) {
      std::cout << count << ". track: " << std::endl;
      std::cout << "params: " << trk << std::endl;
      count++;
      if (count == maxCout) {
        break;
      }
    }
  }

  std::vector<const BoundParameters*> tracksPtr;
  for (const auto& trk : tracks) {
    tracksPtr.push_back(&trk);
  }

  VertexingOptions<BoundParameters> vertexingOptions(geoContext,
                                                     magFieldContext);

  vertexingOptions.vertexConstraint = std::get<0>(csvData);

  auto t1 = std::chrono::system_clock::now();
  auto findResult = finder.find(tracksPtr, vertexingOptions, state);
  auto t2 = std::chrono::system_clock::now();

  auto timediff =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_CHECK(findResult.ok());

  std::vector<Vertex<BoundParameters>> allVertices = *findResult;

  if (debugMode) {
    std::cout << "Time needed: " << timediff << " ms." << std::endl;
    std::cout << "Number of vertices reconstructed: " << allVertices.size()
              << std::endl;

    int count = 0;
    for (const auto& vtx : allVertices) {
      count++;
      std::cout << count << ". Vertex at position: " << vtx.position()[0]
                << ", " << vtx.position()[1] << ", " << vtx.position()[2]
                << std::endl;
      std::cout << count << ". Vertex with cov: " << vtx.covariance()
                << std::endl;
      std::cout << "\t with n tracks: " << vtx.tracks().size() << std::endl;
    }
  }

  // Test expected outcomes from athena implementation
  // Number of reconstructed vertices
  auto verticesInfo = std::get<1>(csvData);
  const int expNRecoVertices = verticesInfo.size();

  BOOST_CHECK_EQUAL(allVertices.size(), expNRecoVertices);

  for(int i = 0; i < expNRecoVertices; i++){
    auto recoVtx = allVertices[i];
    auto expVtx = verticesInfo[i];
    CHECK_CLOSE_ABS(recoVtx.position(), expVtx.position, 0.001_mm);
    CHECK_CLOSE_ABS(recoVtx.covariance(), expVtx.covariance, 0.001_mm);
    BOOST_CHECK_EQUAL(recoVtx.tracks().size(), expVtx.nTracks);
    CHECK_CLOSE_ABS(recoVtx.tracks()[0].trackWeight, expVtx.trk1Weight, 0.003);
    CHECK_CLOSE_ABS(recoVtx.tracks()[0].vertexCompatibility, expVtx.trk1Comp, 0.003);
  }
}

// Dummy user-defined InputTrack type
struct InputTrack {
  InputTrack(const BoundParameters& params, int id)
      : m_parameters(params), m_id(id) {}

  const BoundParameters& parameters() const { return m_parameters; }
  // store e.g. link to original objects here

  int id() const { return m_id; }

 private:
  BoundParameters m_parameters;

  // Some test track ID
  int m_id;
};

/// @brief AMVF test with user-defined input track type
BOOST_AUTO_TEST_CASE(adaptive_multi_vertex_finder_usertype_test) {
  // Set debug mode
  bool debugMode = false;
  // Set up constant B-Field
  ConstantBField bField(Vector3D(0., 0., 2_T));

  // Set up EigenStepper
  // EigenStepper<ConstantBField> stepper(bField);
  EigenStepper<ConstantBField> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  // Create a custom std::function to extract BoundParameters from
  // user-defined InputTrack
  std::function<BoundParameters(InputTrack)> extractParameters =
      [](InputTrack params) { return params.parameters(); };

  // IP 3D Estimator
  using IPEstimator = ImpactPointEstimator<InputTrack, Propagator>;

  IPEstimator::Config ipEstimatorCfg(bField, propagator);
  IPEstimator ipEstimator(ipEstimatorCfg);

  std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  AnnealingUtility::Config annealingConfig(temperatures);
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter<InputTrack, Linearizer>;

  Fitter::Config fitterCfg(ipEstimator);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;

  Fitter fitter(fitterCfg, extractParameters);

  using SeedFinder =
      TrackDensityVertexFinder<Fitter, GaussianTrackDensity<InputTrack>>;

  SeedFinder seedFinder(extractParameters);

  using Finder = AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEstimator,
                              linearizer);
  Finder::State state;

  Finder finder(finderConfig, extractParameters);

  auto csvData = readTracksAndVertexCSV(DATAPATH);
  auto tracks = std::get<2>(csvData);

  std::vector<InputTrack> userTracks;
  int idCount = 0;
  for (const auto& trk : tracks) {
    userTracks.push_back(InputTrack(trk, idCount));
    idCount++;
  }

  if (debugMode) {
    std::cout << "Number of tracks in event: " << tracks.size() << std::endl;
    int maxCout = 10;
    int count = 0;
    for (const auto& trk : tracks) {
      std::cout << count << ". track: " << std::endl;
      std::cout << "params: " << trk << std::endl;
      count++;
      if (count == maxCout) {
        break;
      }
    }
  }

  std::vector<const InputTrack*> userTracksPtr;
  for (const auto& trk : userTracks) {
    userTracksPtr.push_back(&trk);
  }

  VertexingOptions<InputTrack> vertexingOptions(geoContext, magFieldContext);

  Vertex<InputTrack> constraintVtx;
  constraintVtx.setPosition(std::get<0>(csvData).position());
  constraintVtx.setCovariance(std::get<0>(csvData).covariance());

  vertexingOptions.vertexConstraint = constraintVtx;

  auto findResult = finder.find(userTracksPtr, vertexingOptions, state);

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_CHECK(findResult.ok());

  std::vector<Vertex<InputTrack>> allVertices = *findResult;

  if (debugMode) {
    std::cout << "Number of vertices reconstructed: " << allVertices.size()
              << std::endl;

    int count = 0;
    for (const auto& vtx : allVertices) {
      count++;
      std::cout << count << ". Vertex at position: " << vtx.position()[0]
                << ", " << vtx.position()[1] << ", " << vtx.position()[2]
                << std::endl;
      std::cout << count << ". Vertex with cov: " << vtx.covariance()
                << std::endl;
      std::cout << "\t with n tracks: " << vtx.tracks().size() << std::endl;
    }
    for (auto& trk : allVertices[0].tracks()) {
      std::cout << "Track ID at first vertex: " << trk.originalParams->id()
                << std::endl;
    }
  }

  auto verticesInfo = std::get<1>(csvData);
  const int expNRecoVertices = verticesInfo.size();

  BOOST_CHECK_EQUAL(allVertices.size(), expNRecoVertices);

  for(int i = 0; i < expNRecoVertices; i++){
    auto recoVtx = allVertices[i];
    auto expVtx = verticesInfo[i];
    CHECK_CLOSE_ABS(recoVtx.position(), expVtx.position, 0.001_mm);
    CHECK_CLOSE_ABS(recoVtx.covariance(), expVtx.covariance, 0.001_mm);
    BOOST_CHECK_EQUAL(recoVtx.tracks().size(), expVtx.nTracks);
    CHECK_CLOSE_ABS(recoVtx.tracks()[0].trackWeight, expVtx.trk1Weight, 0.003);
    CHECK_CLOSE_ABS(recoVtx.tracks()[0].vertexCompatibility, expVtx.trk1Comp, 0.003);
  }
}

/// @brief AMVF test with grid seed finder
BOOST_AUTO_TEST_CASE(adaptive_multi_vertex_finder_grid_seed_finder_test) {
  // Set debug mode
  bool debugMode = false;
  if (debugMode) {
    std::cout << "Starting AMVF test with grid seed finder..." << std::endl;
  }
  // Set up constant B-Field
  ConstantBField bField(Vector3D(0., 0., 2_T));

  // Set up EigenStepper
  // EigenStepper<ConstantBField> stepper(bField);
  EigenStepper<ConstantBField> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  // IP Estimator
  using IPEstimator = ImpactPointEstimator<BoundParameters, Propagator>;

  IPEstimator::Config ipEstCfg(bField, propagator);
  IPEstimator ipEst(ipEstCfg);

  std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  AnnealingUtility::Config annealingConfig(temperatures);
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter<BoundParameters, Linearizer>;

  Fitter::Config fitterCfg(ipEst);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer for BoundParameters type test
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;

  Fitter fitter(fitterCfg);

  using SeedFinder = GridDensityVertexFinder<4000, 55>;
  SeedFinder::Config seedFinderCfg(250);
  seedFinderCfg.cacheGridStateForTrackRemoval = true;

  SeedFinder seedFinder(seedFinderCfg);

  using Finder = AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEst, linearizer);

  // TODO: test this as well!
  // finderConfig.useBeamSpotConstraint = false;

  Finder finder(finderConfig);
  Finder::State state;

  auto csvData = readTracksAndVertexCSV(DATAPATH);
  auto tracks = std::get<2>(csvData);

  if (debugMode) {
    std::cout << "Number of tracks in event: " << tracks.size() << std::endl;
    int maxCout = 10;
    int count = 0;
    for (const auto& trk : tracks) {
      std::cout << count << ". track: " << std::endl;
      std::cout << "params: " << trk << std::endl;
      count++;
      if (count == maxCout) {
        break;
      }
    }
  }

  std::vector<const BoundParameters*> tracksPtr;
  for (const auto& trk : tracks) {
    tracksPtr.push_back(&trk);
  }

  VertexingOptions<BoundParameters> vertexingOptions(geoContext,
                                                     magFieldContext);

  vertexingOptions.vertexConstraint = std::get<0>(csvData);

  auto t1 = std::chrono::system_clock::now();
  auto findResult = finder.find(tracksPtr, vertexingOptions, state);
  auto t2 = std::chrono::system_clock::now();

  auto timediff =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_CHECK(findResult.ok());

  std::vector<Vertex<BoundParameters>> allVertices = *findResult;

  if (debugMode) {
    std::cout << "Time needed: " << timediff << " ms." << std::endl;
    std::cout << "Number of vertices reconstructed: " << allVertices.size()
              << std::endl;

    int count = 0;
    for (const auto& vtx : allVertices) {
      count++;
      std::cout << count << ". Vertex at position: " << vtx.position()[0]
                << ", " << vtx.position()[1] << ", " << vtx.position()[2]
                << std::endl;
      std::cout << count << ". Vertex with cov: " << vtx.covariance()
                << std::endl;
      std::cout << "\t with n tracks: " << vtx.tracks().size() << std::endl;
    }
  }
  // Test expected outcomes from athena implementation
  // Number of reconstructed vertices
  auto verticesInfo = std::get<1>(csvData);
  const int expNRecoVertices = verticesInfo.size();

  BOOST_CHECK_EQUAL(allVertices.size(), expNRecoVertices);
  std::vector<bool> vtxFound(expNRecoVertices, false);

  for (auto vtx : allVertices) {
    double vtxZ = vtx.position()[2];
    double diffZ = 1e5;
    int foundVtxIdx = -1;
    for (int i = 0; i < expNRecoVertices; i++) {

      if (not vtxFound[i]) {
        if (std::abs(vtxZ - verticesInfo[i].position[2]) < diffZ) {
          diffZ = std::abs(vtxZ - verticesInfo[i].position[2]);
          foundVtxIdx = i;
        }
      }
    }
    if (diffZ < 0.5_mm) {
      vtxFound[foundVtxIdx] = true;
      CHECK_CLOSE_ABS(vtx.tracks().size(), verticesInfo[foundVtxIdx].nTracks, 1);
    }
  }
  for (bool found : vtxFound) {
    BOOST_CHECK_EQUAL(found, true);
  }
}

}  // namespace Test
}  // namespace Acts
