// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsAlignment/Kernel/detail/AlignmentEngine.hpp"

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/KalmanGlobalCovariance.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"

#include <cmath>
#include <random>
#include <string>

using namespace Acts;
using namespace ActsAlignment;

namespace {
using SourceLink = MinimalSourceLink;
using Covariance = BoundSymMatrix;

template <BoundIndices... params>
using MeasurementType = Measurement<SourceLink, BoundIndices, params...>;

std::normal_distribution<double> gauss(0., 1.);
std::default_random_engine generator(42);

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();
CalibrationContext calContext = CalibrationContext();

///
/// @brief Contruct a telescope-like detector
///
struct TelescopeTrackingGeometry {
  /// Default constructor for the Cubit tracking geometry
  ///
  /// @param gctx the geometry context for this geometry at building time
  TelescopeTrackingGeometry(std::reference_wrapper<const GeometryContext> gctx)
      : geoContext(gctx) {
    using namespace UnitLiterals;

    // Construct the rotation
    double rotationAngle = 90_degree;
    Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
    Vector3D yPos(0., 1., 0.);
    Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
    rotation.col(0) = xPos;
    rotation.col(1) = yPos;
    rotation.col(2) = zPos;

    // Boundaries of the surfaces
    rBounds =
        std::make_shared<const RectangleBounds>(RectangleBounds(0.1_m, 0.1_m));

    // Material of the surfaces
    MaterialSlab matProp(Acts::Test::makeSilicon(), 0.5_mm);
    surfaceMaterial = std::make_shared<HomogeneousSurfaceMaterial>(matProp);
  }

  ///
  /// Call operator to build the standard cubic tracking geometry
  ///
  std::shared_ptr<const TrackingGeometry> operator()() {
    using namespace UnitLiterals;

    // Set translation vectors
    std::vector<Vector3D> translations;
    translations.reserve(6);
    translations.push_back({-500_mm, 0., 0.});
    translations.push_back({-300_mm, 0., 0.});
    translations.push_back({-100_mm, 0., 0.});
    translations.push_back({100_mm, 0., 0.});
    translations.push_back({300_mm, 0., 0.});
    translations.push_back({500_mm, 0., 0.});

    // Construct layer configs
    std::vector<CuboidVolumeBuilder::LayerConfig> lConfs;
    lConfs.reserve(6);
    unsigned int i;
    for (i = 0; i < translations.size(); i++) {
      CuboidVolumeBuilder::SurfaceConfig sConf;
      sConf.position = translations[i];
      sConf.rotation = rotation;
      sConf.rBounds = rBounds;
      sConf.surMat = surfaceMaterial;
      // The thickness to construct the associated detector element
      sConf.thickness = 1._um;
      sConf.detElementConstructor =
          [](const Transform3D& trans,
             std::shared_ptr<const RectangleBounds> bounds, double thickness) {
            return new Acts::Test::DetectorElementStub(trans, bounds,
                                                       thickness);
          };
      CuboidVolumeBuilder::LayerConfig lConf;
      lConf.surfaceCfg = sConf;
      lConfs.push_back(lConf);
    }

    // Construct volume config
    CuboidVolumeBuilder::VolumeConfig vConf;
    vConf.position = {0., 0., 0.};
    vConf.length = {1.2_m, 1._m, 1._m};
    vConf.layerCfg = lConfs;
    vConf.name = "Tracker";

    // Construct volume builder config
    CuboidVolumeBuilder::Config conf;
    conf.position = {0., 0., 0.};
    conf.length = {1.2_m, 1._m, 1._m};
    conf.volumeCfg = {vConf};  // one volume

    // Build detector
    CuboidVolumeBuilder cvb;
    cvb.setConfig(conf);
    TrackingGeometryBuilder::Config tgbCfg;
    tgbCfg.trackingVolumeBuilders.push_back(
        [=](const auto& context, const auto& inner, const auto& vb) {
          return cvb.trackingVolume(context, inner, vb);
        });
    TrackingGeometryBuilder tgb(tgbCfg);
    std::shared_ptr<const TrackingGeometry> detector =
        tgb.trackingGeometry(tgContext);

    // Build and return tracking geometry
    return detector;
  }

  RotationMatrix3D rotation = RotationMatrix3D::Identity();
  std::shared_ptr<const RectangleBounds> rBounds = nullptr;
  std::shared_ptr<const ISurfaceMaterial> surfaceMaterial = nullptr;

  std::reference_wrapper<const GeometryContext> geoContext;
};

struct MeasurementCreator {
  /// @brief Constructor
  MeasurementCreator() = default;

  /// The detector resolution
  std::array<double, 2> resolution = {30_um, 50_um};

  struct this_result {
    // The measurements
    std::vector<FittableMeasurement<SourceLink>> measurements;
  };

  using result_type = this_result;

  /// @brief Operater that is callable by an ActionList. The function collects
  /// the surfaces
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagator
  /// @param [out] result Vector of matching surfaces
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    // monitor the current surface
    auto surface = state.navigation.currentSurface;
    if (surface and surface->associatedDetectorElement()) {
      // Apply global to local
      Acts::Vector2D lPos =
          surface
              ->globalToLocal(state.geoContext,
                              stepper.position(state.stepping),
                              stepper.direction(state.stepping))
              .value();
      // 2D measurement
      double dx = resolution[eBoundLoc0] * gauss(generator);
      double dy = resolution[eBoundLoc1] * gauss(generator);
      // Measurment covariance
      ActsSymMatrixD<2> cov2D;
      cov2D << resolution[eBoundLoc0] * resolution[eBoundLoc0], 0., 0.,
          resolution[eBoundLoc1] * resolution[eBoundLoc1];
      // Create a measurement
      MeasurementType<eBoundLoc0, eBoundLoc1> m01(surface->getSharedPtr(), {},
                                                  cov2D, lPos[eBoundLoc0] + dx,
                                                  lPos[eBoundLoc1] + dy);
      result.measurements.push_back(std::move(m01));
    }
  }
};

struct KalmanFitterInputTrajectory {
  // The measurements
  std::vector<FittableMeasurement<SourceLink>> measurements;
  // The source links
  std::vector<SourceLink> sourcelinks;
  // The start parameters
  std::optional<CurvilinearTrackParameters> startParameters;
};

///
/// Function to create trajectories for kalman fitter
///
std::vector<KalmanFitterInputTrajectory> createTrajectories(
    const std::shared_ptr<const TrackingGeometry>& detector,
    size_t nTrajectories,
    const std::array<double, 2>& localSigma = {1000_um, 1000_um},
    const double& pSigma = 0.025_GeV) {
  // Build navigator for the measurement creatoin
  Navigator mNavigator(detector);
  mNavigator.resolvePassive = false;
  mNavigator.resolveMaterial = true;
  mNavigator.resolveSensitive = true;

  // Use straingt line stepper to create the measurements
  StraightLineStepper mStepper;

  // Define the measurement propagator
  using MeasurementPropagator = Propagator<StraightLineStepper, Navigator>;

  // Build propagator for the measurement creation
  MeasurementPropagator mPropagator(mStepper, mNavigator);

  // Create action list for the measurement creation
  using MeasurementActions = ActionList<MeasurementCreator>;
  using MeasurementAborters = AbortList<EndOfWorldReached>;

  // Set options for propagator
  PropagatorOptions<MeasurementActions, MeasurementAborters> mOptions(
      tgContext, mfContext, getDummyLogger());

  std::vector<KalmanFitterInputTrajectory> trajectories;
  trajectories.reserve(nTrajectories);
  for (unsigned int iTrack = 0; iTrack < nTrajectories; iTrack++) {
    if (iTrack % 10 == 0) {
      std::cout << "Processing track: " << iTrack << "..." << std::endl;
    }
    // Set initial parameters for the particle track
    Vector4D mPos4(-1_m, 100_um * gauss(generator), 100_um * gauss(generator),
                   42_ns);
    Vector3D mDir(1_GeV, 0.01_GeV * gauss(generator),
                  0.01_GeV * gauss(generator));

    CurvilinearTrackParameters mStart(mPos4, mDir, 1_e / 1_GeV);
    // Launch and collect - the measurements
    auto mResult = mPropagator.propagate(mStart, mOptions);
    if (not mResult.ok()) {
      continue;
    }

    // Extract measurements from result of propagation.
    KalmanFitterInputTrajectory traj;
    traj.measurements =
        std::move((*mResult)
                      .template get<MeasurementCreator::result_type>()
                      .measurements);

    // Make a vector of source links as input to the KF
    std::transform(traj.measurements.begin(), traj.measurements.end(),
                   std::back_inserter(traj.sourcelinks),
                   [](const auto& m) { return SourceLink{&m}; });

    // Smear the start parameters to be used as input of KF
    Covariance cov;
    cov << std::pow(localSigma[0], 2), 0., 0., 0., 0., 0., 0.,
        std::pow(localSigma[1], 2), 0., 0., 0., 0., 0., 0., pSigma, 0., 0., 0.,
        0., 0., 0., pSigma, 0., 0., 0., 0., 0., 0., 0.01, 0., 0., 0., 0., 0.,
        0., 1.;
    Vector4D rPos4(mPos4.x(), mPos4.y() + localSigma[0] * gauss(generator),
                   mPos4.z() + localSigma[1] * gauss(generator), 42_ns);
    Vector3D rDir(mDir.x(), mDir.y() + pSigma * gauss(generator),
                  mDir.z() + pSigma * gauss(generator));
    CurvilinearTrackParameters rStart(rPos4, rDir, 1_e / 1_GeV, cov);
    traj.startParameters = rStart;

    trajectories.push_back(std::move(traj));
  }
  return trajectories;
}
}  // namespace

///
/// @brief Unit test for KF-based alignment algorithm
///
BOOST_AUTO_TEST_CASE(Alignment_zero_field) {
  // Build detector
  TelescopeTrackingGeometry tGeometry(tgContext);
  auto detector = tGeometry();

  // Create the trajectories
  const auto& trajectories = createTrajectories(detector, 100);

  // The KalmanFitter - we use the eigen stepper for covariance transport
  Navigator rNavigator(detector);
  rNavigator.resolvePassive = false;
  rNavigator.resolveMaterial = true;
  rNavigator.resolveSensitive = true;

  // Configure propagation with deactivated B-field
  ConstantBField bField(Vector3D(0., 0., 0.));
  using RecoStepper = EigenStepper<ConstantBField>;
  RecoStepper rStepper(bField);
  using RecoPropagator = Propagator<RecoStepper, Navigator>;
  RecoPropagator rPropagator(rStepper, rNavigator);

  using Updater = GainMatrixUpdater;
  using Smoother = GainMatrixSmoother;
  using KalmanFitter = KalmanFitter<RecoPropagator, Updater, Smoother>;
  using KalmanFitterOptions =
      KalmanFitterOptions<MinimalSourceLinkCalibrator, VoidOutlierFinder>;

  auto kfLogger = getDefaultLogger("KalmanFilter", Logging::VERBOSE);
  // Construct the KalmanFitter options
  KalmanFitterOptions kfOptions(
      tgContext, mfContext, calContext, MinimalSourceLinkCalibrator(),
      VoidOutlierFinder(), LoggerWrapper{*kfLogger}, PropagatorPlainOptions());

  // Construct the KalmanFitter
  KalmanFitter kFitter(rPropagator);

  // Construct the alignment algorithm
  Alignment alignment(kFitter, getDefaultLogger("Alignment", Logging::VERBOSE));

  // Construct an non-updating alignment updater
  AlignedTransformUpdater voidAlignUpdater =
      [](Acts::DetectorElementBase* /*unused*/,
         const Acts::GeometryContext& /*unused*/,
         const Acts::Transform3D& /*unused*/) { return true; };

  // Construct the alignment options
  AlignmentOptions<KalmanFitterOptions> alignOptions(kfOptions,
                                                     voidAlignUpdater);
  alignOptions.maxIterations = 1;

  // The surfaces to be aligned
  unsigned int iSurface = 0;
  std::unordered_map<const Surface*, size_t> idxedAlignSurfaces;
  detector->visitSurfaces([&](const Surface* surface) {
    // Missing out the forth layer
    if (surface and surface->associatedDetectorElement() and
        surface->geometryId().layer() != 8) {
      alignOptions.alignedDetElements.push_back(
          const_cast<DetectorElementBase*>(
              surface->associatedDetectorElement()));
      idxedAlignSurfaces.emplace(surface, iSurface);
      iSurface++;
    }
  });

  // The alignment mask
  const auto alignMask =
      std::bitset<Acts::eAlignmentSize>(std::string("111111"));

  // Test the method to evaluate alignment state for a single track
  const auto& inputTraj = trajectories.front();
  kfOptions.referenceSurface = &(*inputTraj.startParameters).referenceSurface();
  auto evaluateRes = alignment.evaluateTrackAlignmentState(
      kfOptions.geoContext, inputTraj.sourcelinks, *inputTraj.startParameters,
      kfOptions, idxedAlignSurfaces, alignMask);
  BOOST_CHECK(evaluateRes.ok());
  const auto& alignState = evaluateRes.value();
  std::cout << "Chi2/dof = " << alignState.chi2 / alignState.alignmentDof
            << std::endl;

  // Check the dimensions
  BOOST_CHECK_EQUAL(alignState.measurementDim, 12);
  BOOST_CHECK_EQUAL(alignState.trackParametersDim, 36);
  // Check the alignment dof
  BOOST_CHECK_EQUAL(alignState.alignmentDof, 30);
  BOOST_CHECK_EQUAL(alignState.alignedSurfaces.size(), 5);
  // Check the measurements covariance
  BOOST_CHECK_EQUAL(alignState.measurementCovariance.rows(), 12);
  const ActsSymMatrixD<2> measCov =
      alignState.measurementCovariance.block<2, 2>(2, 2);
  ActsSymMatrixD<2> cov2D;
  cov2D << 30_um * 30_um, 0, 0, 50_um * 50_um;
  CHECK_CLOSE_ABS(measCov, cov2D, 1e-10);
  // Check the track parameters covariance matrix. Its rows/columns scales
  // with the number of measurement states
  BOOST_CHECK_EQUAL(alignState.trackParametersCovariance.rows(), 36);
  // Check the projection matrix
  BOOST_CHECK_EQUAL(alignState.projectionMatrix.rows(), 12);
  BOOST_CHECK_EQUAL(alignState.projectionMatrix.cols(), 36);
  const ActsMatrixD<2, 6> proj = alignState.projectionMatrix.block<2, 6>(0, 0);
  const ActsMatrixD<2, 6> refProj = ActsMatrixD<2, 6>::Identity();
  CHECK_CLOSE_ABS(proj, refProj, 1e-10);
  // Check the residual
  BOOST_CHECK_EQUAL(alignState.residual.size(), 12);
  // Check the residual covariance
  BOOST_CHECK_EQUAL(alignState.residualCovariance.rows(), 12);
  // Check the alignment to residual derivative
  BOOST_CHECK_EQUAL(alignState.alignmentToResidualDerivative.rows(), 12);
  BOOST_CHECK_EQUAL(alignState.alignmentToResidualDerivative.cols(), 30);
  // Check the chi2 derivative
  BOOST_CHECK_EQUAL(alignState.alignmentToChi2Derivative.size(), 30);
  BOOST_CHECK_EQUAL(alignState.alignmentToChi2SecondDerivative.rows(), 30);

  // Test the align method
  std::vector<std::vector<SourceLink>> trajCollection;
  trajCollection.reserve(100);
  std::vector<CurvilinearTrackParameters> sParametersCollection;
  sParametersCollection.reserve(100);
  for (const auto& traj : trajectories) {
    trajCollection.push_back(traj.sourcelinks);
    sParametersCollection.push_back(*traj.startParameters);
  }
  auto alignRes =
      alignment.align(trajCollection, sParametersCollection, alignOptions);

  // BOOST_CHECK(alignRes.ok());
}
