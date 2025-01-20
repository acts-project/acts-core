// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackContainerFrontendConcept.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackProxyConcept.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/TrackFitting/GlobalChiSquareFitterError.hpp"
#include "Acts/TrackFitting/detail/VoidFitterComponents.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"

#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <type_traits>
#include <unordered_map>

namespace Acts::Experimental {

namespace Gx2fConstants {
constexpr std::string_view gx2fnUpdateColumn = "Gx2fnUpdateColumn";

// Mask for the track states. We don't need Predicted and Filtered
constexpr TrackStatePropMask trackStateMask = TrackStatePropMask::Smoothed |
                                              TrackStatePropMask::Jacobian |
                                              TrackStatePropMask::Calibrated;

// A projector used for scattering. By using Jacobian * phiThetaProjector one
// gets only the derivatives for the variables phi and theta.
const Eigen::Matrix<double, eBoundSize, 2> phiThetaProjector = [] {
  Eigen::Matrix<double, eBoundSize, 2> m =
      Eigen::Matrix<double, eBoundSize, 2>::Zero();
  m(eBoundPhi, 0) = 1.0;
  m(eBoundTheta, 1) = 1.0;
  return m;
}();
}  // namespace Gx2fConstants

/// Extension struct which holds delegates to customise the GX2F behaviour
template <typename traj_t>
struct Gx2FitterExtensions {
  using TrackStateProxy = typename MultiTrajectory<traj_t>::TrackStateProxy;
  using ConstTrackStateProxy =
      typename MultiTrajectory<traj_t>::ConstTrackStateProxy;
  using Parameters = typename TrackStateProxy::Parameters;

  using Calibrator =
      Delegate<void(const GeometryContext&, const CalibrationContext&,
                    const SourceLink&, TrackStateProxy)>;

  using Updater = Delegate<Result<void>(const GeometryContext&, TrackStateProxy,
                                        const Logger&)>;

  using OutlierFinder = Delegate<bool(ConstTrackStateProxy)>;

  /// The Calibrator is a dedicated calibration algorithm that allows
  /// to calibrate measurements using track information, this could be
  /// e.g. sagging for wires, module deformations, etc.
  Calibrator calibrator;

  /// The updater incorporates measurement information into the track parameters
  Updater updater;

  /// Determines whether a measurement is supposed to be considered as an
  /// outlier
  OutlierFinder outlierFinder;

  /// Retrieves the associated surface from a source link
  SourceLinkSurfaceAccessor surfaceAccessor;

  /// Default constructor which connects the default void components
  Gx2FitterExtensions() {
    calibrator.template connect<&detail::voidFitterCalibrator<traj_t>>();
    updater.template connect<&detail::voidFitterUpdater<traj_t>>();
    outlierFinder.template connect<&detail::voidOutlierFinder<traj_t>>();
    surfaceAccessor.connect<&detail::voidSurfaceAccessor>();
  }
};

/// Combined options for the Global-Chi-Square fitter.
///
/// @tparam traj_t The trajectory type
template <typename traj_t>
struct Gx2FitterOptions {
  /// PropagatorOptions with context.
  ///
  /// @param gctx The geometry context for this fit
  /// @param mctx The magnetic context for this fit
  /// @param cctx The calibration context for this fit
  /// @param extensions_ The KF extensions
  /// @param pOptions The plain propagator options
  /// @param rSurface The reference surface for the fit to be expressed at
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  /// @param freeToBoundCorrection_ Correction for non-linearity effect during transform from free to bound
  /// @param nUpdateMax_ Max number of iterations for updating the parameters
  /// @param relChi2changeCutOff_ Check for convergence (abort condition). Set to 0 to skip.
  Gx2FitterOptions(const GeometryContext& gctx,
                   const MagneticFieldContext& mctx,
                   std::reference_wrapper<const CalibrationContext> cctx,
                   Gx2FitterExtensions<traj_t> extensions_,
                   const PropagatorPlainOptions& pOptions,
                   const Surface* rSurface = nullptr, bool mScattering = false,
                   bool eLoss = false,
                   const FreeToBoundCorrection& freeToBoundCorrection_ =
                       FreeToBoundCorrection(false),
                   const std::size_t nUpdateMax_ = 5,
                   double relChi2changeCutOff_ = 1e-5)
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        extensions(extensions_),
        propagatorPlainOptions(pOptions),
        referenceSurface(rSurface),
        multipleScattering(mScattering),
        energyLoss(eLoss),
        freeToBoundCorrection(freeToBoundCorrection_),
        nUpdateMax(nUpdateMax_),
        relChi2changeCutOff(relChi2changeCutOff_) {}

  /// Contexts are required and the options must not be default-constructible.
  Gx2FitterOptions() = delete;

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  Gx2FitterExtensions<traj_t> extensions;

  /// The trivial propagator options
  PropagatorPlainOptions propagatorPlainOptions;

  /// The reference Surface
  const Surface* referenceSurface = nullptr;

  /// Whether to consider multiple scattering
  bool multipleScattering = false;

  /// Whether to consider energy loss
  bool energyLoss = false;

  /// Whether to include non-linear correction during global to local
  /// transformation
  FreeToBoundCorrection freeToBoundCorrection;

  /// Max number of iterations during the fit (abort condition)
  std::size_t nUpdateMax = 5;

  /// Check for convergence (abort condition). Set to 0 to skip.
  double relChi2changeCutOff = 1e-7;
};

template <typename traj_t>
struct Gx2FitterResult {
  // Fitted states that the actor has handled.
  traj_t* fittedStates{nullptr};

  // This is the index of the 'tip' of the track stored in multitrajectory.
  // This corresponds to the last measurement state in the multitrajectory.
  // Since this KF only stores one trajectory, it is unambiguous.
  // Acts::MultiTrajectoryTraits::kInvalid is the start of a trajectory.
  std::size_t lastMeasurementIndex = Acts::MultiTrajectoryTraits::kInvalid;

  // This is the index of the 'tip' of the states stored in multitrajectory.
  // This corresponds to the last state in the multitrajectory.
  // Since this KF only stores one trajectory, it is unambiguous.
  // Acts::MultiTrajectoryTraits::kInvalid is the start of a trajectory.
  std::size_t lastTrackIndex = Acts::MultiTrajectoryTraits::kInvalid;

  // The optional Parameters at the provided surface
  std::optional<BoundTrackParameters> fittedParameters;

  // Counter for states with non-outlier measurements
  std::size_t measurementStates = 0;

  // Counter for measurements holes
  // A hole correspond to a surface with an associated detector element with no
  // associated measurement. Holes are only taken into account if they are
  // between the first and last measurements.
  std::size_t measurementHoles = 0;

  // Counter for handled states
  std::size_t processedStates = 0;

  // Counter for handled measurements
  std::size_t processedMeasurements = 0;

  // Indicator if track fitting has been done
  bool finished = false;

  // Measurement surfaces without hits
  std::vector<const Surface*> missedActiveSurfaces;

  // Measurement surfaces handled in both forward and
  // backward filtering
  std::vector<const Surface*> passedAgainSurfaces;

  Result<void> result{Result<void>::success()};

  // Count how many surfaces have been hit
  std::size_t surfaceCount = 0;
};

/// @brief A container to store scattering properties for each material surface
///
/// This struct holds the scattering angles, the inverse covariance of the
/// material, and a validity flag indicating whether the material is valid for
/// the scattering process.
struct ScatteringProperties {
 public:
  /// @brief Constructor to initialize scattering properties.
  ///
  /// @param scatteringAngles_ The vector of scattering angles.
  /// @param invCovarianceMaterial_ The inverse covariance of the material.
  /// @param materialIsValid_ A boolean flag indicating whether the material is valid.
  ScatteringProperties(const BoundVector& scatteringAngles_,
                       const double invCovarianceMaterial_,
                       const bool materialIsValid_)
      : m_scatteringAngles(scatteringAngles_),
        m_invCovarianceMaterial(invCovarianceMaterial_),
        m_materialIsValid(materialIsValid_) {}

  // Accessor for the scattering angles.
  const BoundVector& scatteringAngles() const { return m_scatteringAngles; }

  // Accessor for a modifiable reference to the scattering angles
  BoundVector& scatteringAngles() { return m_scatteringAngles; }

  // Accessor for the inverse covariance of the material.
  double invCovarianceMaterial() const { return m_invCovarianceMaterial; }

  // Accessor for the material validity flag.
  bool materialIsValid() const { return m_materialIsValid; }

 private:
  /// Vector of scattering angles. The vector is usually all zeros except for
  /// eBoundPhi and eBoundTheta.
  BoundVector m_scatteringAngles;

  /// Inverse covariance of the material. Compute with e.g. the Highland
  /// formula.
  double m_invCovarianceMaterial;

  /// Flag indicating whether the material is valid. Commonly vacuum and zero
  /// thickness material will be ignored.
  bool m_materialIsValid;
};

/// @brief A container to manage all properties of a gx2f system
///
/// This struct manages the mathematical infrastructure for the gx2f. It
/// initializes and maintains the extended aMatrix and extended bVector.
struct Gx2fSystem {
 public:
  /// @brief Constructor to initialize matrices and vectors to zero based on specified dimensions.
  ///
  /// @param nDims Number of dimensions for the extended matrix and vector.
  explicit Gx2fSystem(std::size_t nDims)
      : m_nDims{nDims},
        m_aMatrix{Eigen::MatrixXd::Zero(nDims, nDims)},
        m_bVector{Eigen::VectorXd::Zero(nDims)} {}

  // Accessor for nDims (const reference).
  std::size_t nDims() const { return m_nDims; }

  // Accessor for chi2
  double chi2() const { return m_chi2; }

  // Modifier for chi2
  double& chi2() { return m_chi2; }

  // Accessor for the matrix.
  const Eigen::MatrixXd& aMatrix() const { return m_aMatrix; }

  // Accessor for a modifiable reference to the matrix.
  Eigen::MatrixXd& aMatrix() { return m_aMatrix; }

  // Accessor for the vector.
  const Eigen::VectorXd& bVector() const { return m_bVector; }

  // Accessor for a modifiable reference to the vector.
  Eigen::VectorXd& bVector() { return m_bVector; }

  // Accessor for NDF
  std::size_t ndf() const { return m_ndf; }

  // Modifier for NDF
  std::size_t& ndf() { return m_ndf; }

  //  It automatically deduces if we want to fit e.g. q/p and adjusts itself
  //  later. We have only 3 cases, because we always have l0, l1, phi, theta:
  // - 4: no magnetic field -> q/p is empty
  // - 5: no time measurement -> time is not fittable
  // - 6: full fit
  std::size_t findRequiredNdf() {
    std::size_t ndfSystem = 0;
    if (m_aMatrix(4, 4) == 0) {
      ndfSystem = 4;
    } else if (m_aMatrix(5, 5) == 0) {
      ndfSystem = 5;
    } else {
      ndfSystem = 6;
    }

    return ndfSystem;
  }

  bool isWellDefined() { return m_ndf > findRequiredNdf(); }

 private:
  /// Number of dimensions of the (extended) system
  std::size_t m_nDims;

  /// Sum of chi-squared values.
  double m_chi2 = 0.;

  /// Extended matrix for accumulation.
  Eigen::MatrixXd m_aMatrix;

  /// Extended vector for accumulation.
  Eigen::VectorXd m_bVector;

  /// Number of degrees of freedom of the system
  std::size_t m_ndf = 0u;
};

/// @brief Adds a measurement to the GX2F equation system in a modular backend function.
///
/// This function processes measurement data and integrates it into the GX2F
/// system.
///
/// @param extendedSystem All parameters of the current equation system to update.
/// @param jacobianFromStart The Jacobian matrix from the start to the current state.
/// @param covarianceMeasurement The covariance matrix of the measurement.
/// @param predicted The predicted state vector based on the track state.
/// @param measurement The measurement vector.
/// @param projector The projection matrix.
/// @param logger A logger instance.
///
/// @note The dynamic Eigen matrices are suboptimal. We could think of
/// templating again in the future on kMeasDims. We currently use dynamic
/// matrices to reduce the memory during compile time.
void addMeasurementToGx2fSumsBackend(
    Gx2fSystem& extendedSystem,
    const std::vector<BoundMatrix>& jacobianFromStart,
    const Eigen::MatrixXd& covarianceMeasurement, const BoundVector& predicted,
    const Eigen::VectorXd& measurement, const Eigen::MatrixXd& projector,
    const Logger& logger);

/// @brief Process measurements and fill the aMatrix and bVector
///
/// The function processes each measurement for the GX2F Actor fitting process.
/// It extracts the information from the track state and adds it to aMatrix,
/// bVector, and chi2sum.
///
/// @tparam kMeasDim Number of dimensions of the measurement
/// @tparam track_state_t The type of the track state
///
/// @param extendedSystem All parameters of the current equation system to update
/// @param jacobianFromStart The Jacobian matrix from start to the current state
/// @param trackState The track state to analyse
/// @param logger A logger instance
template <std::size_t kMeasDim, typename track_state_t>
void addMeasurementToGx2fSums(Gx2fSystem& extendedSystem,
                              const std::vector<BoundMatrix>& jacobianFromStart,
                              const track_state_t& trackState,
                              const Logger& logger) {
  const ActsSquareMatrix<kMeasDim> covarianceMeasurement =
      trackState.template calibratedCovariance<kMeasDim>();

  const BoundVector predicted = trackState.smoothed();

  const ActsVector<kMeasDim> measurement =
      trackState.template calibrated<kMeasDim>();

  const ActsMatrix<kMeasDim, eBoundSize> projector =
      trackState.template projectorSubspaceHelper<kMeasDim>().projector();

  addMeasurementToGx2fSumsBackend(extendedSystem, jacobianFromStart,
                                  covarianceMeasurement, predicted, measurement,
                                  projector, logger);
}

/// @brief Process material and fill the aMatrix and bVector
///
/// The function processes each material for the GX2F Actor fitting process.
/// It extracts the information from the track state and adds it to aMatrix,
/// bVector, and chi2sum.
///
/// @tparam track_state_t The type of the track state
///
/// @param extendedSystem All parameters of the current equation system
/// @param nMaterialsHandled How many materials we already handled. Used for the offset.
/// @param scatteringMap The scattering map, containing all scattering angles and covariances
/// @param trackState The track state to analyse
/// @param logger A logger instance
template <typename track_state_t>
void addMaterialToGx2fSums(
    Gx2fSystem& extendedSystem, const std::size_t nMaterialsHandled,
    const std::unordered_map<GeometryIdentifier, ScatteringProperties>&
        scatteringMap,
    const track_state_t& trackState, const Logger& logger) {
  // Get and store geoId for the current material surface
  const GeometryIdentifier geoId = trackState.referenceSurface().geometryId();
  const auto scatteringMapId = scatteringMap.find(geoId);
  if (scatteringMapId == scatteringMap.end()) {
    ACTS_ERROR("No scattering angles found for material surface " << geoId);
    throw std::runtime_error(
        "No scattering angles found for material surface.");
  }

  const double sinThetaLoc = std::sin(trackState.smoothed()[eBoundTheta]);

  // The position, where we need to insert the values in aMatrix and bVector
  const std::size_t deltaPosition = eBoundSize + 2 * nMaterialsHandled;

  const BoundVector& scatteringAngles =
      scatteringMapId->second.scatteringAngles();

  const double invCov = scatteringMapId->second.invCovarianceMaterial();

  // Phi contribution
  extendedSystem.aMatrix()(deltaPosition, deltaPosition) +=
      invCov * sinThetaLoc * sinThetaLoc;
  extendedSystem.bVector()(deltaPosition, 0) -=
      invCov * scatteringAngles[eBoundPhi] * sinThetaLoc;
  extendedSystem.chi2() += invCov * scatteringAngles[eBoundPhi] * sinThetaLoc *
                           scatteringAngles[eBoundPhi] * sinThetaLoc;

  // Theta Contribution
  extendedSystem.aMatrix()(deltaPosition + 1, deltaPosition + 1) += invCov;
  extendedSystem.bVector()(deltaPosition + 1, 0) -=
      invCov * scatteringAngles[eBoundTheta];
  extendedSystem.chi2() +=
      invCov * scatteringAngles[eBoundTheta] * scatteringAngles[eBoundTheta];

  ACTS_VERBOSE(
      "Contributions in addMaterialToGx2fSums:\n"
      << "    invCov:        " << invCov << "\n"
      << "    sinThetaLoc:   " << sinThetaLoc << "\n"
      << "    deltaPosition: " << deltaPosition << "\n"
      << "    Phi:\n"
      << "        scattering angle:     " << scatteringAngles[eBoundPhi] << "\n"
      << "        aMatrix contribution: " << invCov * sinThetaLoc * sinThetaLoc
      << "\n"
      << "        bVector contribution: "
      << invCov * scatteringAngles[eBoundPhi] * sinThetaLoc << "\n"
      << "        chi2sum contribution: "
      << invCov * scatteringAngles[eBoundPhi] * sinThetaLoc *
             scatteringAngles[eBoundPhi] * sinThetaLoc
      << "\n"
      << "    Theta:\n"
      << "        scattering angle:     " << scatteringAngles[eBoundTheta]
      << "\n"
      << "        aMatrix contribution: " << invCov << "\n"
      << "        bVector contribution: "
      << invCov * scatteringAngles[eBoundTheta] << "\n"
      << "        chi2sum contribution: "
      << invCov * scatteringAngles[eBoundTheta] * scatteringAngles[eBoundTheta]
      << "\n");

  return;
}

/// @brief Fill the GX2F system with data from a track
///
/// This function processes a track proxy and updates the aMatrix, bVector, and
/// chi2 values for the GX2F fitting system. It considers material only if
/// multiple scattering is enabled.
///
/// @tparam track_proxy_t The type of the track proxy
///
/// @param track A constant track proxy to inspect
/// @param extendedSystem All parameters of the current equation system
/// @param multipleScattering Flag to consider multiple scattering in the calculation
/// @param scatteringMap Map of geometry identifiers to scattering properties,
///        containing scattering angles and validation status
/// @param geoIdVector A vector to store geometry identifiers for tracking processed elements
/// @param logger A logger instance
template <TrackProxyConcept track_proxy_t>
void fillGx2fSystem(
    const track_proxy_t track, Gx2fSystem& extendedSystem,
    const bool multipleScattering,
    const std::unordered_map<GeometryIdentifier, ScatteringProperties>&
        scatteringMap,
    std::vector<GeometryIdentifier>& geoIdVector, const Logger& logger) {
  std::vector<BoundMatrix> jacobianFromStart;
  jacobianFromStart.emplace_back(BoundMatrix::Identity());

  for (const auto& trackState : track.trackStates()) {
    // Get and store geoId for the current surface
    const GeometryIdentifier geoId = trackState.referenceSurface().geometryId();
    ACTS_DEBUG("Start to investigate trackState on surface " << geoId);
    const auto typeFlags = trackState.typeFlags();
    const bool stateHasMeasurement =
        typeFlags.test(TrackStateFlag::MeasurementFlag);
    const bool stateHasMaterial = typeFlags.test(TrackStateFlag::MaterialFlag);

    // First we figure out, if we would need to look into material
    // surfaces at all. Later, we also check, if the material slab is
    // valid, otherwise we modify this flag to ignore the material
    // completely.
    bool doMaterial = multipleScattering && stateHasMaterial;
    if (doMaterial) {
      const auto scatteringMapId = scatteringMap.find(geoId);
      assert(scatteringMapId != scatteringMap.end() &&
             "No scattering angles found for material surface.");
      doMaterial = doMaterial && scatteringMapId->second.materialIsValid();
    }

    // We only consider states with a measurement (and/or material)
    if (!stateHasMeasurement && !doMaterial) {
      ACTS_DEBUG("    Skip state.");
      continue;
    }

    // update all Jacobians from start
    for (auto& jac : jacobianFromStart) {
      jac = trackState.jacobian() * jac;
    }

    // Handle measurement
    if (stateHasMeasurement) {
      ACTS_DEBUG("    Handle measurement.");

      const auto measDim = trackState.calibratedSize();

      if (measDim < 1 || 6 < measDim) {
        ACTS_ERROR("Can not process state with measurement with "
                   << measDim << " dimensions.");
        throw std::domain_error(
            "Found measurement with less than 1 or more than 6 dimension(s).");
      }

      extendedSystem.ndf() += measDim;

      visit_measurement(measDim, [&](auto N) {
        addMeasurementToGx2fSums<N>(extendedSystem, jacobianFromStart,
                                    trackState, logger);
      });
    }

    // Handle material
    if (doMaterial) {
      ACTS_DEBUG("    Handle material");
      // Add for this material a new Jacobian, starting from this surface.
      jacobianFromStart.emplace_back(BoundMatrix::Identity());

      // Add the material contribution to the system
      addMaterialToGx2fSums(extendedSystem, geoIdVector.size(), scatteringMap,
                            trackState, logger);

      geoIdVector.emplace_back(geoId);
    }
  }
}

/// @brief Count the valid material states in a track for scattering calculations.
///
/// This function counts the valid material surfaces encountered in a track
/// by examining each track state. The count is based on the presence of
/// material flags and the availability of scattering information for each
/// surface.
///
/// @tparam track_proxy_t The type of the track proxy
///
/// @param track A constant track proxy to inspect
/// @param scatteringMap Map of geometry identifiers to scattering properties,
///        containing scattering angles and validation status
/// @param logger A logger instance
template <TrackProxyConcept track_proxy_t>
std::size_t countMaterialStates(
    const track_proxy_t track,
    const std::unordered_map<GeometryIdentifier, ScatteringProperties>&
        scatteringMap,
    const Logger& logger) {
  std::size_t nMaterialSurfaces = 0;
  ACTS_DEBUG("Count the valid material surfaces.");
  for (const auto& trackState : track.trackStates()) {
    const auto typeFlags = trackState.typeFlags();
    const bool stateHasMaterial = typeFlags.test(TrackStateFlag::MaterialFlag);

    if (!stateHasMaterial) {
      continue;
    }

    // Get and store geoId for the current material surface
    const GeometryIdentifier geoId = trackState.referenceSurface().geometryId();

    const auto scatteringMapId = scatteringMap.find(geoId);
    assert(scatteringMapId != scatteringMap.end() &&
           "No scattering angles found for material surface.");
    if (!scatteringMapId->second.materialIsValid()) {
      continue;
    }

    nMaterialSurfaces++;
  }

  return nMaterialSurfaces;
}

/// @brief Solve the gx2f system to get the delta parameters for the update
///
/// This function computes the delta parameters for the GX2F Actor fitting
/// process by solving the linear equation system [a] * delta = b. It uses the
/// column-pivoting Householder QR decomposition for numerical stability.
///
/// @param extendedSystem All parameters of the current equation system
Eigen::VectorXd computeGx2fDeltaParams(const Gx2fSystem& extendedSystem);

/// @brief Update parameters (and scattering angles if applicable)
///
/// @param params Parameters to be updated
/// @param deltaParamsExtended Delta parameters for bound parameter and scattering angles
/// @param nMaterialSurfaces Number of material surfaces in the track
/// @param scatteringMap Map of geometry identifiers to scattering properties,
///        containing all scattering angles and covariances
/// @param geoIdVector Vector of geometry identifiers corresponding to material surfaces
void updateGx2fParams(
    BoundTrackParameters& params, const Eigen::VectorXd& deltaParamsExtended,
    const std::size_t nMaterialSurfaces,
    std::unordered_map<GeometryIdentifier, ScatteringProperties>& scatteringMap,
    const std::vector<GeometryIdentifier>& geoIdVector);

/// @brief Calculate and update the covariance of the fitted parameters
///
/// This function calculates the covariance of the fitted parameters using
/// cov = inv([a])
/// It then updates the first square block of size ndfSystem. This ensures,
/// that we only update the covariance for fitted parameters. (In case of
/// no qop/time fit)
///
/// @param fullCovariancePredicted The covariance matrix to update
/// @param extendedSystem All parameters of the current equation system
void updateGx2fCovarianceParams(BoundMatrix& fullCovariancePredicted,
                                Gx2fSystem& extendedSystem);

/// Global Chi Square fitter (GX2F) implementation.
///
/// @tparam propagator_t Type of the propagation class
///
/// TODO Write description
template <typename propagator_t, typename traj_t>
class Gx2Fitter {
  /// The navigator type
  using Gx2fNavigator = typename propagator_t::Navigator;

  /// The navigator has DirectNavigator type or not
  static constexpr bool isDirectNavigator =
      std::is_same_v<Gx2fNavigator, DirectNavigator>;

 public:
  Gx2Fitter(propagator_t pPropagator,
            std::unique_ptr<const Logger> _logger =
                getDefaultLogger("Gx2Fitter", Logging::INFO))
      : m_propagator(std::move(pPropagator)),
        m_logger{std::move(_logger)},
        m_actorLogger{m_logger->cloneWithSuffix("Actor")},
        m_addToSumLogger{m_logger->cloneWithSuffix("AddToSum")} {}

 private:
  /// The propagator for the transport and material update
  propagator_t m_propagator;

  /// The logger instance
  std::unique_ptr<const Logger> m_logger;
  std::unique_ptr<const Logger> m_actorLogger;
  std::unique_ptr<const Logger> m_addToSumLogger;

  const Logger& logger() const { return *m_logger; }

  /// @brief Propagator Actor plugin for the GX2F
  ///
  /// @tparam parameters_t The type of parameters used for "local" parameters.
  /// @tparam calibrator_t The type of calibrator
  /// @tparam outlier_finder_t Type of the outlier finder class
  ///
  /// The GX2F Actor does not rely on the measurements to be sorted along the
  /// track.
  template <typename parameters_t>
  class Actor {
   public:
    /// Broadcast the result_type
    using result_type = Gx2FitterResult<traj_t>;

    /// The target surface
    const Surface* targetSurface = nullptr;

    /// Allows retrieving measurements for a surface
    const std::map<GeometryIdentifier, SourceLink>* inputMeasurements = nullptr;

    /// Whether to consider multiple scattering.
    bool multipleScattering = false;

    /// Whether to consider energy loss.
    bool energyLoss = false;  /// TODO implement later

    /// Whether to include non-linear correction during global to local
    /// transformation
    FreeToBoundCorrection freeToBoundCorrection;

    /// Input MultiTrajectory
    std::shared_ptr<MultiTrajectory<traj_t>> outputStates;

    /// The logger instance
    const Logger* actorLogger{nullptr};

    /// Logger helper
    const Logger& logger() const { return *actorLogger; }

    Gx2FitterExtensions<traj_t> extensions;

    /// The Surface being
    SurfaceReached targetReached;

    /// Calibration context for the fit
    const CalibrationContext* calibrationContext{nullptr};

    /// The particle hypothesis is needed for estimating scattering angles
    const parameters_t* parametersWithHypothesis = nullptr;

    /// The scatteringMap stores for each visited surface their scattering
    /// properties
    std::unordered_map<GeometryIdentifier, ScatteringProperties>*
        scatteringMap = nullptr;

    /// @brief Gx2f actor operation
    ///
    /// @tparam propagator_state_t is the type of Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param navigator The navigator in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void act(propagator_state_t& state, const stepper_t& stepper,
             const navigator_t& navigator, result_type& result,
             const Logger& /*logger*/) const {
      assert(result.fittedStates && "No MultiTrajectory set");

      // Check if we can stop to propagate
      if (result.measurementStates == inputMeasurements->size()) {
        ACTS_DEBUG("Actor: finish: All measurements have been found.");
        result.finished = true;
      } else if (state.navigation.navigationBreak) {
        ACTS_DEBUG("Actor: finish: navigationBreak.");
        result.finished = true;
      }

      // End the propagation and return to the fitter
      if (result.finished || !result.result.ok()) {
        // Remove the missing surfaces that occur after the last measurement
        if (result.measurementStates > 0) {
          result.missedActiveSurfaces.resize(result.measurementHoles);
        }

        return;
      }

      // We are only interested in surfaces. If we are not on a surface, we
      // continue the navigation
      auto surface = navigator.currentSurface(state.navigation);
      if (surface == nullptr) {
        return;
      }

      ++result.surfaceCount;
      const GeometryIdentifier geoId = surface->geometryId();
      ACTS_DEBUG("Surface " << geoId << " detected.");

      const bool surfaceIsSensitive =
          (surface->associatedDetectorElement() != nullptr);
      const bool surfaceHasMaterial = (surface->surfaceMaterial() != nullptr);
      // First we figure out, if we would need to look into material surfaces at
      // all. Later, we also check, if the material slab is valid, otherwise we
      // modify this flag to ignore the material completely.
      bool doMaterial = multipleScattering && surfaceHasMaterial;

      // Found material - add a scatteringAngles entry if not done yet.
      // Handling will happen later
      if (doMaterial) {
        ACTS_DEBUG("    The surface contains material, ...");

        auto scatteringMapId = scatteringMap->find(geoId);
        if (scatteringMapId == scatteringMap->end()) {
          ACTS_DEBUG("    ... create entry in scattering map.");

          detail::PointwiseMaterialInteraction interaction(surface, state,
                                                           stepper);
          // We need to evaluate the material to create the correct slab
          const bool slabIsValid = interaction.evaluateMaterialSlab(
              state, navigator, MaterialUpdateStage::FullUpdate);
          double invSigma2 = 0.;
          if (slabIsValid) {
            const auto& particle =
                parametersWithHypothesis->particleHypothesis();

            const double sigma =
                static_cast<double>(Acts::computeMultipleScatteringTheta0(
                    interaction.slab, particle.absolutePdg(), particle.mass(),
                    static_cast<float>(
                        parametersWithHypothesis->parameters()[eBoundQOverP]),
                    particle.absoluteCharge()));
            ACTS_VERBOSE(
                "        The Highland formula gives sigma = " << sigma);
            invSigma2 = 1. / std::pow(sigma, 2);
          } else {
            ACTS_VERBOSE("        Material slab is not valid.");
          }

          scatteringMap->emplace(
              geoId, ScatteringProperties{BoundVector::Zero(), invSigma2,
                                          slabIsValid});
          scatteringMapId = scatteringMap->find(geoId);
        } else {
          ACTS_DEBUG("    ... found entry in scattering map.");
        }

        doMaterial = doMaterial && scatteringMapId->second.materialIsValid();
      }

      // Here we handle all measurements
      if (auto sourceLinkIt = inputMeasurements->find(geoId);
          sourceLinkIt != inputMeasurements->end()) {
        ACTS_DEBUG("    The surface contains a measurement.");

        // Transport the covariance to the surface
        stepper.transportCovarianceToBound(state.stepping, *surface,
                                           freeToBoundCorrection);

        // TODO generalize the update of the currentTrackIndex
        auto& fittedStates = *result.fittedStates;

        // Add a <trackStateMask> TrackState entry multi trajectory. This
        // allocates storage for all components, which we will set later.
        typename traj_t::TrackStateProxy trackStateProxy =
            fittedStates.makeTrackState(Gx2fConstants::trackStateMask,
                                        result.lastTrackIndex);
        const std::size_t currentTrackIndex = trackStateProxy.index();

        // Set the trackStateProxy components with the state from the ongoing
        // propagation
        {
          trackStateProxy.setReferenceSurface(surface->getSharedPtr());
          // Bind the transported state to the current surface
          auto res = stepper.boundState(state.stepping, *surface, false,
                                        freeToBoundCorrection);
          if (!res.ok()) {
            result.result = res.error();
            return;
          }
          // Not const since, we might need to update with scattering angles
          auto& [boundParams, jacobian, pathLength] = *res;

          // For material surfaces, we also update the angles with the
          // available scattering information
          if (doMaterial) {
            ACTS_DEBUG("    Update parameters with scattering angles.");
            const auto scatteringMapId = scatteringMap->find(geoId);
            ACTS_VERBOSE(
                "        scatteringAngles: "
                << scatteringMapId->second.scatteringAngles().transpose());
            ACTS_VERBOSE("        boundParams before the update: "
                         << boundParams.parameters().transpose());
            boundParams.parameters() +=
                scatteringMapId->second.scatteringAngles();
            ACTS_VERBOSE("        boundParams after the update: "
                         << boundParams.parameters().transpose());
          }

          // Fill the track state
          trackStateProxy.smoothed() = boundParams.parameters();
          trackStateProxy.smoothedCovariance() = state.stepping.cov;

          trackStateProxy.jacobian() = jacobian;
          trackStateProxy.pathLength() = pathLength;

          if (doMaterial) {
            stepper.update(state.stepping,
                           transformBoundToFreeParameters(
                               trackStateProxy.referenceSurface(),
                               state.geoContext, trackStateProxy.smoothed()),
                           trackStateProxy.smoothed(),
                           trackStateProxy.smoothedCovariance(), *surface);
          }
        }

        // We have smoothed parameters, so calibrate the uncalibrated input
        // measurement
        extensions.calibrator(state.geoContext, *calibrationContext,
                              sourceLinkIt->second, trackStateProxy);

        // Get and set the type flags
        auto typeFlags = trackStateProxy.typeFlags();
        typeFlags.set(TrackStateFlag::ParameterFlag);
        if (surfaceHasMaterial) {
          typeFlags.set(TrackStateFlag::MaterialFlag);
        }

        // Set the measurement type flag
        typeFlags.set(TrackStateFlag::MeasurementFlag);
        // We count the processed measurement
        ++result.processedMeasurements;

        result.lastMeasurementIndex = currentTrackIndex;
        result.lastTrackIndex = currentTrackIndex;

        // TODO check for outlier first
        // We count the state with measurement
        ++result.measurementStates;

        // We count the processed state
        ++result.processedStates;

        // Update the number of holes count only when encountering a
        // measurement
        result.measurementHoles = result.missedActiveSurfaces.size();

        return;
      }

      if (doMaterial) {
        // Here we handle material for multipleScattering. If holes exist, we
        // also handle them already. We create a full trackstate (unlike for
        // simple holes), since we need to evaluate the material later
        ACTS_DEBUG(
            "    The surface contains no measurement, but material and maybe "
            "a hole.");

        // Transport the covariance to the surface
        stepper.transportCovarianceToBound(state.stepping, *surface,
                                           freeToBoundCorrection);

        // TODO generalize the update of the currentTrackIndex
        auto& fittedStates = *result.fittedStates;

        // Add a <trackStateMask> TrackState entry multi trajectory. This
        // allocates storage for all components, which we will set later.
        typename traj_t::TrackStateProxy trackStateProxy =
            fittedStates.makeTrackState(Gx2fConstants::trackStateMask,
                                        result.lastTrackIndex);
        const std::size_t currentTrackIndex = trackStateProxy.index();

        // Set the trackStateProxy components with the state from the ongoing
        // propagation
        {
          trackStateProxy.setReferenceSurface(surface->getSharedPtr());
          // Bind the transported state to the current surface
          auto res = stepper.boundState(state.stepping, *surface, false,
                                        freeToBoundCorrection);
          if (!res.ok()) {
            result.result = res.error();
            return;
          }
          // Not const since, we might need to update with scattering angles
          auto& [boundParams, jacobian, pathLength] = *res;

          // For material surfaces, we also update the angles with the
          // available scattering information
          // We can skip the if here, since we already know, that we do
          // multipleScattering and have material
          ACTS_DEBUG("    Update parameters with scattering angles.");
          const auto scatteringMapId = scatteringMap->find(geoId);
          ACTS_VERBOSE(
              "        scatteringAngles: "
              << scatteringMapId->second.scatteringAngles().transpose());
          ACTS_VERBOSE("        boundParams before the update: "
                       << boundParams.parameters().transpose());
          boundParams.parameters() +=
              scatteringMapId->second.scatteringAngles();
          ACTS_VERBOSE("        boundParams after the update: "
                       << boundParams.parameters().transpose());

          // Fill the track state
          trackStateProxy.smoothed() = boundParams.parameters();
          trackStateProxy.smoothedCovariance() = state.stepping.cov;

          trackStateProxy.jacobian() = jacobian;
          trackStateProxy.pathLength() = pathLength;

          stepper.update(state.stepping,
                         transformBoundToFreeParameters(
                             trackStateProxy.referenceSurface(),
                             state.geoContext, trackStateProxy.smoothed()),
                         trackStateProxy.smoothed(),
                         trackStateProxy.smoothedCovariance(), *surface);
        }

        // Get and set the type flags
        auto typeFlags = trackStateProxy.typeFlags();
        typeFlags.set(TrackStateFlag::ParameterFlag);
        typeFlags.set(TrackStateFlag::MaterialFlag);

        // Set hole only, if we are on a sensitive surface and had
        // measurements before (no holes before the first measurement)
        const bool precedingMeasurementExists = (result.measurementStates > 0);
        if (surfaceIsSensitive && precedingMeasurementExists) {
          ACTS_DEBUG("    Surface is also sensitive. Marked as hole.");
          typeFlags.set(TrackStateFlag::HoleFlag);

          // Count the missed surface
          result.missedActiveSurfaces.push_back(surface);
        }

        result.lastTrackIndex = currentTrackIndex;

        ++result.processedStates;

        return;
      }

      if (surfaceIsSensitive || surfaceHasMaterial) {
        // Here we handle holes. If material hasn't been handled before
        // (because multipleScattering is turned off), we will also handle it
        // here
        if (multipleScattering) {
          ACTS_DEBUG(
              "    The surface contains no measurement, but maybe a hole.");
        } else {
          ACTS_DEBUG(
              "    The surface contains no measurement, but maybe a hole "
              "and/or material.");
        }

        // We only create track states here if there is already a measurement
        // detected (no holes before the first measurement) or if we encounter
        // material
        const bool precedingMeasurementExists = (result.measurementStates > 0);
        if (!precedingMeasurementExists && !surfaceHasMaterial) {
          ACTS_DEBUG(
              "    Ignoring hole, because there are no preceding "
              "measurements.");
          return;
        }

        auto& fittedStates = *result.fittedStates;

        // Add a <trackStateMask> TrackState entry multi trajectory. This
        // allocates storage for all components, which we will set later.
        typename traj_t::TrackStateProxy trackStateProxy =
            fittedStates.makeTrackState(Gx2fConstants::trackStateMask,
                                        result.lastTrackIndex);
        const std::size_t currentTrackIndex = trackStateProxy.index();

        // Set the trackStateProxy components with the state from the
        // ongoing propagation
        {
          trackStateProxy.setReferenceSurface(surface->getSharedPtr());
          // Bind the transported state to the current surface
          auto res = stepper.boundState(state.stepping, *surface, false,
                                        freeToBoundCorrection);
          if (!res.ok()) {
            result.result = res.error();
            return;
          }
          const auto& [boundParams, jacobian, pathLength] = *res;

          // Fill the track state
          trackStateProxy.smoothed() = boundParams.parameters();
          trackStateProxy.smoothedCovariance() = state.stepping.cov;

          trackStateProxy.jacobian() = jacobian;
          trackStateProxy.pathLength() = pathLength;
        }

        // Get and set the type flags
        auto typeFlags = trackStateProxy.typeFlags();
        typeFlags.set(TrackStateFlag::ParameterFlag);
        if (surfaceHasMaterial) {
          ACTS_DEBUG("    It is material.");
          typeFlags.set(TrackStateFlag::MaterialFlag);
        }

        // Set hole only, if we are on a sensitive surface
        if (surfaceIsSensitive && precedingMeasurementExists) {
          ACTS_DEBUG("    It is a hole.");
          typeFlags.set(TrackStateFlag::HoleFlag);
          // Count the missed surface
          result.missedActiveSurfaces.push_back(surface);
        }

        result.lastTrackIndex = currentTrackIndex;

        ++result.processedStates;

        return;
      }

      ACTS_DEBUG("    The surface contains no measurement/material/hole.");
      return;
    }

    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t, typename result_t>
    bool checkAbort(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/, const result_t& result,
                    const Logger& /*logger*/) const {
      if (!result.result.ok() || result.finished) {
        return true;
      }
      return false;
    }
  };

 public:
  /// Fit implementation
  ///
  /// @tparam source_link_iterator_t Iterator type used to pass source links
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam parameters_t Type of parameters used for local parameters
  /// @tparam track_container_t Type of the track container backend
  /// @tparam holder_t Type defining track container backend ownership
  ///
  /// @param it Begin iterator for the fittable uncalibrated measurements
  /// @param end End iterator for the fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param gx2fOptions Gx2FitterOptions steering the fit
  /// @param trackContainer Input track container storage to append into
  /// @note The input measurements are given in the form of @c SourceLink s.
  /// It's the calibrators job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename source_link_iterator_t, typename start_parameters_t,
            typename parameters_t = BoundTrackParameters,
            TrackContainerFrontend track_container_t>
  Result<typename track_container_t::TrackProxy> fit(
      source_link_iterator_t it, source_link_iterator_t end,
      const start_parameters_t& sParameters,
      const Gx2FitterOptions<traj_t>& gx2fOptions,
      track_container_t& trackContainer) const
    requires(!isDirectNavigator)
  {
    // Preprocess Measurements (SourceLinks -> map)
    // To be able to find measurements later, we put them into a map.
    // We need to copy input SourceLinks anyway, so the map can own them.
    ACTS_VERBOSE("Preparing " << std::distance(it, end)
                              << " input measurements");
    std::map<GeometryIdentifier, SourceLink> inputMeasurements;

    for (; it != end; ++it) {
      SourceLink sl = *it;
      auto geoId = gx2fOptions.extensions.surfaceAccessor(sl)->geometryId();
      inputMeasurements.emplace(geoId, std::move(sl));
    }

    // Store, if we want to do multiple scattering. We still need to pass this
    // option to the Actor.
    const bool multipleScattering = gx2fOptions.multipleScattering;

    // Create the ActorList
    using GX2FActor = Actor<parameters_t>;

    using GX2FResult = typename GX2FActor::result_type;
    using Actors = Acts::ActorList<GX2FActor>;

    using PropagatorOptions = typename propagator_t::template Options<Actors>;

    start_parameters_t params = sParameters;
    double chi2sum = 0;
    double oldChi2sum = std::numeric_limits<double>::max();

    // We need to create a temporary track container. We create several times a
    // new track and delete it after updating the parameters. However, if we
    // would work on the externally provided track container, it would be
    // difficult to remove the correct track, if it contains more than one.
    typename track_container_t::TrackContainerBackend trackContainerTempBackend;
    traj_t trajectoryTempBackend;
    TrackContainer trackContainerTemp{trackContainerTempBackend,
                                      trajectoryTempBackend};

    // Create an index of the 'tip' of the track stored in multitrajectory. It
    // is needed outside the update loop. It will be updated with each iteration
    // and used for the final track
    std::size_t tipIndex = Acts::MultiTrajectoryTraits::kInvalid;

    // The scatteringMap stores for each visited surface their scattering
    // properties
    std::unordered_map<GeometryIdentifier, ScatteringProperties> scatteringMap;

    // This will be filled during the updates with the final covariance of the
    // track parameters.
    BoundMatrix fullCovariancePredicted = BoundMatrix::Identity();

    ACTS_VERBOSE("Initial parameters: " << params.parameters().transpose());

    /// Actual Fitting /////////////////////////////////////////////////////////
    ACTS_DEBUG("Start to iterate");

    // Iterate the fit and improve result. Abort after n steps or after
    // convergence.
    // nUpdate is initialized outside to save its state for the track.
    std::size_t nUpdate = 0;
    for (nUpdate = 0; nUpdate < gx2fOptions.nUpdateMax; nUpdate++) {
      ACTS_DEBUG("nUpdate = " << nUpdate + 1 << "/" << gx2fOptions.nUpdateMax);

      // set up propagator and co
      Acts::GeometryContext geoCtx = gx2fOptions.geoContext;
      Acts::MagneticFieldContext magCtx = gx2fOptions.magFieldContext;
      // Set options for propagator
      PropagatorOptions propagatorOptions(geoCtx, magCtx);

      // Add the measurement surface as external surface to the navigator.
      // We will try to hit those surface by ignoring boundary checks.
      for (const auto& [surfaceId, _] : inputMeasurements) {
        propagatorOptions.navigation.insertExternalSurface(surfaceId);
      }

      auto& gx2fActor = propagatorOptions.actorList.template get<GX2FActor>();
      gx2fActor.inputMeasurements = &inputMeasurements;
      gx2fActor.multipleScattering = false;
      gx2fActor.extensions = gx2fOptions.extensions;
      gx2fActor.calibrationContext = &gx2fOptions.calibrationContext.get();
      gx2fActor.actorLogger = m_actorLogger.get();
      gx2fActor.scatteringMap = &scatteringMap;
      gx2fActor.parametersWithHypothesis = &params;

      auto propagatorState = m_propagator.makeState(params, propagatorOptions);

      auto propagatorInitResult = m_propagator.initialize(propagatorState);
      if (!propagatorInitResult.ok()) {
        ACTS_ERROR("Propagation initialization failed: "
                   << propagatorInitResult.error());
        return propagatorInitResult.error();
      }

      auto& r = propagatorState.template get<Gx2FitterResult<traj_t>>();
      r.fittedStates = &trajectoryTempBackend;

      // Clear the track container. It could be more performant to update the
      // existing states, but this needs some more thinking.
      trackContainerTemp.clear();

      auto propagationResult = m_propagator.propagate(propagatorState);

      // Run the fitter
      auto result =
          m_propagator.makeResult(std::move(propagatorState), propagationResult,
                                  propagatorOptions, false);

      if (!result.ok()) {
        ACTS_ERROR("Propagation failed: " << result.error());
        return result.error();
      }

      // TODO Improve Propagator + Actor [allocate before loop], rewrite
      // makeMeasurements
      auto& propRes = *result;
      GX2FResult gx2fResult = std::move(propRes.template get<GX2FResult>());

      if (!gx2fResult.result.ok()) {
        ACTS_INFO("GlobalChiSquareFitter failed in actor: "
                  << gx2fResult.result.error() << ", "
                  << gx2fResult.result.error().message());
        return gx2fResult.result.error();
      }

      auto track = trackContainerTemp.makeTrack();
      tipIndex = gx2fResult.lastMeasurementIndex;

      // It could happen, that no measurements were found. Then the track would
      // be empty and the following operations would be invalid. Usually, this
      // only happens during the first iteration, due to bad initial parameters.
      if (tipIndex == Acts::MultiTrajectoryTraits::kInvalid) {
        ACTS_INFO("Did not find any measurements in nUpdate "
                  << nUpdate + 1 << "/" << gx2fOptions.nUpdateMax);
        return Experimental::GlobalChiSquareFitterError::NotEnoughMeasurements;
      }

      track.tipIndex() = tipIndex;
      track.linkForward();

      // Count the material surfaces, to set up the system. In the multiple
      // scattering case, we need to extend our system.
      const std::size_t nMaterialSurfaces = 0u;

      // We need 6 dimensions for the bound parameters and 2 * nMaterialSurfaces
      // dimensions for the scattering angles.
      const std::size_t dimsExtendedParams = eBoundSize + 2 * nMaterialSurfaces;

      // System that we fill with the information gathered by the actor and
      // evaluate later
      Gx2fSystem extendedSystem{dimsExtendedParams};

      // This vector stores the IDs for each visited material in order. We use
      // it later for updating the scattering angles. We cannot use
      // scatteringMap directly, since we cannot guarantee, that we will visit
      // all stored material in each propagation.
      std::vector<GeometryIdentifier> geoIdVector;

      fillGx2fSystem(track, extendedSystem, false, scatteringMap, geoIdVector,
                     *m_addToSumLogger);

      chi2sum = extendedSystem.chi2();

      // This check takes into account the evaluated dimensions of the
      // measurements. To fit, we need at least NDF+1 measurements. However, we
      // count n-dimensional measurements for n measurements, reducing the
      // effective number of needed measurements. We might encounter the case,
      // where we cannot use some (parts of a) measurements, maybe if we do not
      // support that kind of measurement. This is also taken into account here.
      // We skip the check during the first iteration, since we cannot guarantee
      // to hit all/enough measurement surfaces with the initial parameter
      // guess.
      if ((nUpdate > 0) && !extendedSystem.isWellDefined()) {
        ACTS_INFO("Not enough measurements. Require "
                  << extendedSystem.findRequiredNdf() + 1 << ", but only "
                  << extendedSystem.ndf() << " could be used.");
        return Experimental::GlobalChiSquareFitterError::NotEnoughMeasurements;
      }

      Eigen::VectorXd deltaParamsExtended =
          computeGx2fDeltaParams(extendedSystem);

      ACTS_VERBOSE("aMatrix:\n"
                   << extendedSystem.aMatrix() << "\n"
                   << "bVector:\n"
                   << extendedSystem.bVector() << "\n"
                   << "deltaParamsExtended:\n"
                   << deltaParamsExtended << "\n"
                   << "oldChi2sum = " << oldChi2sum << "\n"
                   << "chi2sum = " << extendedSystem.chi2());

      if ((gx2fOptions.relChi2changeCutOff != 0) && (nUpdate > 0) &&
          (std::abs(extendedSystem.chi2() / oldChi2sum - 1) <
           gx2fOptions.relChi2changeCutOff)) {
        ACTS_DEBUG("Abort with relChi2changeCutOff after "
                   << nUpdate + 1 << "/" << gx2fOptions.nUpdateMax
                   << " iterations.");
        updateGx2fCovarianceParams(fullCovariancePredicted, extendedSystem);
        break;
      }

      if (extendedSystem.chi2() > oldChi2sum + 1e-5) {
        ACTS_DEBUG("chi2 not converging monotonically in update " << nUpdate);
      }

      // If this is the final iteration, update the covariance and break.
      // Otherwise, we would update the scattering angles too much.
      if (nUpdate == gx2fOptions.nUpdateMax - 1) {
        // Since currently most of our tracks converge in 4-5 updates, we want
        // to set nUpdateMax higher than that to guarantee convergence for most
        // tracks. In cases, where we set a smaller nUpdateMax, it's because we
        // want to investigate the behaviour of the fitter before it converges,
        // like in some unit-tests.
        if (gx2fOptions.nUpdateMax > 5) {
          ACTS_INFO("Did not converge in " << gx2fOptions.nUpdateMax
                                           << " updates.");
          return Experimental::GlobalChiSquareFitterError::DidNotConverge;
        }

        updateGx2fCovarianceParams(fullCovariancePredicted, extendedSystem);
        break;
      }

      updateGx2fParams(params, deltaParamsExtended, nMaterialSurfaces,
                       scatteringMap, geoIdVector);
      ACTS_VERBOSE("Updated parameters: " << params.parameters().transpose());

      oldChi2sum = extendedSystem.chi2();
    }
    ACTS_DEBUG("Finished to iterate");
    ACTS_VERBOSE("Final parameters: " << params.parameters().transpose());
    /// Finish Fitting /////////////////////////////////////////////////////////

    /// Actual MATERIAL Fitting ////////////////////////////////////////////////
    ACTS_DEBUG("Start to evaluate material");
    if (multipleScattering) {
      // set up propagator and co
      Acts::GeometryContext geoCtx = gx2fOptions.geoContext;
      Acts::MagneticFieldContext magCtx = gx2fOptions.magFieldContext;
      // Set options for propagator
      PropagatorOptions propagatorOptions(geoCtx, magCtx);

      // Add the measurement surface as external surface to the navigator.
      // We will try to hit those surface by ignoring boundary checks.
      for (const auto& [surfaceId, _] : inputMeasurements) {
        propagatorOptions.navigation.insertExternalSurface(surfaceId);
      }

      auto& gx2fActor = propagatorOptions.actorList.template get<GX2FActor>();
      gx2fActor.inputMeasurements = &inputMeasurements;
      gx2fActor.multipleScattering = true;
      gx2fActor.extensions = gx2fOptions.extensions;
      gx2fActor.calibrationContext = &gx2fOptions.calibrationContext.get();
      gx2fActor.actorLogger = m_actorLogger.get();
      gx2fActor.scatteringMap = &scatteringMap;
      gx2fActor.parametersWithHypothesis = &params;

      auto propagatorState = m_propagator.makeState(params, propagatorOptions);

      auto propagatorInitResult = m_propagator.initialize(propagatorState);
      if (!propagatorInitResult.ok()) {
        ACTS_ERROR("Propagation initialization failed: "
                   << propagatorInitResult.error());
        return propagatorInitResult.error();
      }

      auto& r = propagatorState.template get<Gx2FitterResult<traj_t>>();
      r.fittedStates = &trajectoryTempBackend;

      // Clear the track container. It could be more performant to update the
      // existing states, but this needs some more thinking.
      trackContainerTemp.clear();

      auto propagationResult = m_propagator.propagate(propagatorState);

      // Run the fitter
      auto result =
          m_propagator.makeResult(std::move(propagatorState), propagationResult,
                                  propagatorOptions, false);

      if (!result.ok()) {
        ACTS_ERROR("Propagation failed: " << result.error());
        return result.error();
      }

      // TODO Improve Propagator + Actor [allocate before loop], rewrite
      // makeMeasurements
      auto& propRes = *result;
      GX2FResult gx2fResult = std::move(propRes.template get<GX2FResult>());

      if (!gx2fResult.result.ok()) {
        ACTS_INFO("GlobalChiSquareFitter failed in actor: "
                  << gx2fResult.result.error() << ", "
                  << gx2fResult.result.error().message());
        return gx2fResult.result.error();
      }

      auto track = trackContainerTemp.makeTrack();
      tipIndex = gx2fResult.lastMeasurementIndex;

      // It could happen, that no measurements were found. Then the track would
      // be empty and the following operations would be invalid. Usually, this
      // only happens during the first iteration, due to bad initial parameters.
      if (tipIndex == Acts::MultiTrajectoryTraits::kInvalid) {
        ACTS_INFO("Did not find any measurements in material fit.");
        return Experimental::GlobalChiSquareFitterError::NotEnoughMeasurements;
      }

      track.tipIndex() = tipIndex;
      track.linkForward();

      // Count the material surfaces, to set up the system. In the multiple
      // scattering case, we need to extend our system.
      const std::size_t nMaterialSurfaces =
          countMaterialStates(track, scatteringMap, *m_addToSumLogger);

      // We need 6 dimensions for the bound parameters and 2 * nMaterialSurfaces
      // dimensions for the scattering angles.
      const std::size_t dimsExtendedParams = eBoundSize + 2 * nMaterialSurfaces;

      // System that we fill with the information gathered by the actor and
      // evaluate later
      Gx2fSystem extendedSystem{dimsExtendedParams};

      // This vector stores the IDs for each visited material in order. We use
      // it later for updating the scattering angles. We cannot use
      // scatteringMap directly, since we cannot guarantee, that we will visit
      // all stored material in each propagation.
      std::vector<GeometryIdentifier> geoIdVector;

      fillGx2fSystem(track, extendedSystem, true, scatteringMap, geoIdVector,
                     *m_addToSumLogger);

      chi2sum = extendedSystem.chi2();

      // This check takes into account the evaluated dimensions of the
      // measurements. To fit, we need at least NDF+1 measurements. However, we
      // count n-dimensional measurements for n measurements, reducing the
      // effective number of needed measurements. We might encounter the case,
      // where we cannot use some (parts of a) measurements, maybe if we do not
      // support that kind of measurement. This is also taken into account here.
      // We skip the check during the first iteration, since we cannot guarantee
      // to hit all/enough measurement surfaces with the initial parameter
      // guess.
      if ((nUpdate > 0) && !extendedSystem.isWellDefined()) {
        ACTS_INFO("Not enough measurements. Require "
                  << extendedSystem.findRequiredNdf() + 1 << ", but only "
                  << extendedSystem.ndf() << " could be used.");
        return Experimental::GlobalChiSquareFitterError::NotEnoughMeasurements;
      }

      Eigen::VectorXd deltaParamsExtended =
          computeGx2fDeltaParams(extendedSystem);

      ACTS_VERBOSE("aMatrix:\n"
                   << extendedSystem.aMatrix() << "\n"
                   << "bVector:\n"
                   << extendedSystem.bVector() << "\n"
                   << "deltaParamsExtended:\n"
                   << deltaParamsExtended << "\n"
                   << "oldChi2sum = " << oldChi2sum << "\n"
                   << "chi2sum = " << extendedSystem.chi2());

      chi2sum = extendedSystem.chi2();

      updateGx2fParams(params, deltaParamsExtended, nMaterialSurfaces,
                       scatteringMap, geoIdVector);
      ACTS_VERBOSE("Updated parameters: " << params.parameters().transpose());

      updateGx2fCovarianceParams(fullCovariancePredicted, extendedSystem);
    }
    ACTS_DEBUG("Finished to evaluate material");
    ACTS_VERBOSE(
        "Final parameters after material: " << params.parameters().transpose());
    /// Finish MATERIAL Fitting ////////////////////////////////////////////////

    ACTS_VERBOSE("Final scattering angles:");
    for (const auto& [key, value] : scatteringMap) {
      if (!value.materialIsValid()) {
        continue;
      }
      const auto& angles = value.scatteringAngles();
      ACTS_VERBOSE("    ( " << angles[eBoundTheta] << " | " << angles[eBoundPhi]
                            << " )");
    }

    ACTS_VERBOSE("Final covariance:\n" << fullCovariancePredicted);

    // Propagate again with the final covariance matrix. This is necessary to
    // obtain the propagated covariance for each state.
    // We also need to recheck the result and find the tipIndex, because at this
    // step, we will not ignore the boundary checks for measurement surfaces. We
    // want to create trackstates only on surfaces, that we actually hit.
    if (gx2fOptions.nUpdateMax > 0) {
      ACTS_VERBOSE("Propagate with the final covariance.");
      // update covariance
      params.covariance() = fullCovariancePredicted;

      // set up propagator and co
      Acts::GeometryContext geoCtx = gx2fOptions.geoContext;
      Acts::MagneticFieldContext magCtx = gx2fOptions.magFieldContext;
      // Set options for propagator
      PropagatorOptions propagatorOptions(geoCtx, magCtx);
      auto& gx2fActor = propagatorOptions.actorList.template get<GX2FActor>();
      gx2fActor.inputMeasurements = &inputMeasurements;
      gx2fActor.multipleScattering = multipleScattering;
      gx2fActor.extensions = gx2fOptions.extensions;
      gx2fActor.calibrationContext = &gx2fOptions.calibrationContext.get();
      gx2fActor.actorLogger = m_actorLogger.get();
      gx2fActor.scatteringMap = &scatteringMap;
      gx2fActor.parametersWithHypothesis = &params;

      auto propagatorState = m_propagator.makeState(params, propagatorOptions);

      auto propagatorInitResult = m_propagator.initialize(propagatorState);
      if (!propagatorInitResult.ok()) {
        ACTS_ERROR("Propagation initialization failed: "
                   << propagatorInitResult.error());
        return propagatorInitResult.error();
      }

      auto& r = propagatorState.template get<Gx2FitterResult<traj_t>>();
      r.fittedStates = &trackContainer.trackStateContainer();

      auto propagationResult = m_propagator.propagate(propagatorState);

      // Run the fitter
      auto result =
          m_propagator.makeResult(std::move(propagatorState), propagationResult,
                                  propagatorOptions, false);

      if (!result.ok()) {
        ACTS_ERROR("Propagation failed: " << result.error());
        return result.error();
      }

      auto& propRes = *result;
      GX2FResult gx2fResult = std::move(propRes.template get<GX2FResult>());

      if (!gx2fResult.result.ok()) {
        ACTS_INFO("GlobalChiSquareFitter failed in actor: "
                  << gx2fResult.result.error() << ", "
                  << gx2fResult.result.error().message());
        return gx2fResult.result.error();
      }

      if (tipIndex != gx2fResult.lastMeasurementIndex) {
        ACTS_INFO("Final fit used unreachable measurements.");
        return Experimental::GlobalChiSquareFitterError::
            UsedUnreachableMeasurements;
      }
    }

    if (!trackContainer.hasColumn(
            Acts::hashString(Gx2fConstants::gx2fnUpdateColumn))) {
      trackContainer.template addColumn<std::uint32_t>("Gx2fnUpdateColumn");
    }

    // Prepare track for return
    auto track = trackContainer.makeTrack();
    track.tipIndex() = tipIndex;
    track.parameters() = params.parameters();
    track.covariance() = fullCovariancePredicted;
    track.setReferenceSurface(params.referenceSurface().getSharedPtr());

    if (trackContainer.hasColumn(
            Acts::hashString(Gx2fConstants::gx2fnUpdateColumn))) {
      ACTS_DEBUG("Add nUpdate to track");
      track.template component<std::uint32_t>("Gx2fnUpdateColumn") =
          static_cast<std::uint32_t>(nUpdate);
    }

    // TODO write test for calculateTrackQuantities
    calculateTrackQuantities(track);

    // Set the chi2sum for the track summary manually, since we don't calculate
    // it for each state
    track.chi2() = chi2sum;

    // Return the converted Track
    return track;
  }
};

}  // namespace Acts::Experimental
