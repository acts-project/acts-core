// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>
#include <map>
#include <queue>
#include <vector>

#include "ActsAlignment/Kernel/AlignmentError.hpp"
#include "ActsAlignment/Kernel/detail/AlignmentEngine.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/KalmanGlobalCovariance.hpp"
#include "Acts/Utilities/AlignmentDefinitions.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

namespace ActsAlignment {
using AlignedTransformUpdater =
    std::function<bool(Acts::DetectorElementBase*, const Acts::GeometryContext&,
                       const Acts::Transform3D&)>;
///
/// @brief Options for align() call
///
/// @tparam fit_options_t The fit options type
///
template <typename fit_options_t>
struct AlignmentOptions {
  /// Deleted default constructor
  AlignmentOptions() = delete;

  /// AlignmentOptions
  ///
  /// @param fOptions The fit options
  /// @param aTransformUpdater The updater to update aligned transform
  /// @param aDetElements The alignable detector elements
  /// @param chi2CufOff The alignment chi2 tolerance
  /// @param maxIters The alignment maximum iterations
  //  @param iterState The alignment mask for each
  //  iteration @Todo: use a json file to handle this
  AlignmentOptions(
      const fit_options_t& fOptions,
      const AlignedTransformUpdater& aTransformUpdater,
      const std::vector<Acts::DetectorElementBase*>& aDetElements = {},
      double chi2CutOff = 0.05,
      const std::pair<size_t, double>& deltaChi2CutOff = {10, 0.00001},
      size_t maxIters = 5,
      const std::map<unsigned int, std::bitset<Acts::eAlignmentSize>>&
          iterState = {})
      : fitOptions(fOptions),
        alignedTransformUpdater(aTransformUpdater),
        alignedDetElements(aDetElements),
        averageChi2ONdfCutOff(chi2CutOff),
        deltaAverageChi2ONdfCutOff(deltaChi2CutOff),
        maxIterations(maxIters),
        iterationState(iterState) {}

  // The fit options
  fit_options_t fitOptions;

  /// The updater to the aligned transform
  AlignedTransformUpdater alignedTransformUpdater = nullptr;

  // The detector elements to be aligned
  std::vector<Acts::DetectorElementBase*> alignedDetElements;

  // The alignment tolerance
  double averageChi2ONdfCutOff = 0.05;

  // The delta of average chi2/ndf within a couple of iterations to determine if
  // alignment is converged
  std::pair<size_t, double> deltaAverageChi2ONdfCutOff = {10, 0.00001};

  // The maximum number of iterations to run alignment
  size_t maxIterations = 5;

  // The covariance of the source links and alignment mask for different
  // iterations
  std::map<unsigned int, std::bitset<Acts::eAlignmentSize>> iterationState;
};

/// @brief Alignment result struct
///
struct AlignmentResult {
  // The change of alignment parameters
  Acts::ActsVectorX<Acts::BoundScalar> deltaAlignmentParameters;

  // The aligned parameters
  std::unordered_map<Acts::DetectorElementBase*, Acts::Transform3D>
      alignedParameters;

  // The covariance of alignment parameters
  Acts::ActsMatrixX<Acts::BoundScalar> alignmentCovariance;

  // The avarage chi2/ndf (ndf is the measurement dim)
  double averageChi2ONdf = std::numeric_limits<double>::max();

  // The delta chi2
  double deltaChi2 = std::numeric_limits<double>::max();

  // The chi2
  double chi2 = 0;

  // The measurement dim from all track
  size_t measurementDim = 0;

  // The number of alignment dof
  size_t alignmentDof = 0;

  // The number of tracks used for alignment
  size_t numTracks = 0;

  // The indexed alignable surfaces
  std::unordered_map<const Acts::Surface*, size_t> idxedAlignSurfaces;

  Acts::Result<void> result{Acts::Result<void>::success()};
};

/// @brief KalmanFitter-based alignment implementation
///
/// @tparam fitter_t Type of the fitter class
template <typename fitter_t>
struct Alignment {
  /// Default constructor is deleted
  Alignment() = delete;

  /// Constructor from arguments
  Alignment(fitter_t fitter,
            std::unique_ptr<const Acts::Logger> logger =
                Acts::getDefaultLogger("Alignment", Acts::Logging::INFO))
      : m_fitter(std::move(fitter)), m_logger(logger.release()) {}

  /// @brief evaluate alignment state for a single track
  ///
  /// @tparam source_link_t Source link type identifying uncalibrated input
  /// measurements.
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam fit_options_t The fit options type
  ///
  /// @param gctx The current geometry context object
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param fitOptions The fit Options steering the fit
  /// @param idxedAlignSurfaces The idxed surfaces to be aligned
  /// @param alignMask The alignment mask (same for all detector element for the
  /// moment)
  ///
  /// @result The alignment state for a single track
  template <typename source_link_t, typename start_parameters_t,
            typename fit_options_t>
  Acts::Result<detail::TrackAlignmentState> evaluateTrackAlignmentState(
      const Acts::GeometryContext& gctx,
      const std::vector<source_link_t>& sourcelinks,
      const start_parameters_t& sParameters, const fit_options_t& fitOptions,
      const std::unordered_map<const Acts::Surface*, size_t>&
          idxedAlignSurfaces,
      const std::bitset<Acts::eAlignmentSize>& alignMask) const {
    // Perform the fit
    auto fitRes = m_fitter.fit(sourcelinks, sParameters, fitOptions);
    if (not fitRes.ok()) {
      ACTS_WARNING("Fit failure");
      return fitRes.error();
    }
    // The fit results
    const auto& fitOutput = fitRes.value();
    // Calculate the global track parameters covariance with the fitted track
    const auto& globalTrackParamsCov =
        Acts::detail::globalTrackParametersCovariance(fitOutput.fittedStates,
                                                      fitOutput.trackTip);
    // Calculate the alignment state
    const auto alignState = detail::trackAlignmentState(
        gctx, fitOutput.fittedStates, fitOutput.trackTip, globalTrackParamsCov,
        idxedAlignSurfaces, alignMask);
    if (alignState.alignmentDof == 0) {
      ACTS_VERBOSE("No alignment dof on track");
      return AlignmentError::NoAlignmentDofOnTrack;
    }
    return alignState;
  }

  /// @brief calculate the alignment parameters delta
  ///
  /// @tparam trajectory_container_t The trajectories container type
  /// @tparam start_parameters_t The initial parameters container type
  /// @tparam fit_options_t The fit options type
  ///
  /// @param trajectoryCollection The collection of trajectories as input of
  /// fitting
  /// @param startParametersCollection The collection of starting parameters as
  /// input of fitting
  /// @param fitOptions The fit Options steering the fit
  /// @param alignResult [in, out] The aligned result
  /// @param alignMask The alignment mask (same for all measurements now)
  template <typename trajectory_container_t,
            typename start_parameters_container_t, typename fit_options_t>
  void calculateAlignmentParameters(
      trajectory_container_t& trajectoryCollection,
      const start_parameters_container_t& startParametersCollection,
      const fit_options_t& fitOptions, AlignmentResult& alignResult,
      const std::bitset<Acts::eAlignmentSize>& alignMask =
          std::bitset<Acts::eAlignmentSize>(std::string("111111"))) const {
    // The number of trajectories must be eual to the number of starting
    // parameters
    assert(trajectoryCollection.size() == startParametersCollection.size());

    // The total alignment degree of freedom
    alignResult.alignmentDof =
        alignResult.idxedAlignSurfaces.size() * Acts::eAlignmentSize;
    // Initialize derivative of chi2 w.r.t. aligment parameters for all tracks
    Acts::ActsVectorX<Acts::BoundScalar> sumChi2Derivative =
        Acts::ActsVectorX<Acts::BoundScalar>::Zero(alignResult.alignmentDof);
    Acts::ActsMatrixX<Acts::BoundScalar> sumChi2SecondDerivative =
        Acts::ActsMatrixX<Acts::BoundScalar>::Zero(alignResult.alignmentDof,
                                                   alignResult.alignmentDof);
    // Copy the fit options
    fit_options_t fitOptionsWithRefSurface = fitOptions;
    // Calculate contribution to chi2 derivatives from all input trajectories
    // @Todo: How to update the source link error iteratively?
    alignResult.chi2 = 0;
    alignResult.measurementDim = 0;
    alignResult.numTracks = trajectoryCollection.size();
    double sumChi2ONdf = 0;
    for (unsigned int iTraj = 0; iTraj < trajectoryCollection.size(); iTraj++) {
      const auto& sourcelinks = trajectoryCollection.at(iTraj);
      const auto& sParameters = startParametersCollection.at(iTraj);
      // Set the target surface
      fitOptionsWithRefSurface.referenceSurface =
          &sParameters.referenceSurface();
      // The result for one single track
      auto evaluateRes = evaluateTrackAlignmentState(
          fitOptions.geoContext, sourcelinks, sParameters,
          fitOptionsWithRefSurface, alignResult.idxedAlignSurfaces, alignMask);
      if (not evaluateRes.ok()) {
        ACTS_WARNING("Evaluation of alignment state for track " << iTraj
                                                                << " failed");
        continue;
      }
      const auto& alignState = evaluateRes.value();
      for (const auto& [rowSurface, rows] : alignState.alignedSurfaces) {
        const auto& [dstRow, srcRow] = rows;
        // Fill the results into full chi2 derivative matrixs
        sumChi2Derivative.segment<Acts::eAlignmentSize>(dstRow *
                                                        Acts::eAlignmentSize) +=
            alignState.alignmentToChi2Derivative.segment(
                srcRow * Acts::eAlignmentSize, Acts::eAlignmentSize);

        for (const auto& [colSurface, cols] : alignState.alignedSurfaces) {
          const auto& [dstCol, srcCol] = cols;
          sumChi2SecondDerivative
              .block<Acts::eAlignmentSize, Acts::eAlignmentSize>(
                  dstRow * Acts::eAlignmentSize,
                  dstCol * Acts::eAlignmentSize) +=
              alignState.alignmentToChi2SecondDerivative.block(
                  srcRow * Acts::eAlignmentSize, srcCol * Acts::eAlignmentSize,
                  Acts::eAlignmentSize, Acts::eAlignmentSize);
        }
      }
      alignResult.chi2 += alignState.chi2;
      alignResult.measurementDim += alignState.measurementDim;
      sumChi2ONdf += alignState.chi2 / alignState.measurementDim;
    }
    alignResult.averageChi2ONdf = sumChi2ONdf / alignResult.numTracks;

    // Get the inverse of chi2 second derivative matrix (we need this to
    // calculate the covariance of the alignment parameters)
    // @Todo: use more stable method for solving the inverse
    size_t alignDof = alignResult.alignmentDof;
    Acts::ActsMatrixX<Acts::BoundScalar> sumChi2SecondDerivativeInverse =
        Acts::ActsMatrixX<Acts::BoundScalar>::Zero(alignDof, alignDof);
    sumChi2SecondDerivativeInverse = sumChi2SecondDerivative.inverse();
    if (sumChi2SecondDerivativeInverse.hasNaN()) {
      ACTS_WARNING("Chi2 second derivative inverse has NaN");
      // return AlignmentError::AlignmentParametersUpdateFailure;
    }

    // Initialize the alignment results
    alignResult.deltaAlignmentParameters =
        Acts::ActsVectorX<Acts::BoundScalar>::Zero(alignDof);
    alignResult.alignmentCovariance =
        Acts::ActsMatrixX<Acts::BoundScalar>::Zero(alignDof, alignDof);
    // Solve the linear equation to get alignment parameters change
    alignResult.deltaAlignmentParameters =
        -sumChi2SecondDerivative.fullPivLu().solve(sumChi2Derivative);
    ACTS_INFO("The solved delta of alignmentParameters = \n "
              << alignResult.deltaAlignmentParameters);

    // Alignment parameters covariance
    alignResult.alignmentCovariance = 2 * sumChi2SecondDerivativeInverse;
    // chi2 change
    alignResult.deltaChi2 = 0.5 * sumChi2Derivative.transpose() *
                            alignResult.deltaAlignmentParameters;
  }

  /// @brief update the detector element alignment parameters
  ///
  /// @param gctx The geometry context
  /// @param alignedDetElements The detector elements to be aligned
  /// @param alignedTransformUpdater The updater for updating the aligned
  /// @param alignResult [in, out] The aligned result
  Acts::Result<void> updateAlignmentParameters(
      const Acts::GeometryContext& gctx,
      const std::vector<Acts::DetectorElementBase*>& alignedDetElements,
      const AlignedTransformUpdater& alignedTransformUpdater,
      AlignmentResult& alignResult) const {
    // Update the aligned transform
    Acts::AlignmentVector deltaAlignmentParam = Acts::AlignmentVector::Zero();
    for (const auto& [surface, index] : alignResult.idxedAlignSurfaces) {
      // (1) The original transform
      const Acts::Vector3D& oldCenter = surface->center(gctx);
      const Acts::Transform3D& oldTransform = surface->transform(gctx);
      const Acts::RotationMatrix3D& oldRotation = oldTransform.rotation();
      // The elements stored below is (rotZ, rotY, rotX)
      const Acts::Vector3D& oldEulerAngles = oldRotation.eulerAngles(2, 1, 0);

      // (2) The delta transform
      deltaAlignmentParam = alignResult.deltaAlignmentParameters.segment(
          Acts::eAlignmentSize * index, Acts::eAlignmentSize);
      // The delta translation
      Acts::Vector3D deltaCenter =
          deltaAlignmentParam.segment<3>(Acts::eAlignmentCenter0);
      // The delta Euler angles
      Acts::Vector3D deltaEulerAngles =
          deltaAlignmentParam.segment<3>(Acts::eAlignmentRotation0);

      // (3) The new transform
      const Acts::Vector3D newCenter = oldCenter + deltaCenter;
      // The rotation around global z axis
      Acts::AngleAxis3D rotZ(oldEulerAngles(0) + deltaEulerAngles(2),
                             Acts::Vector3D::UnitZ());
      // The rotation around global y axis
      Acts::AngleAxis3D rotY(oldEulerAngles(1) + deltaEulerAngles(1),
                             Acts::Vector3D::UnitY());
      // The rotation around global x axis
      Acts::AngleAxis3D rotX(oldEulerAngles(2) + deltaEulerAngles(0),
                             Acts::Vector3D::UnitX());
      Acts::Rotation3D newRotation = rotZ * rotY * rotX;
      const Acts::Transform3D newTransform =
          Acts::Translation3D(newCenter) * newRotation;

      // Update the aligned transform
      //@Todo: use a better way to handle this(need dynamic cast to inherited
      // detector element type)
      ACTS_VERBOSE("Delta of alignment parameters at element "
                   << index << "= \n"
                   << deltaAlignmentParam);
      bool updated = alignedTransformUpdater(alignedDetElements.at(index), gctx,
                                             newTransform);
      if (not updated) {
        ACTS_ERROR("Update alignment parameters for detector element failed");
        return AlignmentError::AlignmentParametersUpdateFailure;
      }
    }

    return Acts::Result<void>::success();
  }

  /// @brief Alignment implementation
  ///
  /// @tparam trajectory_container_t The trajectories container type
  /// @tparam start_parameters_t The initial parameters container type
  /// @tparam fit_options_t The fit options type
  ///
  /// @param trajectoryCollection The collection of trajectories as input of
  /// fitting
  /// @param startParametersCollection The collection of starting parameters as
  /// input of fitting
  /// @param alignOptions The alignment options
  ///
  /// @result The alignment result
  template <typename trajectory_container_t,
            typename start_parameters_container_t, typename fit_options_t>
  Acts::Result<AlignmentResult> align(
      trajectory_container_t& trajectoryCollection,
      const start_parameters_container_t& startParametersCollection,
      const AlignmentOptions<fit_options_t>& alignOptions) const {
    // Construct an AlignmentResult object
    AlignmentResult alignResult;

    // Assign index to the alignable surface
    for (unsigned int iDetElement = 0;
         iDetElement < alignOptions.alignedDetElements.size(); iDetElement++) {
      alignResult.idxedAlignSurfaces.emplace(
          &alignOptions.alignedDetElements.at(iDetElement)->surface(),
          iDetElement);
    }

    // Start the iteration to minimize the chi2
    bool converged = false;
    std::queue<double> recentChi2ONdf;
    ACTS_INFO("Max number of iterations: " << alignOptions.maxIterations);
    for (unsigned int iIter = 0; iIter < alignOptions.maxIterations; iIter++) {
      // Perform the fit to the trajectories and update alignment parameters
      std::bitset<Acts::eAlignmentSize> alignmentMask(std::string("111111"));
      auto iter_it = alignOptions.iterationState.find(iIter);
      if (iter_it != alignOptions.iterationState.end()) {
        alignmentMask = iter_it->second;
      }
      // Calculate the alignment parameters delta etc.
      calculateAlignmentParameters(
          trajectoryCollection, startParametersCollection,
          alignOptions.fitOptions, alignResult, alignmentMask);
      // Screen out the information
      ACTS_INFO("iIter = " << iIter << ", total chi2 = " << alignResult.chi2
                           << ", total measurementDim = "
                           << alignResult.measurementDim
                           << "\n Average chi2/ndf = "
                           << alignResult.averageChi2ONdf);
      // Check if it has converged against the provided precision
      // (1) firstly check the average chi2/ndf (is this correct?)
      if (alignResult.averageChi2ONdf <= alignOptions.averageChi2ONdfCutOff) {
        ACTS_INFO("Alignment has converaged with average chi2/ndf smaller than "
                  << alignOptions.averageChi2ONdfCutOff);
        converged = true;
        break;
      }
      // (2) secondly check if the delta average chi2/ndf in the last few
      // iterations is within tolerance
      if (recentChi2ONdf.size() >=
          alignOptions.deltaAverageChi2ONdfCutOff.first) {
        if (std::abs(recentChi2ONdf.front() - alignResult.averageChi2ONdf) <=
            alignOptions.deltaAverageChi2ONdfCutOff.second) {
          ACTS_INFO(
              "Alignment has converaged with change of chi2/ndf smaller than "
              << alignOptions.deltaAverageChi2ONdfCutOff.second
              << " in the latest "
              << alignOptions.deltaAverageChi2ONdfCutOff.first
              << " iterations");
          converged = true;
          break;
        }
        // Remove the first element
        recentChi2ONdf.pop();
      }
      // Store the result in the queue
      recentChi2ONdf.push(alignResult.averageChi2ONdf);

      // Not coveraged yet, update the detector element alignment parameters
      auto updateRes = updateAlignmentParameters(
          alignOptions.fitOptions.geoContext, alignOptions.alignedDetElements,
          alignOptions.alignedTransformUpdater, alignResult);
      if (not updateRes.ok()) {
        ACTS_ERROR("Update alignment parameters failed: " << updateRes.error());
        return updateRes.error();
      }
    }  // end of all iterations

    // Alignment failure if not converged
    if (not converged) {
      ACTS_ERROR("Alignment is not converged.");
      alignResult.result = AlignmentError::ConvergeFailure;
    }

    // Screen out the final aligned parameters
    unsigned int iDetElement = 0;
    for (const auto& det : alignOptions.alignedDetElements) {
      const auto& surface = &det->surface();
      const auto& transform =
          det->transform(alignOptions.fitOptions.geoContext);
      // write it to the result
      alignResult.alignedParameters.emplace(det, transform);
      const auto& translation = transform.translation();
      const auto& rotation = transform.rotation();
      const Acts::Vector3D rotAngles = rotation.eulerAngles(2, 1, 0);
      ACTS_INFO("Detector element with surface "
                << surface->geometryId()
                << " has aligned geometry position as below:");
      ACTS_INFO("Center (cenX, cenY, cenZ) = " << translation.transpose());
      ACTS_INFO("Euler angles (rotZ, rotY, rotX) = " << rotAngles.transpose());
      ACTS_INFO("Rotation marix = \n" << rotation);
      iDetElement++;
    }

    return alignResult;
  }

 private:
  // The fitter
  fitter_t m_fitter;

  /// Logger getter to support macros
  const Acts::Logger& logger() const { return *m_logger; }

  /// Owned logging instance
  std::shared_ptr<const Acts::Logger> m_logger;
};
}  // namespace ActsAlignment
