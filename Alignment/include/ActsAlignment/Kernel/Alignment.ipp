// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename fitter_t>
template <typename source_link_t, typename start_parameters_t,
          typename fit_options_t>
Acts::Result<ActsAlignment::detail::TrackAlignmentState>
ActsAlignment::Alignment<fitter_t>::evaluateTrackAlignmentState(
    const Acts::GeometryContext& gctx,
    const std::vector<source_link_t>& sourcelinks,
    const start_parameters_t& sParameters, const fit_options_t& fitOptions,
    const std::unordered_map<const Acts::Surface*, size_t>& idxedAlignSurfaces,
    const ActsAlignment::AlignmentMask& alignMask,
    Acts::LoggerWrapper logger) const {
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
      Acts::detail::globalTrackParametersCovariance(
          fitOutput.fittedStates, fitOutput.lastMeasurementIndex);
  // Calculate the alignment state
  const auto alignState = detail::trackAlignmentState(
      gctx, fitOutput.fittedStates, fitOutput.lastMeasurementIndex,
      globalTrackParamsCov, idxedAlignSurfaces, alignMask);
  if (alignState.alignmentDof == 0) {
    ACTS_VERBOSE("No alignment dof on track!");
    return AlignmentError::NoAlignmentDofOnTrack;
  }
  return alignState;
}

template <typename fitter_t>
template <typename trajectory_container_t,
          typename start_parameters_container_t, typename fit_options_t>
void ActsAlignment::Alignment<fitter_t>::calculateAlignmentParameters(
    const trajectory_container_t& trajectoryCollection,
    const start_parameters_container_t& startParametersCollection,
    const fit_options_t& fitOptions,
    ActsAlignment::AlignmentResult& alignResult,
    const ActsAlignment::AlignmentMask& alignMask,
    Acts::LoggerWrapper logger) const {
  // The number of trajectories must be eual to the number of starting
  // parameters
  assert(trajectoryCollection.size() == startParametersCollection.size());

  // The total alignment degree of freedom
  alignResult.alignmentDof =
      alignResult.idxedAlignSurfaces.size() * Acts::eAlignmentSize;
  // Initialize derivative of chi2 w.r.t. aligment parameters for all tracks
  Acts::ActsDynamicVector sumChi2Derivative =
      Acts::ActsDynamicVector::Zero(alignResult.alignmentDof);
  Acts::ActsDynamicMatrix sumChi2SecondDerivative =
      Acts::ActsDynamicMatrix::Zero(alignResult.alignmentDof,
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
    fitOptionsWithRefSurface.referenceSurface = &sParameters.referenceSurface();
    // The result for one single track
    auto evaluateRes = evaluateTrackAlignmentState(
        fitOptions.geoContext, sourcelinks, sParameters,
        fitOptionsWithRefSurface, alignResult.idxedAlignSurfaces, alignMask,
        logger);
    if (not evaluateRes.ok()) {
      ACTS_DEBUG("Evaluation of alignment state for track " << iTraj
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
                dstRow * Acts::eAlignmentSize, dstCol * Acts::eAlignmentSize) +=
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
  Acts::ActsDynamicMatrix sumChi2SecondDerivativeInverse =
      Acts::ActsDynamicMatrix::Zero(alignDof, alignDof);
  sumChi2SecondDerivativeInverse = sumChi2SecondDerivative.inverse();
  if (sumChi2SecondDerivativeInverse.hasNaN()) {
    ACTS_DEBUG("Chi2 second derivative inverse has NaN");
    // return AlignmentError::AlignmentParametersUpdateFailure;
  }

  // Initialize the alignment results
  alignResult.deltaAlignmentParameters =
      Acts::ActsDynamicVector::Zero(alignDof);
  alignResult.alignmentCovariance =
      Acts::ActsDynamicMatrix::Zero(alignDof, alignDof);
  // Solve the linear equation to get alignment parameters change
  alignResult.deltaAlignmentParameters =
      -sumChi2SecondDerivative.fullPivLu().solve(sumChi2Derivative);
  ACTS_VERBOSE("sumChi2SecondDerivative = \n" << sumChi2SecondDerivative);
  ACTS_VERBOSE("sumChi2Derivative = \n" << sumChi2Derivative);
  ACTS_VERBOSE("alignResult.deltaAlignmentParameters \n");

  // Alignment parameters covariance
  alignResult.alignmentCovariance = 2 * sumChi2SecondDerivativeInverse;
  // chi2 change
  alignResult.deltaChi2 = 0.5 * sumChi2Derivative.transpose() *
                          alignResult.deltaAlignmentParameters;
}

template <typename fitter_t>
Acts::Result<void>
ActsAlignment::Alignment<fitter_t>::updateAlignmentParameters(
    const Acts::GeometryContext& gctx,
    const std::vector<Acts::DetectorElementBase*>& alignedDetElements,
    const ActsAlignment::AlignedTransformUpdater& alignedTransformUpdater,
    ActsAlignment::AlignmentResult& alignResult,
    Acts::LoggerWrapper logger) const {
  // Update the aligned transform
  Acts::AlignmentVector deltaAlignmentParam = Acts::AlignmentVector::Zero();
  for (const auto& [surface, index] : alignResult.idxedAlignSurfaces) {
    // 1. The original transform
    const Acts::Vector3& oldCenter = surface->center(gctx);
    const Acts::Transform3& oldTransform = surface->transform(gctx);
    const Acts::RotationMatrix3& oldRotation = oldTransform.rotation();
    // The elements stored below is (rotZ, rotY, rotX)
    const Acts::Vector3& oldEulerAngles = oldRotation.eulerAngles(2, 1, 0);

    // 2. The delta transform
    deltaAlignmentParam = alignResult.deltaAlignmentParameters.segment(
        Acts::eAlignmentSize * index, Acts::eAlignmentSize);
    // The delta translation
    Acts::Vector3 deltaCenter =
        deltaAlignmentParam.segment<3>(Acts::eAlignmentCenter0);
    // The delta Euler angles
    Acts::Vector3 deltaEulerAngles =
        deltaAlignmentParam.segment<3>(Acts::eAlignmentRotation0);

    // 3. The new transform
    const Acts::Vector3 newCenter = oldCenter + deltaCenter;
    // The rotation around global z axis
    Acts::AngleAxis3 rotZ(oldEulerAngles(0) + deltaEulerAngles(2),
                          Acts::Vector3::UnitZ());
    // The rotation around global y axis
    Acts::AngleAxis3 rotY(oldEulerAngles(1) + deltaEulerAngles(1),
                          Acts::Vector3::UnitY());
    // The rotation around global x axis
    Acts::AngleAxis3 rotX(oldEulerAngles(2) + deltaEulerAngles(0),
                          Acts::Vector3::UnitX());
    Eigen::Quaternion<Acts::ActsScalar> newRotation = rotZ * rotY * rotX;
    const Acts::Transform3 newTransform =
        Acts::Translation3(newCenter) * newRotation;

    // 4. Update the aligned transform
    //@Todo: use a better way to handle this (need dynamic cast to inherited
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

template <typename fitter_t>
template <typename trajectory_container_t,
          typename start_parameters_container_t, typename fit_options_t>
Acts::Result<ActsAlignment::AlignmentResult>
ActsAlignment::Alignment<fitter_t>::align(
    const trajectory_container_t& trajectoryCollection,
    const start_parameters_container_t& startParametersCollection,
    const ActsAlignment::AlignmentOptions<fit_options_t>& alignOptions) const {
  const auto& logger = alignOptions.logger;

  // Construct an AlignmentResult object
  AlignmentResult alignResult;

  // Assign index to the alignable surface
  for (unsigned int iDetElement = 0;
       iDetElement < alignOptions.alignedDetElements.size(); iDetElement++) {
    alignResult.idxedAlignSurfaces.emplace(
        &alignOptions.alignedDetElements.at(iDetElement)->surface(),
        iDetElement);
  }
  ACTS_VERBOSE("There are " << alignResult.idxedAlignSurfaces.size()
                            << " detector elements to be aligned");

  // Start the iteration to minimize the chi2
  bool converged = false;
  bool alignmentParametersUpdated = false;
  std::queue<double> recentChi2ONdf;
  ACTS_INFO("Max number of iterations: " << alignOptions.maxIterations);
  for (unsigned int iIter = 0; iIter < alignOptions.maxIterations; iIter++) {
    // Perform the fit to the trajectories and update alignment parameters
    // Initialize the alignment mask (all dof in default)
    AlignmentMask alignMask = AlignmentMask::All;
    // Set the alignment mask
    auto iter_it = alignOptions.iterationState.find(iIter);
    if (iter_it != alignOptions.iterationState.end()) {
      alignMask = iter_it->second;
    }
    // Calculate the alignment parameters delta etc.
    calculateAlignmentParameters(
        trajectoryCollection, startParametersCollection,
        alignOptions.fitOptions, alignResult, alignMask, logger);
    // Screen out the information
    ACTS_INFO("iIter = " << iIter << ", total chi2 = " << alignResult.chi2
                         << ", total measurementDim = "
                         << alignResult.measurementDim
                         << " and average chi2/ndf = "
                         << alignResult.averageChi2ONdf);
    // Check if it has converged against the provided precision
    // 1. either the delta average chi2/ndf in the last few
    // iterations is within tolerance
    if (recentChi2ONdf.size() >=
        alignOptions.deltaAverageChi2ONdfCutOff.first) {
      if (std::abs(recentChi2ONdf.front() - alignResult.averageChi2ONdf) <=
          alignOptions.deltaAverageChi2ONdfCutOff.second) {
        ACTS_INFO(
            "Alignment has converaged with change of chi2/ndf < "
            << alignOptions.deltaAverageChi2ONdfCutOff.second << " in the last "
            << alignOptions.deltaAverageChi2ONdfCutOff.first << " iterations"
            << " after " << iIter << " iteration(s)");
        converged = true;
        break;
      }
      recentChi2ONdf.pop();
    }
    // 2. or the average chi2/ndf (is this correct?)
    if (alignResult.averageChi2ONdf <= alignOptions.averageChi2ONdfCutOff) {
      ACTS_INFO("Alignment has converaged with average chi2/ndf < "
                << alignOptions.averageChi2ONdfCutOff << " after " << iIter
                << " iteration(s)");
      converged = true;
      break;
    }
    // Remove the first element
    // Store the result in the queue
    recentChi2ONdf.push(alignResult.averageChi2ONdf);

    ACTS_INFO("The solved delta of alignmentParameters = \n "
              << alignResult.deltaAlignmentParameters);
    // Not coveraged yet, update the detector element alignment parameters
    auto updateRes = updateAlignmentParameters(
        alignOptions.fitOptions.geoContext, alignOptions.alignedDetElements,
        alignOptions.alignedTransformUpdater, alignResult, logger);
    if (not updateRes.ok()) {
      ACTS_ERROR("Update alignment parameters failed: " << updateRes.error());
      return updateRes.error();
    }
    alignmentParametersUpdated = true;
  }  // end of all iterations

  // Alignment failure if not converged
  if (not converged) {
    ACTS_ERROR("Alignment is not converged.");
    alignResult.result = AlignmentError::ConvergeFailure;
  }

  // Screen out the final aligned parameters
  // @todo
  if (alignmentParametersUpdated) {
    unsigned int iDetElement = 0;
    for (const auto& det : alignOptions.alignedDetElements) {
      const auto& surface = &det->surface();
      const auto& transform =
          det->transform(alignOptions.fitOptions.geoContext);
      // write it to the result
      alignResult.alignedParameters.emplace(det, transform);
      const auto& translation = transform.translation();
      const auto& rotation = transform.rotation();
      const Acts::Vector3 rotAngles = rotation.eulerAngles(2, 1, 0);
      ACTS_VERBOSE("Detector element with surface "
                   << surface->geometryId()
                   << " has aligned geometry position as below:");
      ACTS_VERBOSE("Center (cenX, cenY, cenZ) = " << translation.transpose());
      ACTS_VERBOSE(
          "Euler angles (rotZ, rotY, rotX) = " << rotAngles.transpose());
      ACTS_VERBOSE("Rotation marix = \n" << rotation);
      iDetElement++;
    }
  } else {
    ACTS_DEBUG("Alignment parameters is not updated.");
  }

  return alignResult;
}
