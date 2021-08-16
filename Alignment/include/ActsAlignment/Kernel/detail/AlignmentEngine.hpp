// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"

#include <unordered_map>

namespace ActsAlignment {
namespace detail {

using namespace Acts;
///
///@brief struct to store info needed for track-based alignment
///
struct TrackAlignmentState {
  // The dimension of measurements
  size_t measurementDim = 0;

  // The dimension of track parameters
  size_t trackParametersDim = 0;

  // The contributed alignment degree of freedom
  size_t alignmentDof = 0;

  // The measurements covariance
  ActsDynamicMatrix measurementCovariance;

  // The track parameters covariance
  ActsDynamicMatrix trackParametersCovariance;

  // The projection matrix
  ActsDynamicMatrix projectionMatrix;

  // The residual
  ActsDynamicVector residual;

  // The covariance of residual
  ActsDynamicMatrix residualCovariance;

  // The chi2
  double chi2 = 0;

  // The derivative of residual w.r.t. alignment parameters
  ActsDynamicMatrix alignmentToResidualDerivative;

  // The derivative of chi2 w.r.t. alignment parameters
  ActsDynamicVector alignmentToChi2Derivative;

  // The second derivative of chi2 w.r.t. alignment parameters
  ActsDynamicMatrix alignmentToChi2SecondDerivative;

  // The alignable surfaces on the track and their indices in both the global
  // alignable surfaces pool and those relevant with this track
  std::unordered_map<const Surface*, std::pair<size_t, size_t>> alignedSurfaces;
};

/// Reset some columns of the alignment to bound derivative to zero if the
/// relevant degree of freedom is fixed
///
/// @param alignToBound The alignment to bound parameters derivative
/// @param mask The alignment mask
void resetAlignmentDerivative(Acts::AlignmentToBoundMatrix& alignToBound,
                              AlignmentMask mask);

///
/// Calculate the first and second derivative of chi2 w.r.t. alignment
/// parameters for a single track
///
/// Suppose there are n measurements on the track, and m (m<=n) of them are on
/// alignable surface, then (eAlignmentSize*m) alignment parameters
/// will be involved for this particular track, i.e. this track will contribute
/// to at most (eAlignmentSize*m*2) elements of the full chi2
/// second derivative matrix
///
/// @tparam source_link_t The source link type of the trajectory
/// @tparam parameters_t The track parameters type
///
/// @param gctx The current geometry context object
/// @param multiTraj The MultiTrajectory containing the trajectory to be
/// investigated
/// @param entryIndex The trajectory entry index
/// @param globalTrackParamsCov The global track parameters covariance for a
/// single track and the starting row/column for smoothed states. This contains
/// all smoothed track states including those non-measurement states. Selection
/// of certain rows/columns for measurement states is needed.
/// @param idxedAlignSurfaces The indexed surfaces to be aligned
///
/// @return The track alignment state containing fundamental alignment
/// ingredients
template <typename source_link_t, typename parameters_t = BoundTrackParameters>
TrackAlignmentState trackAlignmentState(
    const GeometryContext& gctx,
    const Acts::MultiTrajectory<source_link_t>& multiTraj,
    const size_t& entryIndex,
    const std::pair<ActsDynamicMatrix, std::unordered_map<size_t, size_t>>&
        globalTrackParamsCov,
    const std::unordered_map<const Surface*, size_t>& idxedAlignSurfaces,
    const AlignmentMask& alignMask) {
  using CovMatrix = typename parameters_t::CovarianceMatrix;

  // Construct an alignment state
  TrackAlignmentState alignState;

  // Remember the index within the trajectory and whether it's alignable
  std::vector<std::pair<size_t, bool>> measurementStates;
  measurementStates.reserve(15);
  // Number of smoothed states on the track
  size_t nSmoothedStates = 0;
  // Number of alignable surfaces on the track
  size_t nAlignSurfaces = 0;

  // Visit the track states on the track
  multiTraj.visitBackwards(entryIndex, [&](const auto& ts) {
    // Remember the number of smoothed states
    if (ts.hasSmoothed()) {
      nSmoothedStates++;
    } else {
      // @note: this should in principle never happen now. But still keep it as a note
      return true;
    }

    // Only measurement states matter (we can't align non-measurement states,
    // no?)
    if (not ts.typeFlags().test(TrackStateFlag::MeasurementFlag)) {
      return true;
    }
    // Check if the reference surface is to be aligned
    bool isAlignable = false;
    const auto surface = &ts.referenceSurface();
    auto it = idxedAlignSurfaces.find(surface);
    if (it != idxedAlignSurfaces.end()) {
      isAlignable = true;
      // Remember the surface and its index
      alignState.alignedSurfaces[surface].first = it->second;
      nAlignSurfaces++;
    }
    // Rember the index of the state within the trajectory and whether it's
    // alignable
    measurementStates.push_back({ts.index(), isAlignable});
    // Add up measurement dimension
    alignState.measurementDim += ts.calibratedSize();
    return true;
  });

  // Return now if the track contains no alignable surfaces
  if (nAlignSurfaces == 0) {
    return alignState;
  }

  // The alignment degree of freedom
  alignState.alignmentDof = eAlignmentSize * nAlignSurfaces;
  // Dimension of global track parameters (from only measurement states)
  alignState.trackParametersDim = eBoundSize * measurementStates.size();

  // Initialize the alignment matrixs with components from the measurement
  // states
  // The measurement covariance
  alignState.measurementCovariance = ActsDynamicMatrix::Zero(
      alignState.measurementDim, alignState.measurementDim);
  // The bound parameters to measurement projection matrix
  alignState.projectionMatrix = ActsDynamicMatrix::Zero(
      alignState.measurementDim, alignState.trackParametersDim);
  // The derivative of residual w.r.t. alignment parameters
  alignState.alignmentToResidualDerivative = ActsDynamicMatrix::Zero(
      alignState.measurementDim, alignState.alignmentDof);
  // The track parameters covariance
  alignState.trackParametersCovariance = ActsDynamicMatrix::Zero(
      alignState.trackParametersDim, alignState.trackParametersDim);
  // The residual
  alignState.residual = ActsDynamicVector::Zero(alignState.measurementDim);

  // Unpack global track parameters covariance and the starting row/column for
  // all smoothed states.
  // Note that the dimension of provided global track parameters covariance
  // should be same as eBoundSize * nSmoothedStates
  const auto& [sourceTrackParamsCov, stateRowIndices] = globalTrackParamsCov;

  // Loop over the measurement states to fill those alignment matrixs
  // This is done in reverse order
  size_t iMeasurement = alignState.measurementDim;
  size_t iParams = alignState.trackParametersDim;
  size_t iSurface = nAlignSurfaces;
  for (const auto& [rowStateIndex, isAlignable] : measurementStates) {
    const auto& state = multiTraj.getTrackState(rowStateIndex);
    const size_t measdim = state.calibratedSize();
    // Update index of current measurement and parameter
    iMeasurement -= measdim;
    iParams -= eBoundSize;
    // (a) Get and fill the measurement covariance matrix
    const ActsDynamicMatrix measCovariance =
        state.effectiveCalibratedCovariance();
    alignState.measurementCovariance.block(iMeasurement, iMeasurement, measdim,
                                           measdim) = measCovariance;

    // (b) Get and fill the bound parameters to measurement projection matrix
    const ActsDynamicMatrix H = state.effectiveProjector();
    alignState.projectionMatrix.block(iMeasurement, iParams, measdim,
                                      eBoundSize) = H;
    // (c) Get and fill the residual
    alignState.residual.segment(iMeasurement, measdim) =
        state.effectiveCalibrated() - H * state.smoothed();

    // (d) Get the derivative of alignment parameters w.r.t. measurement
    // or residual
    if (isAlignable) {
      iSurface -= 1;
      const auto surface = &state.referenceSurface();
      alignState.alignedSurfaces.at(surface).second = iSurface;
      // The free parameters transformed from the smoothed parameters
      const FreeVector freeParams =
          Acts::MultiTrajectoryHelpers::freeSmoothed(gctx, state);
      // The direction
      const Vector3 direction = freeParams.segment<3>(eFreeDir0);
      // The derivative of free parameters w.r.t. path length. @note Here, we
      // assumes a linear track model, i.e. negecting the change of track
      // direction. Otherwise, we need to know the magnetic field at the free
      // parameters
      FreeVector pathDerivative = FreeVector::Zero();
      pathDerivative.head<3>() = direction;
      // Get the derivative of bound parameters w.r.t. alignment parameters
      AlignmentToBoundMatrix alignToBound =
          surface->alignmentToBoundDerivative(gctx, freeParams, pathDerivative);
      // Set the degree of freedom per surface.
      // @Todo: don't allocate memory for fixed degree of freedom and consider surface/layer/volume wise align mask (instead of using global mask as now)
      resetAlignmentDerivative(alignToBound, alignMask);

      // Residual is calculated as the (measurement - parameters), thus we need
      // a minus sign below
      alignState.alignmentToResidualDerivative.block(
          iMeasurement, iSurface * eAlignmentSize, measdim, eAlignmentSize) =
          -H * (alignToBound);
    }

    // (e) Extract and fill the track parameters covariance matrix for only
    // measurement states
    // @Todo: add helper function to select rows/columns of a matrix
    for (unsigned int iColState = 0; iColState < measurementStates.size();
         iColState++) {
      size_t colStateIndex = measurementStates.at(iColState).first;
      // Retrieve the block from the source covariance matrix
      CovMatrix correlation =
          sourceTrackParamsCov.block<eBoundSize, eBoundSize>(
              stateRowIndices.at(rowStateIndex),
              stateRowIndices.at(colStateIndex));
      // Fill the block of the target covariance matrix
      size_t iCol =
          alignState.trackParametersDim - (iColState + 1) * eBoundSize;
      alignState.trackParametersCovariance.block<eBoundSize, eBoundSize>(
          iParams, iCol) = correlation;
    }
  }

  // Calculate the chi2 and chi2 derivatives based on the alignment matrixs
  alignState.chi2 = alignState.residual.transpose() *
                    alignState.measurementCovariance.inverse() *
                    alignState.residual;
  alignState.alignmentToChi2Derivative =
      ActsDynamicVector::Zero(alignState.alignmentDof);
  alignState.alignmentToChi2SecondDerivative =
      ActsDynamicMatrix::Zero(alignState.alignmentDof, alignState.alignmentDof);
  // The covariance of residual
  alignState.residualCovariance = ActsDynamicMatrix::Zero(
      alignState.measurementDim, alignState.measurementDim);
  alignState.residualCovariance = alignState.measurementCovariance -
                                  alignState.projectionMatrix *
                                      alignState.trackParametersCovariance *
                                      alignState.projectionMatrix.transpose();

  alignState.alignmentToChi2Derivative =
      2 * alignState.alignmentToResidualDerivative.transpose() *
      alignState.measurementCovariance.inverse() *
      alignState.residualCovariance *
      alignState.measurementCovariance.inverse() * alignState.residual;
  alignState.alignmentToChi2SecondDerivative =
      2 * alignState.alignmentToResidualDerivative.transpose() *
      alignState.measurementCovariance.inverse() *
      alignState.residualCovariance *
      alignState.measurementCovariance.inverse() *
      alignState.alignmentToResidualDerivative;

  return alignState;
}

}  // namespace detail
}  // namespace ActsAlignment
