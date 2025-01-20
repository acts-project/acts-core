// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/NeuralCalibrator.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"

#include <TFile.h>

namespace detail {

template <typename Array>
std::size_t fillChargeMatrix(Array& arr, const ActsExamples::Cluster& cluster,
                             std::size_t size0 = 7u, std::size_t size1 = 7u) {
  // First, rescale the activations to sum to unity. This promotes
  // numerical stability in the index computation
  double totalAct = 0;
  for (const ActsExamples::Cluster::Cell& cell : cluster.channels) {
    totalAct += cell.activation;
  }
  std::vector<double> weights;
  for (const ActsExamples::Cluster::Cell& cell : cluster.channels) {
    weights.push_back(cell.activation / totalAct);
  }

  double acc0 = 0;
  double acc1 = 0;
  for (std::size_t i = 0; i < cluster.channels.size(); i++) {
    acc0 += cluster.channels.at(i).bin[0] * weights.at(i);
    acc1 += cluster.channels.at(i).bin[1] * weights.at(i);
  }

  // By convention, put the center pixel in the middle cell.
  // Achieved by translating the cluster --> compute the offsets
  int offset0 = static_cast<int>(acc0) - size0 / 2;
  int offset1 = static_cast<int>(acc1) - size1 / 2;

  // Zero the charge matrix first, to guard against leftovers
  arr = Eigen::ArrayXXf::Zero(1, size0 * size1);
  // Fill the matrix
  for (const ActsExamples::Cluster::Cell& cell : cluster.channels) {
    // Translate each pixel
    int iMat = cell.bin[0] - offset0;
    int jMat = cell.bin[1] - offset1;
    if (iMat >= 0 && iMat < static_cast<int>(size0) && jMat >= 0 &&
        jMat < static_cast<int>(size1)) {
      typename Array::Index index = iMat * size0 + jMat;
      if (index < arr.size()) {
        arr(index) = cell.activation;
      }
    }
  }
  return size0 * size1;
}

}  // namespace detail

ActsExamples::NeuralCalibrator::NeuralCalibrator(
    const std::filesystem::path& modelPath, std::size_t nComponents,
    std::vector<std::size_t> volumeIds)
    : m_env(ORT_LOGGING_LEVEL_WARNING, "NeuralCalibrator"),
      m_model(m_env, modelPath.c_str()),
      m_nComponents{nComponents},
      m_volumeIds{std::move(volumeIds)} {}

void ActsExamples::NeuralCalibrator::calibrate(
    const MeasurementContainer& measurements, const ClusterContainer* clusters,
    const Acts::GeometryContext& gctx, const Acts::CalibrationContext& cctx,
    const Acts::SourceLink& sourceLink,
    Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy&
        trackState) const {
  trackState.setUncalibratedSourceLink(Acts::SourceLink{sourceLink});
  const IndexSourceLink& idxSourceLink = sourceLink.get<IndexSourceLink>();
  assert((idxSourceLink.index() < measurements.size()) and
         "Source link index is outside the container bounds");

  if (!Acts::rangeContainsValue(m_volumeIds,
                                idxSourceLink.geometryId().volume())) {
    m_fallback.calibrate(measurements, clusters, gctx, cctx, sourceLink,
                         trackState);
    return;
  }

  Acts::NetworkBatchInput inputBatch(1, m_nInputs);
  auto input = inputBatch(0, Eigen::all);

  // TODO: Matrix size should be configurable perhaps?
  std::size_t matSize0 = 7u;
  std::size_t matSize1 = 7u;
  std::size_t iInput = ::detail::fillChargeMatrix(
      input, (*clusters)[idxSourceLink.index()], matSize0, matSize1);

  input[iInput++] = idxSourceLink.geometryId().volume();
  input[iInput++] = idxSourceLink.geometryId().layer();

  const Acts::Surface& referenceSurface = trackState.referenceSurface();
  auto trackParameters = trackState.parameters();

  const ConstVariableBoundMeasurementProxy measurement =
      measurements.getMeasurement(idxSourceLink.index());

  assert(measurement.contains(Acts::eBoundLoc0) &&
         "Measurement does not contain the required bound loc0");
  assert(measurement.contains(Acts::eBoundLoc1) &&
         "Measurement does not contain the required bound loc1");

  auto boundLoc0 = measurement.indexOf(Acts::eBoundLoc0);
  auto boundLoc1 = measurement.indexOf(Acts::eBoundLoc1);

  Acts::Vector2 localPosition{measurement.parameters()[boundLoc0],
                              measurement.parameters()[boundLoc1]};
  Acts::Vector2 localCovariance{measurement.covariance()(boundLoc0, boundLoc0),
                                measurement.covariance()(boundLoc1, boundLoc1)};

  Acts::Vector3 dir = Acts::makeDirectionFromPhiTheta(
      trackParameters[Acts::eBoundPhi], trackParameters[Acts::eBoundTheta]);
  Acts::Vector3 globalPosition =
      referenceSurface.localToGlobal(gctx, localPosition, dir);

  // Rotation matrix. When applied to global coordinates, they
  // are rotated into the local reference frame of the
  // surface. Note that this such a rotation can be found by
  // inverting a matrix whose columns correspond to the
  // coordinate axes of the local coordinate system.
  Acts::RotationMatrix3 rot =
      referenceSurface.referenceFrame(gctx, globalPosition, dir).inverse();
  std::pair<double, double> angles =
      Acts::VectorHelpers::incidentAngles(dir, rot);

  input[iInput++] = angles.first;
  input[iInput++] = angles.second;
  input[iInput++] = localPosition[0];
  input[iInput++] = localPosition[1];
  input[iInput++] = localCovariance[0];
  input[iInput++] = localCovariance[1];
  if (iInput != m_nInputs) {
    throw std::runtime_error("Expected input size of " +
                             std::to_string(m_nInputs) +
                             ", got: " + std::to_string(iInput));
  }

  // Input is a single row, hence .front()
  std::vector<float> output = m_model.runONNXInference(inputBatch).front();
  // Assuming 2-D measurements, the expected params structure is:
  // [           0,    nComponent[ --> priors
  // [  nComponent,  3*nComponent[ --> means
  // [3*nComponent,  5*nComponent[ --> variances
  std::size_t nParams = 5 * m_nComponents;
  if (output.size() != nParams) {
    throw std::runtime_error("Got output vector of size " +
                             std::to_string(output.size()) +
                             ", expected size " + std::to_string(nParams));
  }

  // Most probable value computation of mixture density
  std::size_t iMax = 0;
  if (m_nComponents > 1) {
    iMax = std::distance(
        output.begin(),
        std::max_element(output.begin(), output.begin() + m_nComponents));
  }
  std::size_t iLoc0 = m_nComponents + iMax * 2;
  std::size_t iVar0 = 3 * m_nComponents + iMax * 2;

  Acts::visit_measurement(measurement.size(), [&](auto N) -> void {
    constexpr std::size_t kMeasurementSize = decltype(N)::value;
    const ConstFixedBoundMeasurementProxy<kMeasurementSize> fixedMeasurement =
        static_cast<ConstFixedBoundMeasurementProxy<kMeasurementSize>>(
            measurement);

    Acts::ActsVector<kMeasurementSize> calibratedParameters =
        fixedMeasurement.parameters();
    Acts::ActsSquareMatrix<kMeasurementSize> calibratedCovariance =
        fixedMeasurement.covariance();

    calibratedParameters[boundLoc0] = output[iLoc0];
    calibratedParameters[boundLoc1] = output[iLoc0 + 1];
    calibratedCovariance(boundLoc0, boundLoc0) = output[iVar0];
    calibratedCovariance(boundLoc1, boundLoc1) = output[iVar0 + 1];

    trackState.allocateCalibrated(calibratedParameters, calibratedCovariance);
    trackState.setProjectorSubspaceIndices(fixedMeasurement.subspaceIndices());
  });
}
