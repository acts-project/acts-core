// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <ActsExamples/EventData/NeuralCalibrator.hpp>

#include <TFile.h>

namespace detail {

// TODO make size configurable
std::array<float, 7 * 7> chargeMatrix(const ActsExamples::Cluster& cl) {
  // Compute the centroid index
  double acc0 = 0;
  double acc1 = 0;
  double tot = 0;
  for (const ActsExamples::Cluster::Cell& cell : cl.channels) {
    acc0 += cell.bin[0] * cell.activation;
    acc1 += cell.bin[1] * cell.activation;
    tot += cell.activation;
  }
  int icntr0 = static_cast<int>(acc0 / tot);
  int icntr1 = static_cast<int>(acc1 / tot);

  // By convention, put centroid in (3, 3) cell
  int diff0 = icntr0 - 3;
  int diff1 = icntr1 - 3;

  // Fill the matrix
  std::array<float, 7 * 7> mat{};
  for (const ActsExamples::Cluster::Cell& cell : cl.channels) {
    int im = cell.bin[0] - diff0;
    int jm = cell.bin[1] - diff1;
    if (im >= 0 && im < 7 && jm >= 0 && jm < 7)
      mat.at(im * 7 + jm) = cell.activation;
  }

  return mat;
}

}  // namespace detail

ActsExamples::NeuralCalibrator::NeuralCalibrator(
    const std::filesystem::path& modelPath,
    const std::filesystem::path& normalizationPath)
    : m_env(ORT_LOGGING_LEVEL_WARNING, "NeuralCalibrator"),
      m_model(m_env, modelPath.c_str()) {
  readNormalization(normalizationPath);
}

void ActsExamples::NeuralCalibrator::readNormalization(
    const std::filesystem::path& normalizationPath) {
  TFile ifile(normalizationPath.c_str(), "READ");
  if (ifile.IsZombie()) {
    throw std::runtime_error("Unable to read input normalization from " +
                             normalizationPath.string());
  }

  std::vector<float>* offsets = ifile.Get<std::vector<float>>("offsets");
  if (offsets == nullptr || offsets->size() != n_inputs) {
    throw std::runtime_error("Unable to normalization constants (offsets)");
  }

  std::vector<float>* scales = ifile.Get<std::vector<float>>("scales");
  if (scales == nullptr || scales->size() != n_inputs) {
    throw std::runtime_error("Unable to normalization constants (scales)");
  }

  m_offsets = *offsets;
  m_scales = *scales;
}

void ActsExamples::NeuralCalibrator::calibrate(
    const MeasurementContainer& measurements, const ClusterContainer* clusters,
    const Acts::GeometryContext& /*gctx*/,
    Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy&
        trackState) const {
  const IndexSourceLink& sourceLink =
      trackState.getUncalibratedSourceLink().get<IndexSourceLink>();
  assert((sourceLink.index() < measurements.size()) and
         "Source link index is outside the container bounds");

  std::array<float, 7 * 7> matrix =
      ::detail::chargeMatrix((*clusters)[sourceLink.index()]);
  std::vector<float> input;
  input.assign(matrix.begin(), matrix.end());
  input.push_back(sourceLink.geometryId().volume());
  input.push_back(sourceLink.geometryId().layer());
  input.push_back(trackState.parameters()[Acts::eBoundPhi]);
  input.push_back(trackState.parameters()[Acts::eBoundTheta]);

  std::visit(
      [&](const auto& meas) {
        // Evaluate the MDN
        auto E = meas.expander();
        auto P = meas.projector();
        Acts::ActsVector<Acts::eBoundSize> fpar = E * meas.parameters();

        input.push_back(fpar[Acts::eBoundLoc0]);
        input.push_back(fpar[Acts::eBoundLoc1]);
        assert(input.size() == n_inputs && "Missing some inputs!");

        // normalize
        for (size_t i = 0; i < n_inputs; i++) {
          input[i] += m_offsets[i];
          input[i] *= m_scales[i];
        }

        std::vector<float> output = m_model.runONNXInference(input);
        if (output.size() != 4) {
          throw std::runtime_error("Expected output size of 4, got: " +
                                   std::to_string(output.size()));
        }

        // Save the calibrated positions

        Acts::ActsSymMatrix<Acts::eBoundSize> fcov =
            E * meas.covariance() * E.transpose();

        fpar[Acts::eBoundLoc0] = output[0];
        fpar[Acts::eBoundLoc1] = output[1];
        fcov(Acts::eBoundLoc0, Acts::eBoundLoc0) = 1.0 / output[2];
        fcov(Acts::eBoundLoc1, Acts::eBoundLoc1) = 1.0 / output[3];

        constexpr size_t kSize = std::decay_t<decltype(meas)>::size();
        std::array<Acts::BoundIndices, kSize> indices = meas.indices();
        Acts::ActsVector<kSize> cpar = P * fpar;
        Acts::ActsSymMatrix<kSize> ccov = P * fcov * P.transpose();

        Acts::SourceLink sl{sourceLink.geometryId(), sourceLink};

        Acts::Measurement<Acts::BoundIndices, kSize> cmeas(std::move(sl),
                                                           indices, cpar, ccov);

        trackState.allocateCalibrated(cmeas.size());
        trackState.setCalibrated(cmeas);
      },
      (measurements)[sourceLink.index()]);
}
