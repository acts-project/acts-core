// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <unsupported/Eigen/Splines>

namespace Acts::Interpolation3D {

/// @brief Helper function to interpolate points using a spline
/// from Eigen
///
/// The only requirement is that the input trajectory type has
/// a method empty() and size() and that the elements can be
/// accessed with operator[] and have themselves a operator[] to
/// access the coordinates.
///
/// @tparam input_trajectory_type input trajectory type
///
/// @param inputsRaw input vector points
/// @param nPoints number of interpolation points
/// @param keepOriginalHits keep the original hits in the trajectory
///
/// @return std::vector<Acts::Vector3> interpolated points
template <typename trajectory_type>
trajectory_type spline(const trajectory_type& inputsRaw, std::size_t nPoints,
                       bool keepOriginalHits = false) {
  trajectory_type output;
  if (inputsRaw.empty()) {
    return output;
  }

  using InputVectorType = typename trajectory_type::value_type;

  std::vector<Vector3> inputs;
  // If input type is a vector of Vector3 we can use it directly
  if constexpr (std::is_same_v<trajectory_type, std::vector<Vector3>>) {
    inputs = inputsRaw;
  } else {
    inputs.reserve(inputsRaw.size());
    for (const auto& input : inputsRaw) {
      inputs.push_back(Vector3(input[0], input[1], input[2]));
    }
  }

  // Don't do anything if we have less than 3 points or less interpolation
  // points than input points
  if (inputsRaw.size() < 3 || nPoints <= inputsRaw.size()) {
    return inputsRaw;
  } else {
    Eigen::MatrixXd points(3, inputs.size());
    for (std::size_t i = 0; i < inputs.size(); ++i) {
      points.col(i) = inputs[i].transpose();
    }

    // MARK: fpeMaskBegin(FLTDIV, 1, #4024)
    // MARK: fpeMaskBegin(FLTINV, 1, #4024)
    Eigen::Spline<double, 3> spline3D =
        Eigen::SplineFitting<Eigen::Spline<double, 3>>::Interpolate(points, 2);
    // MARK: fpeMaskEnd(FLTDIV)
    // MARK: fpeMaskEnd(FLTINV)

    double step = 1. / (nPoints - 1);
    for (std::size_t i = 0; i < nPoints; ++i) {
      double t = i * step;
      InputVectorType point;
      point[0] = spline3D(t)[0];
      point[1] = spline3D(t)[1];
      point[2] = spline3D(t)[2];
      output.push_back(point);
    }
  }
  // If we want to keep the original hits, we add them to the output
  // (first and last are there anyway)
  if (keepOriginalHits) {
    output.insert(output.begin(), inputsRaw.begin() + 1, inputsRaw.end() - 1);
    // We need to sort the output in distance to first
    std::sort(output.begin(), output.end(),
              [&inputs](const auto& a, const auto& b) {
                const auto ifront = inputs.front();
                double da2 = (a[0] - ifront[0]) * (a[0] - ifront[0]) +
                             (a[1] - ifront[1]) * (a[1] - ifront[1]) +
                             (a[2] - ifront[2]) * (a[2] - ifront[2]);
                double db2 = (b[0] - ifront[0]) * (b[0] - ifront[0]) +
                             (b[1] - ifront[1]) * (b[1] - ifront[1]) +
                             (b[2] - ifront[2]) * (b[2] - ifront[2]);
                return da2 < db2;
              });
  }

  return output;
}

}  // namespace Acts::Interpolation3D
