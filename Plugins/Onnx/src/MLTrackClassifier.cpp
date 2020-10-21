// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cassert>
#include <stdexcept>

#include "Acts/Plugins/Onnx/MLTrackClassifier.hpp"

// prediction function
Acts::MLTrackClassifier::TrackLabels Acts::MLTrackClassifier::predictTrackLabel(
    std::vector<float>& inputFeatures, const double& decisionThreshProb) const {
  // check that the decision threshold is a probability
  if ((decisionThreshProb < 0.) || (decisionThreshProb > 1.)) {
    throw std::invalid_argument(
        "predictTrackLabel: Decision threshold "
        "probability is not in [0, 1].");
  }

  // run the model over the input
  float* outputTensor = runONNXInference(inputFeatures);
  // this is binary classification, so only need first value
  float outputProbability = outputTensor[0];

  // the output layer computes how confident the network is that the track is a
  // duplicate, so need to convert that to a label
  if (outputProbability > decisionThreshProb) {
    return TrackLabels::duplicate;
  }
  return TrackLabels::good;
}

// function that checks if the predicted track label is duplicate
bool Acts::MLTrackClassifier::isDuplicate(
    std::vector<float>& inputFeatures, const double& decisionThreshProb) const {
  Acts::MLTrackClassifier::TrackLabels predictedLabel =
      Acts::MLTrackClassifier::predictTrackLabel(inputFeatures,
                                                 decisionThreshProb);
  return predictedLabel == Acts::MLTrackClassifier::TrackLabels::duplicate;
}
