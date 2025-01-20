// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/ExaTrkX/Stages.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

#include <torch/torch.h>

namespace nvinfer1 {
class IRuntime;
class ICudaEngine;
class ILogger;
class IExecutionContext;
}  // namespace nvinfer1

namespace Acts {

class TensorRTEdgeClassifier final : public Acts::EdgeClassificationBase {
 public:
  struct Config {
    std::string modelPath;
    std::vector<int> selectedFeatures = {};
    float cut = 0.21;
    int deviceID = 0;
    bool useEdgeFeatures = false;
    bool doSigmoid = true;
  };

  TensorRTEdgeClassifier(const Config &cfg,
                         std::unique_ptr<const Logger> logger);
  ~TensorRTEdgeClassifier();

  std::tuple<std::any, std::any, std::any, std::any> operator()(
      std::any nodeFeatures, std::any edgeIndex, std::any edgeFeatures = {},
      const ExecutionContext &execContext = {}) override;

  Config config() const { return m_cfg; }
  torch::Device device() const override { return torch::kCUDA; };

 private:
  std::unique_ptr<const Acts::Logger> m_logger;
  const auto &logger() const { return *m_logger; }

  Config m_cfg;

  std::unique_ptr<nvinfer1::IRuntime> m_runtime;
  std::unique_ptr<nvinfer1::ICudaEngine> m_engine;
  std::unique_ptr<nvinfer1::ILogger> m_trtLogger;
  std::unique_ptr<nvinfer1::IExecutionContext> m_context;
};

}  // namespace Acts
