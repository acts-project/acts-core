// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/TensorRTEdgeClassifier.hpp"

#include "Acts/Plugins/ExaTrkX/detail/Utils.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>

#include <NvInfer.h>
#include <NvInferPlugin.h>
#include <NvInferRuntimeBase.h>
#include <c10/cuda/CUDAGuard.h>
#include <cuda_runtime.h>

#include "printCudaMemInfo.hpp"

using namespace torch::indexing;

namespace {

class TensorRTLogger : public nvinfer1::ILogger {
  std::unique_ptr<const Acts::Logger> m_logger;

 public:
  TensorRTLogger(Acts::Logging::Level lvl)
      : m_logger(Acts::getDefaultLogger("TensorRT", lvl)) {}

  void log(Severity severity, const char *msg) noexcept override {
    const auto &logger = *m_logger;
    switch (severity) {
      case Severity::kVERBOSE:
        ACTS_DEBUG(msg);
        break;
      case Severity::kINFO:
        ACTS_INFO(msg);
        break;
      case Severity::kWARNING:
        ACTS_WARNING(msg);
        break;
      case Severity::kERROR:
        ACTS_ERROR(msg);
        break;
      case Severity::kINTERNAL_ERROR:
        ACTS_FATAL(msg);
        break;
    }
  }
};

}  // namespace

namespace Acts {

TensorRTEdgeClassifier::TensorRTEdgeClassifier(
    const Config &cfg, std::unique_ptr<const Logger> _logger)
    : m_logger(std::move(_logger)),
      m_cfg(cfg),
      m_trtLogger(std::make_unique<TensorRTLogger>(m_logger->level())) {
  auto status = initLibNvInferPlugins(m_trtLogger.get(), "");
  assert(status);

  std::size_t fsize =
      std::filesystem::file_size(std::filesystem::path(m_cfg.modelPath));
  std::vector<char> engineData(fsize);

  ACTS_DEBUG("Load '" << m_cfg.modelPath << "' with size " << fsize);

  std::ifstream engineFile(m_cfg.modelPath);
  engineFile.read(engineData.data(), fsize);

  m_runtime.reset(nvinfer1::createInferRuntime(*m_trtLogger));

  m_engine.reset(m_runtime->deserializeCudaEngine(engineData.data(), fsize));

  m_context.reset(m_engine->createExecutionContext());
}

TensorRTEdgeClassifier::~TensorRTEdgeClassifier() {}

auto milliseconds = [](const auto &a, const auto &b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
};

struct TimePrinter {
  const char *name;
  decltype(std::chrono::high_resolution_clock::now()) t0, t1;
  TimePrinter(const char *n) : name(n) {
    t0 = std::chrono::high_resolution_clock::now();
  }
  ~TimePrinter() {
    std::cout << name << ": " << milliseconds(t0, t1) << std::endl;
  }
};

#if 0
#define TIME_BEGIN(name) TimePrinter printer##name(#name);
#define TIME_END(name) \
  printer##name.t1 = std::chrono::high_resolution_clock::now();
#else
#define TIME_BEGIN(name) /*nothing*/
#define TIME_END(name)   /*ǹothing*/
#endif

std::tuple<std::any, std::any, std::any, std::any>
TensorRTEdgeClassifier::operator()(std::any inNodeFeatures,
                                   std::any inEdgeIndex,
                                   std::any inEdgeFeatures,
                                   const ExecutionContext &execContext) {
  decltype(std::chrono::high_resolution_clock::now()) t0, t1, t2, t3, t4, t5;
  t0 = std::chrono::high_resolution_clock::now();

  c10::cuda::CUDAStreamGuard(execContext.stream.value());

  auto nodeFeatures =
      std::any_cast<torch::Tensor>(inNodeFeatures).to(torch::kCUDA);

  auto edgeIndex = std::any_cast<torch::Tensor>(inEdgeIndex).to(torch::kCUDA);
  ACTS_DEBUG("edgeIndex: " << detail::TensorDetails{edgeIndex});

  auto edgeFeatures =
      std::any_cast<torch::Tensor>(inEdgeFeatures).to(torch::kCUDA);
  ACTS_DEBUG("edgeFeatures: " << detail::TensorDetails{edgeFeatures});

  t1 = std::chrono::high_resolution_clock::now();

  m_context->setInputShape(
      "x", nvinfer1::Dims2{nodeFeatures.size(0), nodeFeatures.size(1)});
  m_context->setTensorAddress("x", nodeFeatures.data_ptr());

  m_context->setInputShape(
      "edge_index", nvinfer1::Dims2{edgeIndex.size(0), edgeIndex.size(1)});
  m_context->setTensorAddress("edge_index", edgeIndex.data_ptr());

  m_context->setInputShape(
      "edge_attr", nvinfer1::Dims2{edgeFeatures.size(0), edgeFeatures.size(1)});
  m_context->setTensorAddress("edge_attr", edgeFeatures.data_ptr());

  void *outputMem{nullptr};
  std::size_t outputSize = edgeIndex.size(1) * sizeof(float);
  cudaMalloc(&outputMem, outputSize);
  m_context->setTensorAddress("output", outputMem);

  t2 = std::chrono::high_resolution_clock::now();

  {
    auto stream = execContext.stream.value().stream();
    auto status = m_context->enqueueV3(stream);
    cudaStreamSynchronize(stream);
    ACTS_VERBOSE("TensorRT output status: " << std::boolalpha << status);
  }

  t3 = std::chrono::high_resolution_clock::now();

  auto scores = torch::from_blob(
      outputMem, edgeIndex.size(1), 1, [](void *ptr) { cudaFree(ptr); },
      torch::TensorOptions().device(torch::kCUDA).dtype(torch::kFloat32));

  scores.sigmoid_();

  ACTS_VERBOSE("Size after classifier: " << scores.size(0));
  ACTS_VERBOSE("Slice of classified output:\n"
               << scores.slice(/*dim=*/0, /*start=*/0, /*end=*/9));
  printCudaMemInfo(logger());

  torch::Tensor mask = scores > m_cfg.cut;
  torch::Tensor edgesAfterCut = edgeIndex.index({Slice(), mask});

  scores = scores.masked_select(mask);
  ACTS_VERBOSE("Size after score cut: " << edgesAfterCut.size(1));
  printCudaMemInfo(logger());

  t4 = std::chrono::high_resolution_clock::now();

  ACTS_DEBUG("Time anycast:  " << milliseconds(t0, t1));
  ACTS_DEBUG("Time alloc, set shape " << milliseconds(t1, t2));
  ACTS_DEBUG("Time inference:       " << milliseconds(t2, t3));
  ACTS_DEBUG("Time sigmoid and cut: " << milliseconds(t3, t4));

  return {nodeFeatures, edgesAfterCut, edgeFeatures, scores};
}

}  // namespace Acts
