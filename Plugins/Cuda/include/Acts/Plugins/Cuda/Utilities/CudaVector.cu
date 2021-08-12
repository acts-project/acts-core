// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Cuda/Utilities/CpuVector.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager.hpp"

#include <iostream>
#include <memory>

#include "CudaUtils.cu"
#include "cuda.h"
#include "cuda_runtime.h"

namespace Acts {

template <typename var_t>
class CpuVector;

template <typename var_t>
class CudaVector {
 public:
  CudaVector() = delete;
  CudaVector(size_t size) {
    m_size = size;
    //ACTS_CUDA_ERROR_CHECK(
    //    cudaMalloc((var_t**)&m_devPtr, m_size * sizeof(var_t)));
    m_devPtr = static_cast< var_t* >(getGPUmmr()->allocate(m_size * sizeof(var_t)));
  }

  CudaVector(size_t size, var_t* vector) {
    m_size = size;
    //ACTS_CUDA_ERROR_CHECK(
    //    cudaMalloc((var_t**)&m_devPtr, m_size * sizeof(var_t)));
    m_devPtr = static_cast< var_t* >(getGPUmmr()->allocate(m_size * sizeof(var_t)));
    copyH2D(vector, m_size, 0);
  }

  CudaVector(size_t size, var_t* vector, size_t len, size_t offset) {
    m_size = size;
    //ACTS_CUDA_ERROR_CHECK(
    //    cudaMalloc((var_t**)&m_devPtr, m_size * sizeof(var_t)));
    m_devPtr = static_cast< var_t* >(getGPUmmr()->allocate(m_size * sizeof(var_t)));
    copyH2D(vector, len, offset);
  }

  ~CudaVector() {
    if (m_devPtr){
      //cudaFree(m_devPtr);
      getGPUmmr()->deallocate(static_cast< void* >(m_devPtr), m_size);
    }
  }

  var_t* get(size_t offset = 0) { return m_devPtr + offset; }

  void copyH2D(var_t* vector, size_t len, size_t offset) {
    ACTS_CUDA_ERROR_CHECK(cudaMemcpy(m_devPtr + offset, vector,
                                     len * sizeof(var_t),
                                     cudaMemcpyHostToDevice));
  }
  void copyH2D(var_t* vector, size_t len, size_t offset, cudaStream_t* stream) {
    ACTS_CUDA_ERROR_CHECK(cudaMemcpyAsync(m_devPtr + offset, vector,
                                          len * sizeof(var_t),
                                          cudaMemcpyHostToDevice, *stream));
  }

  void zeros() {
    ACTS_CUDA_ERROR_CHECK(cudaMemset(m_devPtr, 0, m_size * sizeof(var_t)));
  }

 private:
  var_t* m_devPtr = nullptr;
  size_t m_size;
};
}  // namespace Acts
