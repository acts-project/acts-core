// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Cuda/Utilities/CudaVector.cu"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager.hpp"

#include <cstring>

namespace Acts {

template <typename var_t>
class CudaVector;

template <typename var_t>
class CpuVector {
 public:
  CpuVector() = delete;
  CpuVector(size_t size, bool pinned = 0) {
    m_size = size;
    m_pinned = pinned;
    if (pinned == 0) {
      m_hostPtr = new var_t[m_size];
    } else if (pinned == 1) {
      //cudaMallocHost(&m_hostPtr, m_size * sizeof(var_t));
      m_hostPtr = static_cast<var_t*>(Acts::getPinnedmmr()->allocate(m_size * sizeof(var_t)));
    }
  }

  CpuVector(size_t size, CudaVector<var_t>* cuVec, bool pinned = 0) {
    m_size = size;
    m_pinned = pinned;
    if (pinned == 0) {
      m_hostPtr = new var_t[m_size];
    } else if (pinned == 1) {
      // cudaMallocHost(&m_hostPtr, m_size * sizeof(var_t));
      m_hostPtr = static_cast<var_t*>(Acts::getPinnedmmr()->allocate(m_size * sizeof(var_t)));
    }
    cudaMemcpy(m_hostPtr, cuVec->get(), m_size * sizeof(var_t),
               cudaMemcpyDeviceToHost);
  }

  ~CpuVector() {
    if (!m_pinned) {
      delete m_hostPtr;
    } else if (m_pinned && m_hostPtr) {
      // cudaFreeHost(m_hostPtr);
      Acts::getPinnedmmr()->deallocate(m_hostPtr, m_size);
    }
  }

  var_t* get(size_t offset = 0) { return m_hostPtr + offset; }

  void set(size_t offset, var_t val) { m_hostPtr[offset] = val; }

  void copyD2H(var_t* devPtr, size_t len, size_t offset) {
    cudaMemcpy(m_hostPtr + offset, devPtr, len * sizeof(var_t),
               cudaMemcpyDeviceToHost);
  }

  void copyD2H(var_t* devPtr, size_t len, size_t offset, cudaStream_t* stream) {
    cudaMemcpyAsync(m_hostPtr + offset, devPtr, len * sizeof(var_t),
                    cudaMemcpyDeviceToHost, *stream);
  }

  void zeros() { memset(m_hostPtr, 0, m_size * sizeof(var_t)); }

 private:
  var_t* m_hostPtr = nullptr;
  size_t m_size;
  bool m_pinned;
};

}  // namespace Acts
