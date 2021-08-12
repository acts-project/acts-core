// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Cuda/Utilities/CudaMatrix.cu"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager.hpp"

// column-major style Matrix Definition

namespace Acts {

template <typename var_t>
class CudaMatrix;

template <typename var_t>
class CpuMatrix {
 public:
  CpuMatrix() = delete;
  CpuMatrix(size_t nRows, size_t nCols, bool pinned = 0) {
    m_setSize(nRows, nCols);
    m_pinned = pinned;
    if (pinned == 0) {
      m_hostPtr = new var_t[m_size];
    } else if (pinned == 1) {
      // cudaMallocHost(&m_hostPtr, m_size * sizeof(var_t));
      m_hostPtr = static_cast<var_t*>(Acts::getPinnedmmr()->allocate(m_size * sizeof(var_t)));
    }
  }

  CpuMatrix(size_t nRows, size_t nCols, CudaMatrix<var_t>* cuMat,
            bool pinned = 0) {
    m_setSize(nRows, nCols);
    m_pinned = pinned;
    if (pinned == 0) {
      m_hostPtr = new var_t[m_size];
    } else if (pinned == 1) {
      // cudaMallocHost(&m_hostPtr, m_nRows * m_nCols * sizeof(var_t));
      m_hostPtr = static_cast<var_t*>(Acts::getPinnedmmr()->allocate(m_nRows * m_nCols * sizeof(var_t)));
    }
    cudaMemcpy(m_hostPtr, cuMat->get(0, 0), m_size * sizeof(var_t),
               cudaMemcpyDeviceToHost);
  }

  ~CpuMatrix() {
    if (!m_pinned) {
      delete m_hostPtr;
    } else if (m_pinned && m_hostPtr) {
      //cudaFreeHost(m_hostPtr);
      Acts::getPinnedmmr()->deallocate(m_hostPtr, m_size);
    }
  }

  var_t* get(size_t row = 0, size_t col = 0) {
    size_t offset = row + col * m_nRows;
    return m_hostPtr + offset;
  }

  void set(size_t row, size_t col, var_t val) {
    size_t offset = row + col * m_nRows;
    m_hostPtr[offset] = val;
  }

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
  size_t m_nCols;
  size_t m_nRows;
  size_t m_size;
  bool m_pinned;

  void m_setSize(size_t row, size_t col) {
    m_nRows = row;
    m_nCols = col;
    m_size = m_nRows * m_nCols;
  }
};

}  // namespace Acts
