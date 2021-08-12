// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Cuda/Utilities/CudaScalar.cu"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager.hpp"

namespace Acts {

template <typename var_t>
class CudaScalar;

template <typename var_t>
class CpuScalar {
 public:
  CpuScalar(bool pinned = 0) {
    m_pinned = pinned;
    if (pinned == 0) {
      m_hostPtr = new var_t[1];
    } else if (pinned == 1) {
      // cudaMallocHost(&m_hostPtr, sizeof(var_t));
      m_hostPtr = static_cast<var_t*>(Acts::getPinnedmmr()->allocate(sizeof(var_t)));
    }
  }

  CpuScalar(CudaScalar<var_t>* cuScalar, bool pinned = 0) {
    m_pinned = pinned;
    if (pinned == 0) {
      m_hostPtr = new var_t[1];
    } else if (pinned == 1) {
      // cudaMallocHost(&m_hostPtr, sizeof(var_t));
      m_hostPtr = static_cast<var_t*>(Acts::getPinnedmmr()->allocate(sizeof(var_t)));
    }
    cudaMemcpy(m_hostPtr, cuScalar->get(), sizeof(var_t),
               cudaMemcpyDeviceToHost);
  }

  ~CpuScalar() {
    if (!m_pinned) {
      delete m_hostPtr;
    } else if (m_pinned && m_hostPtr) {
      // cudaFreeHost(m_hostPtr);
      Acts::getPinnedmmr()->deallocate(m_hostPtr, m_size);
    }
  }

  var_t* get() { return m_hostPtr; }

  void Set(var_t val) { m_hostPtr[0] = val; }

 private:
  var_t* m_hostPtr = nullptr;
  size_t m_size;
  bool m_pinned;
};

}  // namespace Acts
