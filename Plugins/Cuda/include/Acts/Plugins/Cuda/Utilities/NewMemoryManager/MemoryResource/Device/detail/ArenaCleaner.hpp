// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/detail/Aligned.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/CudaStreamView.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/detail/Arena.hpp"

// CUDA include(s).
#include <cuda_runtime_api.h>

// System include(s).
#include <algorithm>
#include <limits>
#include <memory>
#include <mutex>
#include <set>
#include <unordered_map>

namespace Acts {
namespace Nmm {
namespace MemoryResource {
namespace detail {
namespace Arena {

template <typename Upstream>
class ArenaCleaner {
 public:
  explicit ArenaCleaner(std::shared_ptr<Arena<Upstream>> const& a) : arena_(a) {}

  // Disable copy (and move) semantics.
  ArenaCleaner(const ArenaCleaner&) = delete;
  ArenaCleaner& operator=(const ArenaCleaner&) = delete;

  ~ArenaCleaner()
  {
    if (!arena_.expired()) {
      auto arenaPtr = arena_.lock();
      arenaPtr->clean();
    }
  }

 private:
  /// A non-owning pointer to the Arena that may need cleaning.
  std::weak_ptr<Arena<Upstream>> arena_;
};// class ArenaCleaner

} // namaspace Arena
} // namaspace detail
} // namaspace MemoryResource
} // namaspace Nmm
} // namaspace Acts