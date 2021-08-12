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
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Host/HostMemoryResource.hpp"

// System include(s)
#include <cstddef>
#include <utility>
#include <iostream>

namespace Acts {
namespace Nmm {
namespace MemoryResource {


// @brief A `HostMemoryResource` that uses `cudaMallocHost` to allocate
// pinned/page-locked host memory.
// See https://devblogs.nvidia.com/how-optimize-data-transfers-cuda-cc/
class CPUMemoryResource final : public HostMemoryResource {
 public:
  CPUMemoryResource()                               = default;
  ~CPUMemoryResource()                              = default;
  CPUMemoryResource(CPUMemoryResource const &) = default;
  CPUMemoryResource(CPUMemoryResource &&)      = default;
  CPUMemoryResource &operator=(CPUMemoryResource const &) = default;
  CPUMemoryResource &operator=(CPUMemoryResource &&) = default;

 private:
  // Allocates memory with malloc on the host of size at least `bytes` bytes.
  //
  // The returned storage is aligned to the specified `alignment` if supported,
  // and to `alignof(std::max_align_t)` otherwise.
  //
  // @throws std::bad_alloc When the requested `bytes` and `alignment` cannot be
  // allocated.
  //
  // @param bytes The size of the allocation
  // @param alignment Alignment of the allocation
  // @return void* Pointer to the newly allocated memory
  void *doAllocate(std::size_t bytes, std::size_t alignment = alignof(std::max_align_t)) override {
    // don't allocate anything if the user requested zero bytes
    if (0 == bytes) { return nullptr; }

    // If the requested alignment isn't supported, use default
    alignment = (Nmm::detail::isSupportedAlignment(alignment)) ? alignment : alignof(std::max_align_t);
    return Nmm::detail::alignedAllocate(bytes, alignment, [](std::size_t size) {
      void *p = malloc(size);
      if (nullptr == p) { throw std::bad_alloc{}; }
      return p;
    });
  }

  // Deallocate memory pointed to by `p`.
  //
  // `p` must have been returned by a prior call to `allocate(bytes,alignment)`
  // on a `HostMemoryResource` that compares equal to `*this`, and the storage
  // it points to must not yet have been deallocated, otherwise behavior is
  // undefined.
  //
  // @throws Nothing.
  //
  // @param p Pointer to be deallocated
  // @param bytes The size in bytes of the allocation. This must be equal to the
  // value of `bytes` that was passed to the `allocate` call that returned `p`.
  // @param alignment Alignment of the allocation. This must be equal to the
  // value of `alignment` that was passed to the `allocate` call that returned
  // `p`.
  // @param stream Stream on which to perform deallocation
  //void doDeallocate(void *p, std::size_t bytes, std::size_t alignment = alignof(std::max_align_t)) override {
  void doDeallocate(void *p, std::size_t alignment = alignof(std::max_align_t)) override {
    (void)alignment;
    if (nullptr == p) { return; }
    Nmm::detail::alignedDeallocate(
      p, alignment, [](void *q) { free(q); });
      //p, bytes, alignment, [](void *q) { free(q); });
  }
};// class CPUMemoryResource
} // namespace MemoryResource
} // namespace Nmm
} // namespace Acts
