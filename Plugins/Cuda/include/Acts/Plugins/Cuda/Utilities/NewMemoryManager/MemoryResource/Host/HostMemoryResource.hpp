// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// System include(s)
#include <cstddef>
#include <utility>

namespace Acts {
namespace Nmm {
namespace MemoryResource {


// @brief Base class for host memory allocation.
// 
// This is based on `std::pmr::memory_resource`:
// https:// en.cppreference.com/w/cpp/memory/memory_resource
// 
// This class serves as the interface that all host memory resource
// implementations must satisfy.
// 
// There are two private, pure virtual functions that all derived classes must
// implement: `doAllocate` and `doDeallocate`. Optionally, derived classes may
// also override `isEqual`. By default, `isEqual` simply performs an identity
// comparison.
// 
// The public, non-virtual functions `allocate`, `deallocate`, and `isEqual`
// simply call the private virtual functions. The reason for this is to allow
// implementing shared, default behavior in the base class. For example, the
// base class' `allocate` function may log every allocation, no matter what
// derived class implementation is used.

class HostMemoryResource {
 public:
  virtual ~HostMemoryResource() = default;

  // Allocates memory on the host of size at least `bytes` bytes.
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
  void* allocate(std::size_t bytes, std::size_t alignment = alignof(std::max_align_t)) {
    return doAllocate(bytes, alignment);
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
  //void deallocate(void* p, std::size_t bytes, std::size_t alignment = alignof(std::max_align_t)) {
  //  doDeallocate(p, bytes, alignment);
  //}

  void deallocate(void* p, std::size_t alignment = alignof(std::max_align_t)) {
    doDeallocate(p, alignment);
  }

  // Compare this resource to another.
  // 
  // Two `HostMemoryResource`s compare equal if and only if memory allocated
  // from one `HostMemoryResource` can be deallocated from the other and vice
  // versa.
  // 
  // By default, simply checks if \p *this and \p other refer to the same
  // object, i.e., does not check if they are two objects of the same class.
  // 
  // @param other The other resource to compare to
  // @returns If the two resources are equivalent
  bool isEqual(HostMemoryResource const& other) const noexcept { return doIsEqual(other); }

 private:
  // Allocates memory on the host of size at least `bytes` bytes.
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
  virtual void* doAllocate(std::size_t bytes,
                            std::size_t alignment = alignof(std::max_align_t)) = 0;

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
  //virtual void doDeallocate(void* p, std::size_t bytes, std::size_t alignment = alignof(std::max_align_t)) = 0;
  virtual void doDeallocate(void* p, std::size_t alignment = alignof(std::max_align_t)) = 0;

  // Compare this resource to another.
  // 
  // Two HostMemoryResources compare equal if and only if memory allocated
  // from one HostMemoryResource can be deallocated from the other and vice
  // versa.
  // 
  // By default, simply checks if \p *this and \p other refer to the same
  // object, i.e., does not check if they are two objects of the same class.
  // 
  // @param other The other resource to compare to
  // @return true If the two resources are equivalent
  // @return false If the two resources are not equal
  virtual bool doIsEqual(HostMemoryResource const& other) const noexcept
  {
    return this == &other;
  }
};// class HostMemoryResource
} // namespace MemoryResource
} // namespace Nmm
} // namespace Acts