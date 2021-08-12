// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#include <iostream>

namespace Acts {
namespace Nmm {
namespace detail {

//defines a default alignment used for a host memory allocated
static constexpr std::size_t NMM_DEFAULT_HOST_ALIGNMENT{alignof(std::max_align_t)};

//@return if `n` is a power of 2
constexpr bool isPow2(std::size_t n) { return (0 == (n & (n - 1)));} 

//@return if `alignment` is a valid memory alignment
constexpr bool isSupportedAlignment(std::size_t alignment) { return isPow2(alignment);}

//@param[in] v value to align
//@param[in] alignment amount, in bytes, must be a power of 2
//
//@return the aligned value, as one would except
constexpr std::size_t align_up(std::size_t v, std::size_t align_bytes) noexcept {
	//if the alignment is not support, the program will end
	assert(isSupportedAlignment(align_bytes));
	return (v + (align_bytes - 1)) & ~(align_bytes - 1);
}

//@param[in] v value to align
//@param[in] alignment amount, in bytes, must be a power of 2
//
//@return the aligned value, as one would except
constexpr std::size_t align_down(std::size_t v, std::size_t align_bytes) noexcept {
	//if the alignment is not support, the program will end
	assert(isSupportedAlignment(align_bytes));
	return v & ~(align_bytes - 1);
}

//@param[in] v value to check for alignment
//@param[in] alignment amount, in bytes, must be a power of 2
//
//@return true if aligned
constexpr bool is_aligned(std::size_t v, std::size_t align_bytes) noexcept {
	//if the alignment is not support, the program will end
	assert(isSupportedAlignment(align_bytes));
	return v == align_down(v, align_bytes);
}

// This function allocates sufficiente memory to satisfy the requested
// size `bytes` with alignment `alignment` using the unary callable `alloc`
// to allocate memory
//
// Allocations returned from `aligned_allocate` MUST be free by calling
// `aligned_deallocate` with the same arguments for `bytes` and `alignment` with 
// a compatible unary `dealloc` callable capable of freeing the memory returned
// from `alloc`
//
// @tparam Alloc a unary callabe type that allocate memory
// @param[in] bytes the size of the allocation
// @param[in] alignment desired alignment of allocation
// @param[in] alloc Unary callable given a size `n` will allocate at least 
// `n` bytes of host memory
// 
// @return void* pointer into allocation of at least `bytes` with desired 
// `alignment`
template <typename Alloc>
void *alignedAllocate(std::size_t bytes, std::size_t alignment, Alloc alloc) {
	assert(isPow2(alignment));

	std::size_t padded_allocation_size{bytes + alignment + sizeof(std::ptrdiff_t)};
	char *const original = static_cast<char *>(alloc(padded_allocation_size));
	void *aligned{original + sizeof(std::ptrdiff_t)};
	std::align(alignment, bytes, aligned, padded_allocation_size);
	std::ptrdiff_t offset = static_cast<char *>(aligned) - original;

	*(static_cast<std::ptrdiff_t *>(aligned) - 1) = offset;

	return aligned;
}

// This funciton frees an allocation from `aligned_allocate`, so have to be 
// called with the same arguments for `bytes` and `alignment` with a compatible
// unary `dealloc` callable capable of freeing the memory returned
//
// @tparam Dealloc a unary callable type that deallocates memory
// @param[in] p the aligned pointer to deallocate
// @param[in] bytes the number of bytes requested from `aligned_allocate`
// @param[in] alignment the alignment required from `aligned_allocate`
// @param[in] dealloc a unvary callable capable of freeing memory returned from
// `alloc` in `aligned_allocate` 
template <typename Dealloc>
//void alignedDeallocate(void *p, std::size_t bytes, std::size_t alignment, Dealloc dealloc){
void alignedDeallocate(void *p, std::size_t alignment, Dealloc dealloc){
	(void)alignment;

	std::ptrdiff_t const offset = *(reinterpret_cast<std::ptrdiff_t *>(p) - 1);

	void *const original = static_cast<char *>(p) - offset;

	dealloc(original);
}

} //namespace detail
} //namespace Nmm
} //namespace Acts