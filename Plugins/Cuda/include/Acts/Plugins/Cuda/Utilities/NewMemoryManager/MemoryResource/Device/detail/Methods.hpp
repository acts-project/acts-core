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
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/detail/Block.hpp"

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

// the required allocation alignment
constexpr std::size_t allocationAlignment = 256;

// align up to the allocation alignment
//
// @param[in] value value to align
// @return the aligned value
constexpr std::size_t alignUp(std::size_t value) noexcept {
	return Nmm::detail::align_up(value, allocationAlignment);
} 

// align down to the allocation alignment
//
// @param[in] value value to align
// @return the aligned value
constexpr std::size_t alignDown(std::size_t value) noexcept {
	return Nmm::detail::align_down(value, allocationAlignment);
}

// get the first free block of at least `size` bytes
//
// @param[in] freeBlocks the adress-ordered set of free blocks
// @param[in] size the number of bytes to allocate
// @return a block of memory of at least `size` bytes, or an empty if
// not found
inline Block firstFit(std::set<Block>& freeBlocks, std::size_t size){
	auto const iter = std::find_if(freeBlocks.cbegin(), freeBlocks.cend(), [size](auto const& b) { return b.fits(size); });

	if(iter == freeBlocks.cend()) {
		return {};
	} else {
		// remove the block from the freeList
		auto const b = *iter;
		auto const i = freeBlocks.erase(iter);

		if(b.size() > size) {
			// split the block and put the remainder back.
			auto const split = b.split(size);
			freeBlocks.insert(i, split.second);
			return split.first;
		} else {
			// b.size == size then return b
			return b;
		}
	}
}

// coalesce the given block with other free blocks
//
// @param[in] freeBlocks the address-ordered set of free blocks.
// @param b the block to coalesce.
// @return the coalesced block.
inline Block coalesceBlock(std::set<Block>& freeBlocks, Block const&b){
	// return the given block in case is not valid
	if(!b.isValid()) return b;

	// find the right place (in ascending address order) to insert the block
	auto const next = freeBlocks.lower_bound(b);
	auto const previous = next == freeBlocks.cend() ? next : std::prev(next);

	// coalesce with neighboring blocks
	bool const mergePrev = previous->isContiguousBefore(b);
	bool const mergeNext = next != freeBlocks.cend() && b.isContiguousBefore(*next);

	Block merged{};
	if(mergePrev && mergeNext) {
		// if can merge with prev and next neighbors
		merged = previous->merge(b).merge(*next);

		freeBlocks.erase(previous);

		auto const i = freeBlocks.erase(next);
		freeBlocks.insert(i, merged);
	} else if(mergePrev) {
		// if only can merge with prev neighbor
		merged = previous->merge(b);

		auto const i = freeBlocks.erase(next);
		freeBlocks.insert(i, merged);
	} else if(mergeNext) {
		// if only can merge with next neighbor
		merged = b.merge(*next);

		auto const i = freeBlocks.erase(next);
		freeBlocks.insert(i, merged);
	} else {
		// if can't be merge with either
		freeBlocks.emplace(b);
		merged = b;
	}

	return merged;
}

} // namaspace Arena
} // namaspace detail
} // namaspace MemoryResource
} // namaspace Nmm
} // namaspace Acts