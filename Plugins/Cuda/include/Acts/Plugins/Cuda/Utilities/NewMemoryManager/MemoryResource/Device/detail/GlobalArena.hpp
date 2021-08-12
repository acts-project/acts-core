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
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/detail/Methods.hpp"

#include <vector>

namespace Acts {
namespace Nmm {
namespace MemoryResource {
namespace detail {
namespace Arena {

template <typename Upstream>
class GlobalArena final {
	public:
		// default initial size for the global arena
		static constexpr std::size_t defaultInitialSize = std::numeric_limits<std::size_t>::max();
		// default maximum size for the global arena
		static constexpr std::size_t defaultMaximumSize = std::numeric_limits<std::size_t>::max();
		// reserved memory that should not be allocated (64 MiB)
		static constexpr std::size_t reserverdSize = 1u << 26u;
		
		// Disable copy (and move) semantics
		GlobalArena(const GlobalArena&) = delete;
		GlobalArena& operator=(const GlobalArena&) = delete;

		// construct a global arena
		//
		// @param[in] upstreamMemoryResource the memory resource from which to allocate 
		// blocks for the pool
		// @param[in] initialSize minimum size, in bytes, of the initial global arena. 
		// Defaults to all the available memory on the current device.
		// @param[in] maximumSize maximum size, in bytes, that the global arena can grow to. 
		// Defaults to all of the available memory on the current device
		GlobalArena(Upstream* upstreamMemoryResource, std::size_t initialSize, std::size_t maximumSize)
			: upstreamMemoryResource_{upstreamMemoryResource}, maximumSize_{maximumSize} {
				// assert unexpected null upstream pointer
				// assert initial arena size required to be a multiple of 256 bytes
				// assert maximum arena size required to be a multiple of 256 bytes

				if(initialSize == defaultInitialSize || maximumSize == defaultMaximumSize) {
					std::size_t free{}, total{};
					cudaMemGetInfo(&free, &total);
					if(initialSize == defaultInitialSize) {
						initialSize = alignUp(std::min(free, total / 2));
					}
					if(maximumSize == defaultMaximumSize) {
						maximumSize_ = alignDown(free) - reserverdSize;
					}
				}
				// initial size exceeds the maxium pool size
				freeBlocks_.emplace(expandArena(initialSize));
		}

		// Destroy the gloabl arean and deallocate all memory using the upstream resource
		~GlobalArena() {
			lockGuard lock(mtx_);
			for(auto const& b : upstreamBlocks_) {
				upstreamMemoryResource_->deallocate(b.pointer(), b.size());
			}
		}

		// Allocates memory of size at least `sizeOfbytes`
		//
		// @param[in] sizeOfBytes the size in bytes of the allocation
		// @retyrb block pointer to the newly allocated memory
		Block allocate(std::size_t sizeOfbytes) {
			lockGuard lock(mtx_);
			return getBlock(sizeOfbytes);
		}

		// Deallocate memory pointer to by the block b
		// 
		// @param[in] b pointer of block to be deallocated
		void deallocate(Block const& b) {
			lockGuard lock(mtx_);
			coalesceBlock(freeBlocks_, b);
		}

		// Deallocate memory of a set of blocks
		//
		// @param[in] freeBlocks set of block to be free
		void deallocate(std::set<Block> const& freeBlocks) {
			lockGuard lock(mtx_);
			for(auto const& b : freeBlocks) {
				coalesceBlock(freeBlocks_, b);
			}
		}

	private:
		using lockGuard = std::lock_guard<std::mutex>;

		// Get an available memory block of at least `size` bytes
		//
		// @param[in] size the number of bytes to allocate
		// @return a block of memory of at least `size` bytes
		Block getBlock(std::size_t size) {
			auto const b = firstFit(freeBlocks_, size);
			if(b.isValid()) return b;

			auto const upstreamBlock = expandArena(sizeToGrow(size));
			coalesceBlock(freeBlocks_, upstreamBlock);
			return firstFit(freeBlocks_, size);
		}

		// Get the size to grow the global arena given the requested `size` bytes
		//
		// @param[in] size the number of bytes required
		// @return the size for the arena to grow
		constexpr std::size_t sizeToGrow(std::size_t size) const {
			/*
			case if maximum pool size exceeded
			if(currentSize_ + size > maximumSize_) {
		
			}
			*/
			return maximumSize_ - currentSize_;
		}

		// Allocate space from upstream to supply the arena and return a sufficiently 
		// sized block.
		//
		// @param[in] size the minimum size to allocate
		// @return a bock of at least `size` bytes
		Block expandArena(std::size_t size) {
			upstreamBlocks_.push_back({upstreamMemoryResource_->allocate(size), size});
			currentSize_ += size;
			return upstreamBlocks_.back();
		}

		// The upstream resource to allocate memory from
		Upstream* upstreamMemoryResource_;
		// The maximum size of the global arena
		std::size_t maximumSize_;
		// The current size of the global arena
		std::size_t currentSize_{};
		// Address-ordered set of free blocks
		std::set<Block> freeBlocks_;
		// Blocks allocated from upstream so that they can be quickly freed
		std::vector<Block> upstreamBlocks_;
		// Mutex for exclusive lock
		mutable std::mutex mtx_;
};// class GlobalArena

} // namaspace Arena
} // namaspace detail
} // namaspace MemoryResource
} // namaspace Nmm
} // namaspace Acts