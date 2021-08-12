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
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/detail/GlobalArena.hpp"

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

// An arena for allocating memory for a thread
// An arena is a per-thread or per-non-default-stream memory pool. It allocates
// superblocks from the global arena, and return them when the superblocks become empty
//
// @tparam Upstream Memory resource to use for allocating the global arena. Implements
// MemoryResource::DeviceMemoryResource interface
template <typename Upstream>
class Arena {
	public:
		// Construct an `Arena`
		//
		// @param[in] GlobalArena the global arena from withc to allocate superblocks
		explicit Arena(GlobalArena<Upstream>& globalArena) : globalArena_{globalArena} {}

		// Allocates memory of size at least `bytes`
		//
		// @param[in] bytes the size in bytes of the allocation
		// @return void* pointer to the newly allocated memory
		void* allocate(std::size_t bytes) {
			lockGuard lock(mtx_);

			auto const b = getBlock(bytes);
			allocatedBlocks_.emplace(b.pointer(), b);

			return b.pointer();
		}

		// Deallocate memory pointed to by `p`, and possibly return superblocks to upstream.
		// return the block to the set that have the free blocks
		//
		// @param[in] p the pointer of the memory
		// @param[in] bytes the size in bytes of the deallocation
		// @return if the allocation was found, false otherwise
		bool deallocate(void* p, std::size_t bytes, CudaStreamView stream) {
			lockGuard lock(mtx_);

			auto const b = freeBlock(p, bytes);
			if(b.isValid()) {
				auto const merged = coalesceBlock(freeBlocks_, b);
				shrinkArena(merged, stream);
			}

			return b.isValid();
		}

		// Deallocate memory pointed to by `p`, and keeping all free superblocks.
		// return the block to the set that have the free blocks. This method
		// is used when deallocating from another arena, since we don't have access the
		// respective CudaStreamView of the areana.
		//
		// @param[in] p the pointer of the memory
		// @param[in] bytes the size in bytes of the deallocation
		// @return if the allocation was found, false otherwise
		bool deallocate(void* p, std::size_t bytes) {
			lockGuard lock(mtx_);

			auto const b = freeBlock(p, bytes);
			if(b.isValid()) {
				globalArena_.deallocate(b);
			}

			return b.isValid();
		}

		// Clean the arena and deallocate free blocks from the global arena.
	    //
	    // This is only needed when a per-thread arena is about to die.
		void clean() {
			lockGuard lock(mtx_);
			globalArena_.deallocate(freeBlocks_);
			freeBlocks_.clear();
			allocatedBlocks_.clear();
		}

	private:
		using lockGuard = std::lock_guard<std::mutex>;

		// @brief Get an available memory block of at least `size` bytes.
		//
		// @param[in] size The number of bytes to allocate.
		// @return block A block of memory of at least `size` bytes.
		Block getBlock(std::size_t size) {
			if(size < minimumSuperblockSize) {
				auto const b = firstFit(freeBlocks_, size);
				if(b.isValid()) {
					return b;
				}
			}

			auto const superblock = expandArena(size);
			coalesceBlock(freeBlocks_, superblock);
			return firstFit(freeBlocks_, size);
		}

		// Allocate space from upstream to supply the arena and return a superblock.
		// 
		// @return block A superblock.
		Block expandArena(std::size_t size) {
			auto const superblockSize = std::max(size, minimumSuperblockSize);
			return globalArena_.allocate(superblockSize);
		}

		// Finds, frees and returns the block associated with pointer `p`.
		//
		// @param[in] p The pointer to the memory to free.
		// @param[in] size The size of the memory to free. Must be equal to the original allocation size.
		// return The (now freed) block associated with `p`. The caller is expected to return the block
		// to the arena.
		Block freeBlock(void* p, std::size_t size) noexcept {
			auto const i = allocatedBlocks_.find(p);

			if(i == allocatedBlocks_.end()) {
				return {};
			}

			auto const found = i->second;
			//assert if found.size == size
			allocatedBlocks_.erase(i);

			return found;
		}

		// Shrink this arena by returning free superblocks to upstream.
		//
		// @param b The block that can be used to shrink the arena.
		// @param stream Stream on which to perform shrinking.
		void shrinkArena(Block const& b, CudaStreamView stream) {
			if(!b.isSuperblock()) return;

			stream.synchronize();

			globalArena_.deallocate(b);
			freeBlocks_.erase(b);
		}

		GlobalArena<Upstream>& globalArena_;
		std::set<Block> freeBlocks_;
		std::unordered_map<void*, Block> allocatedBlocks_;
		mutable std::mutex mtx_;
};// class Arena


} // namaspace Arena
} // namaspace detail
} // namaspace MemoryResource
} // namaspace Nmm
} // namaspace Acts