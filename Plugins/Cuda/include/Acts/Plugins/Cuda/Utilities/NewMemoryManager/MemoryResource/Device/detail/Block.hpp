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

// Minimum size of a Superblock (256 KiB)
constexpr std::size_t minimumSuperblockSize = 1u << 18u;

class Block {

	public:
		// construct a default Block.
		Block() = default;

		// construct a Block given a pointer and size.
		// 
		// @param[in] pointer the address for the beginning of the Block.
		// @param[in] size the size of the Block
		Block(char* pointer, std::size_t size) : pointer_(pointer), size_(size) {}

		// construct a Block given a pointer and size.
		// 
		// @param[in] pointer the address for the beginning of the Block.
		// @param[in] size the size of the Block
		Block(void* pointer, std::size_t size) : pointer_(static_cast<char*>(pointer)), size_(size){}

		// returns the underlying pointer
		void* pointer() const { return pointer_; }

		// returns the size of the Block
		std::size_t size() const { return size_; }

		// returns true if this Block is valid (non-null), false otherwise
		bool isValid() const { return pointer_ != nullptr; }

		// returns true if this Block is a Superblock, false otherwise
		bool isSuperblock() const { return size_ >= minimumSuperblockSize; }

		// verifies wheter this Block can be merged to the beginning of Block b
		//
		// @param[in] b the Block to check for contiguity
		// @return true if this Block's `pointer` + `size` == `b.ptr`, and `not b.isHead`,
		// false otherwise 
		bool isContiguousBefore(Block const& b) const { return pointer_ + size_ == b.pointer_; }

		// is this Block large enough to fit that size of bytes?
		//
		// @param[in] sizeOfBytes the size in bytes to check for fit
		// @return true if this Block is at least sizeOfBytes
		bool fits(std::size_t sizeOfBytes) const { return size_ >= sizeOfBytes; }

		// split this Block into two by the given size
		//
		// @param[in] size the size in bytes of the first Block
		// @return std::pair<Block, Block> a pair of Blocks split by size
		std::pair<Block, Block> split(std::size_t size) const {
			//assert condition of size_ >= size
			if(size_ > size) {
				return {{pointer_, size}, {pointer_ + size, size_ - size}};
			} else {
				return {*this, {}};
			}
		}

		// coalesce two contiguos Blocks into one, this->isContiguousBefore(b) 
		// must be true
		// 
		// @param[in] b Block to merge
		// @return Block the merged Block
		Block merge(Block const& b) const {
			//assert condition isContiguousBefore(b)
			return {pointer_, size_ + b.size_};
		}

		// used by std::set to compare Blocks
		bool operator<(Block const& b) const { return pointer_ < b.pointer_; }

	private:
		char* pointer_;      // raw memory pointer
		std::size_t size_{}; // size in bytes
}; // class Block

} // namaspace Arena
} // namaspace detail
} // namaspace MemoryResource
} // namaspace Nmm
} // namaspace Acts