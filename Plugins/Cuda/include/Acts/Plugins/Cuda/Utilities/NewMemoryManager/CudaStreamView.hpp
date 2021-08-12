// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once	

// CUDA plugin include(s).
//#include "../../../../../../src/Utilities/ErrorCheck.cuh"

// CUDA include(s).
#include "cuda.h"
#include "cuda_runtime.h"

// System include(s).
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <iostream>

namespace Acts {
namespace Nmm {

class CudaStreamView {
	public:
		constexpr CudaStreamView()                        = default;
		constexpr CudaStreamView(CudaStreamView const&) = default;
		constexpr CudaStreamView(CudaStreamView&&)      = default;
		constexpr CudaStreamView& operator=(CudaStreamView const&) = default;
		constexpr CudaStreamView& operator=(CudaStreamView&&) = default;
		~CudaStreamView()                                       = default;

		// Implicit conversion from cudaStream_t
		constexpr CudaStreamView(cudaStream_t stream) : stream_{stream} {}

		// Returns the wrppped stream
		constexpr cudaStream_t value() const noexcept {
			return stream_;
		}

		// Explicit conversion to cudaStream_t
		explicit constexpr operator cudaStream_t() const noexcept {
			return value();
		}

		// Return true if the wrapped stream is the CUDA per-thread default stream
		bool is_per_thread_default() const noexcept {
			#ifdef CUDA_API_PER_THREAD_DEFAULT_STREAM
				return value() == cudaStreamPerThread || value() == 0;
			#else
				return value() == cudaStreamPerThread;
			#endif
		}

		// Return true if the wrapped stream is explicitly the CUDA legacy default stream
		bool is_default() const noexcept {
			#ifdef CUDA_API_PER_THREAD_DEFAULT_STREAM
				return value() == cudaStreamLegacy;
			#else
				return value() == cudaStreamLegacy || value() == 0;
			#endif
		}

		// Synchronize the viewed CUDA stream
		void synchronize() const {
			cudaStreamSynchronize(stream_);
		}

		// Synchronize the viewed CUDA stream, don't throw if there is an error
		/*
		void synchronize_no_throw() const noexcept {
			ACTS_CUDA_ERROR_CHECK(cudaStreamSynchronize(stream_));
		}
		*/
	private:
		cudaStream_t stream_{0};
}; // class CudaStreamView

// Static CudaStreamView of the default stream (stream 0), for convenience
static constexpr CudaStreamView cuda_stream_default{};

// Static CudaStreamView of cudaStreamLegacy, for convenience
static CudaStreamView cuda_stream_legacy{cudaStreamLegacy};

// Static CudaStreamView of cudaStreamPerThread, for convenience
static CudaStreamView cuda_stream_per_thread{cudaStreamPerThread};

// Equality ciomparison operator for streams
// 
// @param[in] lhs the first stream view to compare
// @param[in] rhs the second stream view to compare
// @return true if equal, false if unequal
inline bool operator==(CudaStreamView lhs, CudaStreamView rhs) {
	return lhs.value() == rhs.value();
}

// Inequality comparison operator for streams
//
// @param[in] lhs the first stream view to compare
// @param[in] rhs the second stream view to compare
// @return true if unequal, false if equal
inline bool operator!=(CudaStreamView lhs, CudaStreamView rhs) { 
	return !(lhs == rhs); 
}

// Output stream operator for printing / logging streams
//
// @param[in] os the output ostream
// @param[in] sv the CudaStreamView to output
// @return std::ostream& the output ostream
inline std::ostream& operator<<(std::ostream& os, CudaStreamView sv) {
	os << sv.value();
	return os;
}

} // namespace Nmm
} // namespace Acts