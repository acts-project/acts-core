// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Seeding2/Details/CountDublets.hpp"
#include "Acts/Plugins/Cuda/Seeding2/Details/Types.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager.hpp"

#include "../Utilities/ErrorCheck.cuh"

// CUDA include(s).
#include <cuda_runtime.h>

// System include(s).
#include <algorithm>

#include <iostream>

namespace Acts {
namespace Cuda {
namespace Kernels {

/// Kernel performing the dublet countint
///
/// This is pretty much just a copy of "reduce3" kernel from the samples shipped
/// with CUDA. It finds the sums and maxima of a couple of parameters across all
/// the dublets found by @c Acts::Cuda::Details::findDublets.
///
/// @param[in] nMiddleSPs The number of middle spacepoints
/// @param[in] middleBottomCounts 1-D array of the number of middle-bottom
///            dublets found for each middle spacepoint
/// @param[in] middleTopCounts 1-D array of the number of middle-top dublets
///            found for each middle spacepoint
/// @param[out] dubletCounts 1-D array of @c Acts::Cuda::Details::DubletCounts
///             objects for each execution block
///
__global__ void countDublets(std::size_t nMiddleSPs,
                             const unsigned int* middleBottomCounts,
                             const unsigned int* middleTopCounts,
                             Details::DubletCounts* dubletCounts) {
  // Sum shared by all threads in a single execution block.
  extern __shared__ Details::DubletCounts sum[];

  // Get the thread identifier. Note that the kernel launch requests half as
  // many threads than how many elements we have in the arrays.
  const int middleIndex = blockIdx.x * blockDim.x * 2 + threadIdx.x;

  Details::DubletCounts thisSum;
  if (middleIndex < nMiddleSPs) {
    thisSum.nDublets =
        (middleBottomCounts[middleIndex] + middleTopCounts[middleIndex]);
    thisSum.nTriplets =
        (middleBottomCounts[middleIndex] * middleTopCounts[middleIndex]);
    thisSum.maxMBDublets = middleBottomCounts[middleIndex];
    thisSum.maxMTDublets = middleTopCounts[middleIndex];
    thisSum.maxTriplets = thisSum.nTriplets;
  }
  if (middleIndex + blockDim.x < nMiddleSPs) {
    thisSum.nDublets += (middleBottomCounts[middleIndex + blockDim.x] +
                         middleTopCounts[middleIndex + blockDim.x]);
    thisSum.nTriplets += (middleBottomCounts[middleIndex + blockDim.x] *
                          middleTopCounts[middleIndex + blockDim.x]);
    thisSum.maxMBDublets =
        max(middleBottomCounts[middleIndex + blockDim.x], thisSum.maxMBDublets);
    thisSum.maxMTDublets =
        max(middleTopCounts[middleIndex + blockDim.x], thisSum.maxMTDublets);
    thisSum.maxTriplets = max((middleBottomCounts[middleIndex + blockDim.x] *
                               middleTopCounts[middleIndex + blockDim.x]),
                              thisSum.maxTriplets);
  }

  // Load the first sum step into shared memory.
  sum[threadIdx.x] = thisSum;
  __syncthreads();

  // Do the summation in some iterations.
  for (unsigned int i = blockDim.x / 2; i > 0; i >>= 1) {
    if (threadIdx.x < i) {
      const Details::DubletCounts& otherSum = sum[threadIdx.x + i];
      thisSum.nDublets += otherSum.nDublets;
      thisSum.nTriplets += otherSum.nTriplets;
      thisSum.maxMBDublets = max(thisSum.maxMBDublets, otherSum.maxMBDublets);
      thisSum.maxMTDublets = max(thisSum.maxMTDublets, otherSum.maxMTDublets);
      thisSum.maxTriplets = max(thisSum.maxTriplets, otherSum.maxTriplets);
      sum[threadIdx.x] = thisSum;
    }
    __syncthreads();
  }

  // Write the result of this execution block into the global memory.
  if (threadIdx.x == 0) {
    dubletCounts[blockIdx.x] = thisSum;
  }
  return;
}

}  // namespace Kernels

namespace Details {

DubletCounts countDublets(
    std::size_t maxBlockSize, std::size_t nMiddleSP,
    const device_array<unsigned int>& middleBottomCountArray,
    const device_array<unsigned int>& middleTopCountArray) {
  // Calculate the parallelisation for the dublet counting.
  const int numBlocks = (nMiddleSP + maxBlockSize - 1) / maxBlockSize;
  const int sharedMem = maxBlockSize * sizeof(DubletCounts);
  
  // Create the small memory block in which we will get the count back for each
  // execution block.
  auto dubletCountsDevice = make_device_array<DubletCounts>(numBlocks);

  // Run the reduction kernel.
  Kernels::countDublets<<<numBlocks, maxBlockSize, sharedMem>>>(
      nMiddleSP, middleBottomCountArray.get(), middleTopCountArray.get(),
      dubletCountsDevice.get());
  ACTS_CUDA_ERROR_CHECK(cudaGetLastError());
  ACTS_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

  // Copy the sum(s) back to the host.
  auto dubletCountsHost = make_host_array<DubletCounts>(numBlocks);
  copyToHost(dubletCountsHost, dubletCountsDevice, numBlocks);

  // Perform the final summation on the host. Assuming that the number of
  // middle space points is not so large that it would make sense to do the
  // summation iteratively on the device. (We should get one result object per
  // 1024 middle spacepoints on any modern GPU.)
  DubletCounts result;
  for (int i = 0; i < numBlocks; ++i) {
    result.nDublets += dubletCountsHost.get()[i].nDublets;
    result.nTriplets += dubletCountsHost.get()[i].nTriplets;
    result.maxMBDublets =
        std::max(dubletCountsHost.get()[i].maxMBDublets, result.maxMBDublets);
    result.maxMTDublets =
        std::max(dubletCountsHost.get()[i].maxMTDublets, result.maxMTDublets);
    result.maxTriplets =
        std::max(dubletCountsHost.get()[i].maxTriplets, result.maxTriplets);
  }

  /*
  Cuda::Nmm::MemoryResource::HostMemoryResource *hmr = Cuda::getCPUmmr();
  hmr->deallocate(static_cast<void*>(dubletCountsHost.get()), numBlocks);
  */

  return result;
}

}  // namespace Details
}  // namespace Cuda
}  // namespace Acts
