// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Seeding2/Details/FindDublets.hpp"
#include "Acts/Plugins/Cuda/Seeding2/Details/Types.hpp"

#include "../Utilities/ErrorCheck.cuh"
#include "../Utilities/MatrixMacros.hpp"

// System include(s).
#include <cassert>
#include <cmath>

#include <iostream>

namespace {

/// Type of "other spacepoint" passed to the kernel
enum OtherSPType : int {
  BottomSP = 0,  //< The "other" spacepoint is a bottom one
  TopSP = 1      //< The "other" spacepoint is a top one
};

}  // namespace

namespace Acts {
namespace Cuda {
namespace Kernels {

/// @name Helper functions for calculating "deltaR" between spacepoints
/// @{

/// Prototype for the @c Acts::Cuda::Kernels::getDeltaR functions
template <int SPType>
__device__ float getDeltaR(float /*middleR*/, float /*otherR*/) {
  // This function should *never* be called.
  assert(false);
  return 0.0f;
}

/// Function specialisation for middle-bottom spacepoint pairs
template <>
__device__ float getDeltaR<BottomSP>(float middleR, float bottomR) {
  return middleR - bottomR;
}

/// Function specialisation for middle-top spacepoint pairs
template <>
__device__ float getDeltaR<TopSP>(float middleR, float topR) {
  return topR - middleR;
}

/// @}

/// @name Helper functions for calculating "cotTheta" between spacepoints
/// @{

/// Prototype for the @c Acts::Cuda::Kernels::getCotTheta functions
template <int SPType>
__device__ float getCotTheta(float /*middleZ*/, float /*otherZ*/,
                             float /*deltaR*/) {
  // This function should *never* be called.
  assert(false);
  return 0.0f;
}

/// Function specialisation for middle-bottom spacepoint pairs
template <>
__device__ float getCotTheta<BottomSP>(float middleZ, float bottomZ,
                                       float deltaR) {
  return (middleZ - bottomZ) / deltaR;
}

/// Function specialisation for middle-top spacepoint pairs
template <>
__device__ float getCotTheta<TopSP>(float middleZ, float topZ, float deltaR) {
  return (topZ - middleZ) / deltaR;
}

/// @}

/// Kernel performing the spacepoint dublet finding
///
/// @tparam SPType The "other" spacepoint type used in the call
///         (@c ::BottomSP or @c TopSP)
/// @param[in] nMiddleSPs The number of middle spacepoints in @c middleSPs
/// @param[in] middleSPs Properties of all of the middle spacepoints
/// @param[in] nOtherSPs The number of other (bottom or top) spacepoints in
///            @c otherSPs
/// @param[in] otherSPs Properties of all of the other (bottom or top)
///            spacepoints
/// @param[in] deltaRMin Configuration parameter from @c Acts::SeedfinderConfig
/// @param[in] deltaRMax Configuration parameter from @c Acts::SeedfinderConfig
/// @param[in] cotThetaMax Configuration parameter from
///            @c Acts::SeedfinderConfig
/// @param[in] collisionRegionMin Configuration parameter from
///            @c Acts::SeedfinderConfig
/// @param[in] collisionRegionMax Configuration parameter from
///            @c Acts::SeedfinderConfig
/// @param[out] dubletCounts 1-D array of the middle-other dublets found
///             for each middle spacepoint
/// @param[out] dublets 2-D matrix of size @c nMiddleSPs x @c nOtherSPs, holding
///             the "other" spacepoint indices for the identified dublets
///
template <int SPType>
__global__ void findDublets(std::size_t nMiddleSPs,
                            const Details::SpacePoint* middleSPs,
                            std::size_t nOtherSPs,
                            const Details::SpacePoint* otherSPs,
                            float deltaRMin, float deltaRMax, float cotThetaMax,
                            float collisionRegionMin, float collisionRegionMax,
                            unsigned int* dubletCounts, std::size_t* dublets) {
  // Figure out which dublet the kernel operates on.
  const std::size_t middleIndex = blockIdx.x * blockDim.x + threadIdx.x;
  const std::size_t otherIndex = blockIdx.y * blockDim.y + threadIdx.y;

  // If we're outside of bounds, stop here.
  if ((middleIndex >= nMiddleSPs) || (otherIndex >= nOtherSPs)) {
    return;
  }

  // Calculate variables used in the compatibility check.
  const float deltaR = getDeltaR<SPType>(middleSPs[middleIndex].radius,
                                         otherSPs[otherIndex].radius);
  const float cotTheta = getCotTheta<SPType>(middleSPs[middleIndex].z,
                                             otherSPs[otherIndex].z, deltaR);
  const float zOrigin =
      middleSPs[middleIndex].z - middleSPs[middleIndex].radius * cotTheta;

  // Perform the compatibility check.
  const bool isCompatible =
      ((deltaR >= deltaRMin) && (deltaR <= deltaRMax) &&
       (fabs(cotTheta) <= cotThetaMax) && (zOrigin >= collisionRegionMin) &&
       (zOrigin <= collisionRegionMax));

  // If they are compatible, save their indices into the output matrix.
  if (isCompatible) {
    const unsigned int dubletRow = atomicAdd(dubletCounts + middleIndex, 1);
    ACTS_CUDA_MATRIX2D_ELEMENT(dublets, nMiddleSPs, nOtherSPs, middleIndex,
                               dubletRow) = otherIndex;
  }
  return;
}

}  // namespace Kernels

namespace Details {

void findDublets(std::size_t maxBlockSize, std::size_t nBottomSPs,
                 const device_array<SpacePoint>& bottomSPs,
                 std::size_t nMiddleSPs,
                 const device_array<SpacePoint>& middleSPs, std::size_t nTopSPs,
                 const device_array<SpacePoint>& topSPs, float deltaRMin,
                 float deltaRMax, float cotThetaMax, float collisionRegionMin,
                 float collisionRegionMax,
                 device_array<unsigned int>& middleBottomCounts,
                 device_array<std::size_t>& middleBottomDublets,
                 device_array<unsigned int>& middleTopCounts,
                 device_array<std::size_t>& middleTopDublets) {
  // Calculate the parallelisation for the middle<->bottom spacepoint
  // compatibility flagging.
  const dim3 blockSizeMB(1, maxBlockSize);
  const dim3 numBlocksMB((nMiddleSPs + blockSizeMB.x - 1) / blockSizeMB.x,
                         (nBottomSPs + blockSizeMB.y - 1) / blockSizeMB.y);

  // Launch the middle-bottom dublet finding.
  Kernels::findDublets<BottomSP><<<numBlocksMB, blockSizeMB>>>(
      nMiddleSPs, middleSPs.get(), nBottomSPs, bottomSPs.get(), deltaRMin,
      deltaRMax, cotThetaMax, collisionRegionMin, collisionRegionMax,
      middleBottomCounts.get(), middleBottomDublets.get());
  ACTS_CUDA_ERROR_CHECK(cudaGetLastError());

  // Calculate the parallelisation for the middle<->top spacepoint
  // compatibility flagging.
  const dim3 blockSizeMT(1, maxBlockSize);
  const dim3 numBlocksMT((nMiddleSPs + blockSizeMT.x - 1) / blockSizeMT.x,
                         (nTopSPs + blockSizeMT.y - 1) / blockSizeMT.y);

  // Launch the middle-bottom dublet finding.
  Kernels::findDublets<TopSP><<<numBlocksMT, blockSizeMT>>>(
      nMiddleSPs, middleSPs.get(), nTopSPs, topSPs.get(), deltaRMin, deltaRMax,
      cotThetaMax, collisionRegionMin, collisionRegionMax,
      middleTopCounts.get(), middleTopDublets.get());
  ACTS_CUDA_ERROR_CHECK(cudaGetLastError());
  ACTS_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
  return;
}

}  // namespace Details
}  // namespace Cuda
}  // namespace Acts
