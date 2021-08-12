// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Seeding2/Details/FindTriplets.hpp"
#include "Acts/Plugins/Cuda/Seeding2/Details/Types.hpp"
#include "Acts/Plugins/Cuda/Seeding2/TripletFilterConfig.hpp"
#include "Acts/Plugins/Cuda/Utilities/MemoryManager.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager.hpp"

#include "../Utilities/ErrorCheck.cuh"
#include "../Utilities/MatrixMacros.hpp"

// Acts include(s).
#include "Acts/Seeding/SeedFilterConfig.hpp"

// CUDA include(s).
#include <cuda_runtime.h>

// System include(s).
#include <cassert>
#include <cmath>
#include <cstring>

#include <iostream>

namespace Acts {
namespace Cuda {
namespace Kernels {

/// Function performing coordinate transformation for one spacepoint pair
///
/// @param spM    The middle spacepoint to use
/// @param sp     The "other" spacepoint to use
/// @param bottom @c true If the "other" spacepoint is a bottom one, @c false
///               otherwise
__device__ Details::LinCircle transformCoordinates(
    const Details::SpacePoint& spM, const Details::SpacePoint& sp,
    bool bottom) {
  // Create the result object.
  Details::LinCircle result;

  // Parameters of the middle spacepoint.
  const float cosPhiM = spM.x / spM.radius;
  const float sinPhiM = spM.y / spM.radius;

  // (Relative) Parameters of the spacepoint being transformed.
  const float deltaX = sp.x - spM.x;
  const float deltaY = sp.y - spM.y;
  const float deltaZ = sp.z - spM.z;

  // calculate projection fraction of spM->sp vector pointing in same
  // direction as
  // vector origin->spM (x) and projection fraction of spM->sp vector pointing
  // orthogonal to origin->spM (y)
  const float x = deltaX * cosPhiM + deltaY * sinPhiM;
  const float y = deltaY * cosPhiM - deltaX * sinPhiM;
  // 1/(length of M -> SP)
  const float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
  const float iDeltaR = sqrtf(iDeltaR2);
  //
  const int bottomFactor = 1 * (int(!bottom)) - 1 * (int(bottom));
  // cot_theta = (deltaZ/deltaR)
  const float cot_theta = deltaZ * iDeltaR * bottomFactor;
  // VERY frequent (SP^3) access
  result.cotTheta = cot_theta;
  // location on z-axis of this SP-duplet
  result.Zo = spM.z - spM.radius * cot_theta;
  result.iDeltaR = iDeltaR;
  // transformation of circle equation (x,y) into linear equation (u,v)
  // x^2 + y^2 - 2x_0*x - 2y_0*y = 0
  // is transformed into
  // 1 - 2x_0*u - 2y_0*v = 0
  // using the following m_U and m_V
  // (u = A + B*v); A and B are created later on
  result.U = x * iDeltaR2;
  result.V = y * iDeltaR2;
  // error term for sp-pair without correlation of middle space point
  result.Er = ((spM.varianceZ + sp.varianceZ) +
               (cot_theta * cot_theta) * (spM.varianceR + sp.varianceR)) *
              iDeltaR2;
  return result;
}

/// Kernel performing coordinate transformation on all created dublets
///
/// @param[in] nDublets The total number of dublets found
/// @param[in] maxMBDublets The maximal number of middle-bottom dublets found
///            for any middle spacepoint
/// @param[in] maxMTDublets The maximal number of middle-top dublets found for
///            any middle spacepoint
/// @param[in] nBottomSPs The number of bottom spacepoints in @c bottomSPs
/// @param[in] bottomSPs Properties of all of the bottom spacepoints
/// @param[in] nMiddleSPs The number of middle spacepoints in @c middleSPs
/// @param[in] middleSPs Properties of all of the middle spacepoints
/// @param[in] nTopSPs The number of top spacepoints in @c topSPs
/// @param[in] topSPs Properties of all of the top spacepoints
/// @param[in] middleBottomCounts 1-D array of the number of middle-bottom
///            dublets found for each middle spacepoint
/// @param[in] middleBottomDublets 2-D matrix of size
///            @c nMiddleSPs x @c nBottomSPs, holding the bottom spacepoint
///            indices for the identified middle-bottom dublets
/// @param[in] middleTopCounts 1-D array of the number of middle-top dublets
///            found for each middle spacepoint
/// @param[in] middleTopDublets 2-D matrix of size
///            @c nMiddleSPs x @c nTopSPs, holding the top spacepoint
///            indices for the identified middle-top dublets
/// @param[out] bottomSPLinTransArray 2-dimensional matrix indexed the same way
///             as @c middleBottomDublets
/// @param[out] topSPLinTransArray 2-dimensional matrix indexed the same way as
///             @c middleTopDublets
///
__global__ void transformCoordinates(
    unsigned int nDublets, unsigned int maxMBDublets, unsigned int maxMTDublets,
    std::size_t nBottomSPs, const Details::SpacePoint* bottomSPs,
    std::size_t nMiddleSPs, const Details::SpacePoint* middleSPs,
    std::size_t nTopSPs, const Details::SpacePoint* topSPs,
    const unsigned int* middleBottomCounts,
    const std::size_t* middleBottomDublets, const unsigned int* middleTopCounts,
    const std::size_t* middleTopDublets,
    Details::LinCircle* bottomSPLinTransArray,
    Details::LinCircle* topSPLinTransArray) {
  // Get the global index.
  const int dubletIndex = blockIdx.x * blockDim.x + threadIdx.x;

  // If we're out of bounds, finish right away.
  if (dubletIndex >= nDublets) {
    return;
  }

  // Find the dublet to transform.
  std::size_t middleIndex = 0;
  int runningIndex = dubletIndex;
  int tmpValue = 0;
  while (runningIndex >= (tmpValue = (middleBottomCounts[middleIndex] +
                                      middleTopCounts[middleIndex]))) {
    middleIndex += 1;
    assert(middleIndex < nMiddleSPs);
    runningIndex -= tmpValue;
  }
  const bool transformBottom =
      ((runningIndex < middleBottomCounts[middleIndex]) ? true : false);
  const std::size_t bottomMatrixIndex = (transformBottom ? runningIndex : 0);
  const std::size_t topMatrixIndex =
      (transformBottom ? 0 : runningIndex - middleBottomCounts[middleIndex]);

  // Perform the transformation.
  if (transformBottom) {
    const std::size_t bottomIndex =
        ACTS_CUDA_MATRIX2D_ELEMENT(middleBottomDublets, nMiddleSPs, nBottomSPs,
                                   middleIndex, bottomMatrixIndex);
    assert(bottomIndex < nBottomSPs);
    ACTS_CUDA_MATRIX2D_ELEMENT(bottomSPLinTransArray, nMiddleSPs, maxMBDublets,
                               middleIndex, bottomMatrixIndex) =
        transformCoordinates(middleSPs[middleIndex], bottomSPs[bottomIndex],
                             true);
  } else {
    const std::size_t topIndex = ACTS_CUDA_MATRIX2D_ELEMENT(
        middleTopDublets, nMiddleSPs, nTopSPs, middleIndex, topMatrixIndex);
    assert(topIndex < nTopSPs);
    ACTS_CUDA_MATRIX2D_ELEMENT(topSPLinTransArray, nMiddleSPs, maxMTDublets,
                               middleIndex, topMatrixIndex) =
        transformCoordinates(middleSPs[middleIndex], topSPs[topIndex], false);
  }

  return;
}

/// Kernel used for finding all the triplet candidates
///
/// @param[in] middleIndexStart The middle spacepoint index that the kernel was
///            "started from"
/// @param[in] maxMBDublets The maximal number of middle-bottom dublets found
///            for any middle spacepoint
/// @param[in] maxMTDublets The maximal number of middle-top dublets found for
///            any middle spacepoint
/// @param[in] maxTriplets The maximum number of triplets for which memory is
///            booked
/// @param[in] nParallelMiddleSPs The number of middle spacepoints that the
///            "largest" kernels may be started on in parallel
/// @param[in] nMiddleSPsProcessed The number of middle spacepoints that the
///            kernel was started on in parallel
/// @param[in] nBottomSPs The number of bottom spacepoints in @c bottomSPs
/// @param[in] bottomSPs Properties of all of the bottom spacepoints
/// @param[in] nMiddleSPs The number of middle spacepoints in @c middleSPs
/// @param[in] middleSPs Properties of all of the middle spacepoints
/// @param[in] nTopSPs The number of top spacepoints in @c topSPs
/// @param[in] topSPs Properties of all of the top spacepoints
/// @param[in] middleBottomCounts 1-D array of the number of middle-bottom
///            dublets found for each middle spacepoint
/// @param[in] middleBottomDublets 2-D matrix of size
///            @c nMiddleSPs x @c nBottomSPs, holding the bottom spacepoint
///            indices for the identified middle-bottom dublets
/// @param[in] middleTopCounts 1-D array of the number of middle-top dublets
///            found for each middle spacepoint
/// @param[in] middleTopDublets 2-D matrix of size
///            @c nMiddleSPs x @c nTopSPs, holding the top spacepoint
///            indices for the identified middle-top dublets
/// @param[in] bottomSPLinTransArray 2-dimensional matrix indexed the same way
///            as @c middleBottomArray
/// @param[in] topSPLinTransArray 2-dimensional matrix indexed the same way as
///            @c middleTopArray
/// @param[in] maxScatteringAngle2 Parameter from @c Acts::SeedfinderConfig
/// @param[in] sigmaScattering Parameter from @c Acts::SeedfinderConfig
/// @param[in] minHelixDiameter2 Parameter from @c Acts::SeedfinderConfig
/// @param[in] pT2perRadius Parameter from @c Acts::SeedfinderConfig
/// @param[in] impactMax Parameter from @c Acts::SeedfinderConfig
/// @param[in] impactWeightFactor Parameter from @c Acts::SeedfinderConfig
/// @param[out] tripletsPerBottomDublet 1-dimensional array of the triplet
///             counts for each bottom spacepoint
/// @param[out] tripletIndices 2-dimensional matrix of the indices of the
///             triplets created for each middle-bottom spacepoint dublet
/// @param[out] maxTripletsPerSpB Pointer to the scalar outputting the maximum
///             number of triplets found for any bottom spacepoint dublet
/// @param[out] tripletCount Pointer to the scalar counting the total number of
///             triplets created by the kernel
/// @param[out] triplets 1-dimensional array of all reconstructed triplet
///             candidates
///
__global__ void findTriplets(
    std::size_t middleIndexStart, unsigned int maxMBDublets,
    unsigned int maxMTDublets, unsigned int maxTriplets,
    std::size_t nParallelMiddleSPs, std::size_t nMiddleSPsProcessed,
    std::size_t nBottomSPs, const Details::SpacePoint* bottomSPs,
    std::size_t nMiddleSPs, const Details::SpacePoint* middleSPs,
    std::size_t nTopSPs, const Details::SpacePoint* topSPs,
    const unsigned int* middleBottomCounts,
    const std::size_t* middleBottomDublets, const unsigned int* middleTopCounts,
    const std::size_t* middleTopDublets,
    const Details::LinCircle* bottomSPLinTransArray,
    const Details::LinCircle* topSPLinTransArray, float maxScatteringAngle2,
    float sigmaScattering, float minHelixDiameter2, float pT2perRadius,
    float impactMax, float impactWeightFactor,
    unsigned int* tripletsPerBottomDublet, std::size_t* tripletIndices,
    unsigned int* maxTripletsPerSpB, unsigned int* tripletCount,
    Details::Triplet* triplets) {
  // A sanity check.
  assert(middleIndexStart + nMiddleSPsProcessed <= nMiddleSPs);

  // Find the middle spacepoint index to operate on.
  const unsigned int middleIndexOffset = blockIdx.x * blockDim.x + threadIdx.x;
  if (middleIndexOffset >= nMiddleSPsProcessed) {
    return;
  }
  const unsigned int middleIndex = middleIndexStart + middleIndexOffset;
  assert(middleIndex < nMiddleSPs);

  // Counts of middle-bottom and middle-top pairs for this middle spacepoint.
  const unsigned int middleBottomPairCount = middleBottomCounts[middleIndex];
  const unsigned int middleTopPairCount = middleTopCounts[middleIndex];

  // Find the indices of the middle-bottom and middle-top pairs to operate on.
  const unsigned int tripletCandidateIndex =
      blockIdx.y * blockDim.y + threadIdx.y;
  if (tripletCandidateIndex >= middleBottomPairCount * middleTopPairCount) {
    return;
  }
  const unsigned int bottomDubletIndex =
      tripletCandidateIndex / middleTopPairCount;
  assert(bottomDubletIndex < middleBottomPairCount);
  const unsigned int topDubletIndex =
      tripletCandidateIndex - bottomDubletIndex * middleTopPairCount;
  assert(topDubletIndex < middleTopPairCount);

  // Get the indices of the spacepoints to operate on.
  const unsigned int bottomIndex =
      ACTS_CUDA_MATRIX2D_ELEMENT(middleBottomDublets, nMiddleSPs, nBottomSPs,
                                 middleIndex, bottomDubletIndex);
  assert(bottomIndex < nBottomSPs);
  const unsigned int topIndex = ACTS_CUDA_MATRIX2D_ELEMENT(
      middleTopDublets, nMiddleSPs, nTopSPs, middleIndex, topDubletIndex);
  assert(topIndex < nTopSPs);

  // Load the transformed coordinates of the bottom spacepoint into the thread.
  const Details::LinCircle lb =
      ACTS_CUDA_MATRIX2D_ELEMENT(bottomSPLinTransArray, nMiddleSPs,
                                 maxMBDublets, middleIndex, bottomDubletIndex);

  // 1+(cot^2(theta)) = 1/sin^2(theta)
  float iSinTheta2 = (1. + lb.cotTheta * lb.cotTheta);
  // calculate max scattering for min momentum at the seed's theta angle
  // scaling scatteringAngle^2 by sin^2(theta) to convert pT^2 to p^2
  // accurate would be taking 1/atan(thetaBottom)-1/atan(thetaTop) <
  // scattering
  // but to avoid trig functions we approximate cot by scaling by
  // 1/sin^4(theta)
  // resolving with pT to p scaling --> only divide by sin^2(theta)
  // max approximation error for allowed scattering angles of 0.04 rad at
  // eta=infinity: ~8.5%
  float scatteringInRegion2 = maxScatteringAngle2 * iSinTheta2;
  // multiply the squared sigma onto the squared scattering
  scatteringInRegion2 *= sigmaScattering * sigmaScattering;

  // Load the transformed coordinates of the top spacepoint into the thread.
  const Details::LinCircle lt =
      ACTS_CUDA_MATRIX2D_ELEMENT(topSPLinTransArray, nMiddleSPs, maxMTDublets,
                                 middleIndex, topDubletIndex);

  // Load the parameters of the middle spacepoint into the thread.
  const Details::SpacePoint spM = middleSPs[middleIndex];

  // add errors of spB-spM and spM-spT pairs and add the correlation term
  // for errors on spM
  float error2 =
      lt.Er + lb.Er +
      2 * (lb.cotTheta * lt.cotTheta * spM.varianceR + spM.varianceZ) *
          lb.iDeltaR * lt.iDeltaR;

  float deltaCotTheta = lb.cotTheta - lt.cotTheta;
  float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;
  float dCotThetaMinusError2 = 0.0f;

  // if the error is larger than the difference in theta, no need to
  // compare with scattering
  if (deltaCotTheta2 - error2 > 0) {
    deltaCotTheta = fabs(deltaCotTheta);
    // if deltaTheta larger than the scattering for the lower pT cut, skip
    float error = sqrtf(error2);
    dCotThetaMinusError2 = deltaCotTheta2 + error2 - 2 * deltaCotTheta * error;
    // avoid taking root of scatteringInRegion
    // if left side of ">" is positive, both sides of unequality can be
    // squared
    // (scattering is always positive)
    if (dCotThetaMinusError2 > scatteringInRegion2) {
      return;
    }
  }

  // protects against division by 0
  float dU = lt.U - lb.U;
  if (dU == 0.) {
    return;
  }
  // A and B are evaluated as a function of the circumference parameters
  // x_0 and y_0
  float A = (lt.V - lb.V) / dU;
  float S2 = 1. + A * A;
  float B = lb.V - A * lb.U;
  float B2 = B * B;
  // sqrt(S2)/B = 2 * helixradius
  // calculated radius must not be smaller than minimum radius
  if (S2 < B2 * minHelixDiameter2) {
    return;
  }
  // 1/helixradius: (B/sqrt(S2))/2 (we leave everything squared)
  float iHelixDiameter2 = B2 / S2;
  // calculate scattering for p(T) calculated from seed curvature
  float pT2scatter = 4 * iHelixDiameter2 * pT2perRadius;
  // TODO: include upper pT limit for scatter calc
  // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
  // from rad to deltaCotTheta
  float p2scatter = pT2scatter * iSinTheta2;
  // if deltaTheta larger than allowed scattering for calculated pT, skip
  if ((deltaCotTheta2 - error2 > 0) &&
      (dCotThetaMinusError2 > p2scatter * sigmaScattering * sigmaScattering)) {
    return;
  }
  // A and B allow calculation of impact params in U/V plane with linear
  // function
  // (in contrast to having to solve a quadratic function in x/y plane)
  float Im = fabs((A - B * spM.radius) * spM.radius);

  // Check if the triplet candidate should be accepted.
  if (Im > impactMax) {
    return;
  }

  // Reserve elements (positions) in the global matrices/arrays.
  unsigned int* tripletIndexRowPtr = &(ACTS_CUDA_MATRIX2D_ELEMENT(
      tripletsPerBottomDublet, nParallelMiddleSPs, maxMBDublets,
      middleIndexOffset, bottomDubletIndex));
  const unsigned int tripletIndexRow = atomicAdd(tripletIndexRowPtr, 1);
  assert(tripletIndexRow < maxMTDublets);
  const unsigned int tripletIndex = atomicAdd(tripletCount, 1);
  assert(tripletIndex < maxTriplets);

  // Collect the maximal value of tripletIndexRow + 1 (since we want the
  // count, not the index values) for the next kernel.
  atomicMax(maxTripletsPerSpB, tripletIndexRow + 1);

  // Save the index of the triplet candidate, which will be created now.
  ACTS_CUDA_MATRIX3D_ELEMENT(tripletIndices, nParallelMiddleSPs, maxMBDublets,
                             maxMTDublets, middleIndexOffset, bottomDubletIndex,
                             tripletIndexRow) = tripletIndex;

  // Now store the triplet in the above mentioned location.
  Details::Triplet triplet = {bottomIndex,   middleIndex,
                              topIndex,      Im,
                              B / sqrtf(S2), -(Im * impactWeightFactor)};
  triplets[tripletIndex] = triplet;

  return;
}

/// Kernel performing the "2 fixed spacepoint filtering" of the triplets
///
/// @param[in] seedWeight Pointer to the user-provided seed weight calculating
///            function
/// @param[in] singleSeedCut Pointer to the user-provided seed filtering
///            function
/// @param[in] middleIndexStart The middle spacepoint index that the kernel was
///            "started from"
/// @param[in] maxMBDublets The maximal number of middle-bottom dublets found
///            for any middle spacepoint
/// @param[in] maxMTDublets The maximal number of middle-top dublets found for
///            any middle spacepoint
/// @param[in] maxTriplets The maximum number of triplets for which memory is
///            booked
/// @param[in] nAllTriplets The number of triplets that were reconstructed for
///            this middle spacepoint group
/// @param[in] nParallelMiddleSPs The number of middle spacepoints that the
///            "largest" kernels may be started on in parallel
/// @param[in] nMiddleSPsProcessed The number of middle spacepoints that the
///            kernel was started on in parallel
/// @param[in] middleBottomCounts 1-D array of the number of middle-bottom
///            dublets found for each middle spacepoint
/// @param[in] nBottomSPs The number of bottom spacepoints in @c bottomSPs
/// @param[in] bottomSPs Properties of all of the bottom spacepoints
/// @param[in] nMiddleSPs The number of middle spacepoints in @c middleSPs
/// @param[in] middleSPs Properties of all of the middle spacepoints
/// @param[in] nTopSPs The number of top spacepoints in @c topSPs
/// @param[in] topSPs Properties of all of the top spacepoints
/// @param[in] tripletsPerBottomDublet 1-dimensional array of the triplet
///            counts for each bottom spacepoint
/// @param[in] tripletIndices 2-dimensional matrix of the indices of the
///            triplets created for each middle-bottom spacepoint dublet
/// @param[in] allTriplets 1-dimensional array of all the found triplets
/// @param[in] deltaInvHelixDiameter Parameter from @c Acts::SeedFilterConfig
/// @param[in] deltaRMin Parameter from @c Acts::SeedFilterConfig
/// @param[in] compatSeedWeight Parameter from @c Acts::SeedFilterConfig
/// @param[in] compatSeedLimit Parameter from @c Acts::SeedFilterConfig
/// @param[out] nFilteredTriplets Pointer to the scalar counting all triplets
///             that survive this filter
/// @param[out] filteredTriplets 1-dimensional array of triplets that survive
///             this filter
///
__global__ void filterTriplets2Sp(
    TripletFilterConfig::seedWeightFunc_t seedWeight,
    TripletFilterConfig::singleSeedCutFunc_t singleSeedCut,
    std::size_t middleIndexStart, unsigned int maxMBDublets,
    unsigned int maxMTDublets, unsigned int maxTriplets,
    unsigned int nAllTriplets, std::size_t nParallelMiddleSPs,
    std::size_t nMiddleSPsProcessed, unsigned int* middleBottomCounts,
    std::size_t nBottomSPs, const Details::SpacePoint* bottomSPs,
    std::size_t nMiddleSPs, const Details::SpacePoint* middleSPs,
    std::size_t nTopSPs, const Details::SpacePoint* topSPs,
    const unsigned int* tripletsPerBottomDublet,
    const std::size_t* tripletIndices, const Details::Triplet* allTriplets,
    float deltaInvHelixDiameter, float deltaRMin, float compatSeedWeight,
    std::size_t compatSeedLimit, unsigned int* nFilteredTriplets,
    Details::Triplet* filteredTriplets) {
  // Sanity checks.
  assert(seedWeight != nullptr);
  assert(singleSeedCut != nullptr);
  assert(middleIndexStart + nMiddleSPsProcessed <= nMiddleSPs);

  // Find the middle spacepoint index to operate on.
  const unsigned int middleIndexOffset = blockIdx.x * blockDim.x + threadIdx.x;
  if (middleIndexOffset >= nMiddleSPsProcessed) {
    return;
  }
  const unsigned int middleIndex = middleIndexStart + middleIndexOffset;
  assert(middleIndex < nMiddleSPs);

  // Find the middle-bottom dublet to operate on.
  const unsigned int middleBottomPairCount = middleBottomCounts[middleIndex];
  const unsigned int bottomDubletIndex = blockIdx.y * blockDim.y + threadIdx.y;
  if (bottomDubletIndex >= middleBottomPairCount) {
    return;
  }

  // Find the triplet to operate on.
  const unsigned int nTripletsForMiddleBottom = ACTS_CUDA_MATRIX2D_ELEMENT(
      tripletsPerBottomDublet, nParallelMiddleSPs, maxMBDublets,
      middleIndexOffset, bottomDubletIndex);
  const unsigned int tripletCandidateIndex =
      blockIdx.z * blockDim.z + threadIdx.z;
  if (tripletCandidateIndex >= nTripletsForMiddleBottom) {
    return;
  }

  // Get the index of this triplet.
  const std::size_t triplet1Index = ACTS_CUDA_MATRIX3D_ELEMENT(
      tripletIndices, nParallelMiddleSPs, maxMBDublets, maxMTDublets,
      middleIndexOffset, bottomDubletIndex, tripletCandidateIndex);
  assert(triplet1Index < nAllTriplets);

  // Load this triplet into the thread.
  Details::Triplet triplet1 = allTriplets[triplet1Index];

  // Pre-compute some variables.
  float lowerLimitCurv = triplet1.invHelixDiameter - deltaInvHelixDiameter;
  float upperLimitCurv = triplet1.invHelixDiameter + deltaInvHelixDiameter;
  float currentTop_r = topSPs[triplet1.topIndex].radius;

  // Allow only a maximum number of top spacepoints in the filtering. Since a
  // limit is coming from @c compatSeedLimit anyway, this could potentially be
  // re-written with an array alocation, instead of statically defining the
  // array's size.
  static constexpr std::size_t MAX_TOP_SP = 10;
  assert(compatSeedLimit < MAX_TOP_SP);
  float compatibleSeedR[MAX_TOP_SP];
  std::size_t nCompatibleSeedR = 0;

  // Loop over all the other triplets found for this bottom-middle dublet.
  for (std::size_t i = 0; i < nTripletsForMiddleBottom; ++i) {
    // Don't consider the same triplet that the thread is evaluating in the
    // first place.
    if (i == tripletCandidateIndex) {
      continue;
    }
    // Get the index of the second triplet.
    const std::size_t triplet2Index = ACTS_CUDA_MATRIX3D_ELEMENT(
        tripletIndices, nParallelMiddleSPs, maxMBDublets, maxMTDublets,
        middleIndexOffset, bottomDubletIndex, i);
    assert(triplet2Index < nAllTriplets);
    assert(triplet2Index != triplet1Index);

    // Load the second triplet into the thread.
    const Details::Triplet triplet2 = allTriplets[triplet2Index];
    assert(triplet1.bottomIndex == triplet2.bottomIndex);

    // compared top SP should have at least deltaRMin distance
    float otherTop_r = topSPs[triplet2.topIndex].radius;
    float deltaR = currentTop_r - otherTop_r;
    if (fabs(deltaR) < deltaRMin) {
      continue;
    }

    // curvature difference within limits?
    // TODO: how much slower than sorting all vectors by curvature
    // and breaking out of loop? i.e. is vector size large (e.g. in jets?)
    if (triplet2.invHelixDiameter < lowerLimitCurv) {
      continue;
    }
    if (triplet2.invHelixDiameter > upperLimitCurv) {
      continue;
    }

    bool newCompSeed = true;
    for (std::size_t k = 0; k < nCompatibleSeedR; ++k) {
      // original ATLAS code uses higher min distance for 2nd found compatible
      // seed (20mm instead of 5mm)
      // add new compatible seed only if distance larger than rmin to all
      // other compatible seeds
      if (fabs(compatibleSeedR[k] - otherTop_r) < deltaRMin) {
        newCompSeed = false;
        break;
      }
    }
    if (newCompSeed) {
      compatibleSeedR[nCompatibleSeedR++] = otherTop_r;
      assert(nCompatibleSeedR < MAX_TOP_SP);
      triplet1.weight += compatSeedWeight;
    }
    if (nCompatibleSeedR >= compatSeedLimit) {
      break;
    }
  }

  // Decide whether to keep the triplet or not.
  triplet1.weight +=
      seedWeight(bottomSPs[triplet1.bottomIndex], middleSPs[middleIndex],
                 topSPs[triplet1.topIndex]);
  if (!singleSeedCut(triplet1.weight, bottomSPs[triplet1.bottomIndex],
                     middleSPs[middleIndex], topSPs[triplet1.topIndex])) {
    return;
  }

  // Put the triplet into the "filtered list".
  const unsigned int tripletRow = atomicAdd(nFilteredTriplets, 1);
  assert(tripletRow < nAllTriplets);
  filteredTriplets[tripletRow] = triplet1;
  return;
}

}  // namespace Kernels

namespace Details {

std::vector<std::vector<Triplet>> findTriplets(
    const Info::Device& device, std::size_t maxBlockSize,
    const DubletCounts& dubletCounts, const SeedFilterConfig& seedConfig,
    const TripletFilterConfig& filterConfig, std::size_t nBottomSPs,
    const device_array<SpacePoint>& bottomSPs, std::size_t nMiddleSPs,
    const device_array<SpacePoint>& middleSPs, std::size_t nTopSPs,
    const device_array<SpacePoint>& topSPs,
    const device_array<unsigned int>& middleBottomCounts,
    const device_array<std::size_t>& middleBottomDublets,
    const device_array<unsigned int>& middleTopCounts,
    const device_array<std::size_t>& middleTopDublets,
    float maxScatteringAngle2, float sigmaScattering, float minHelixDiameter2,
    float pT2perRadius, float impactMax) {
  //Acts::Cuda::Nmm::MemoryResource::CudaMemoryResource upstream;
  //Acts::Cuda::Nmm::MemoryResource::ArenaMemoryResource< Acts::Cuda::Nmm::MemoryResource::CudaMemoryResource > arenaMr{&upstream};
  //Acts::Cuda::Nmm::MemoryResource::detail::setCurrentDeviceResource(&arenaMr);
  //Acts::Cuda::Nmm::MemoryResource::DeviceMemoryResource *dmr = Acts::Cuda::Nmm::MemoryResource::detail::getCurrentDeviceResource();
  
  // Calculate the parallelisation for the parameter transformation.
  const int numBlocksLT =
      (dubletCounts.nDublets + maxBlockSize - 1) / maxBlockSize;

  // Create the arrays holding the linear transformed spacepoint parameters.
  auto bottomSPLinTransArray =
      make_device_array<LinCircle>(nMiddleSPs * dubletCounts.maxMBDublets);
  auto topSPLinTransArray =
      make_device_array<LinCircle>(nMiddleSPs * dubletCounts.maxMTDublets);

  // Launch the coordinate transformations.
  Kernels::transformCoordinates<<<numBlocksLT, maxBlockSize>>>(
      dubletCounts.nDublets, dubletCounts.maxMBDublets,
      dubletCounts.maxMTDublets, nBottomSPs, bottomSPs.get(), nMiddleSPs,
      middleSPs.get(), nTopSPs, topSPs.get(), middleBottomCounts.get(),
      middleBottomDublets.get(), middleTopCounts.get(), middleTopDublets.get(),
      bottomSPLinTransArray.get(), topSPLinTransArray.get());
  ACTS_CUDA_ERROR_CHECK(cudaGetLastError());
  ACTS_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

  // With the information from @c Acts::Cuda::Details::DubletCounts, figure out
  // how many middle spacepoints we could handle at the same time in the triplet
  // finding/filtering.

  // For one middle spacepoint we need the following amount:
  const std::size_t memorySizePerMiddleSP =
      // First let's consider the storage of the triplet objects themselves.
      2 * dubletCounts.maxTriplets * sizeof(Triplet) +
      // Then the objects holding indices to the triplets per middle-bottom
      // dublet.
      dubletCounts.maxMBDublets * sizeof(unsigned int) +
      dubletCounts.maxMBDublets * dubletCounts.maxMTDublets *
          sizeof(std::size_t) +
      // Finally the array holding the filtered triplet counts per middle
      // spacepoint.
      sizeof(unsigned int);

  // See how many we can fit into the (still) available memory.
  const std::size_t nParallelMiddleSPs =
      std::min(MemoryManager::instance().availableMemory(device.id) /
                   memorySizePerMiddleSP,
               nMiddleSPs);
  assert(nParallelMiddleSPs > 0);

  // Helper variables for handling the various object counts in device memory.
  enum ObjectCountType : int {
    AllTriplets = 0,        ///< All viable triplets
    FilteredTriplets = 1,   ///< Triplets after the "2SpFixed" filtering
    MaxTripletsPerSpB = 2,  ///< Maximal number of triplets found per SpB
    NObjectCountTypes = 3   ///< The number of different object/counter types
  };

  // Set up the object counters in device memory. The host array is only used to
  // reset the device memory before every iteration.
  auto objectCountsHostNull = make_host_array<unsigned int>(NObjectCountTypes);
  memset(objectCountsHostNull.get(), 0,
         NObjectCountTypes * sizeof(unsigned int));
  auto objectCountsHost = make_host_array<unsigned int>(NObjectCountTypes);
  auto objectCounts = make_device_array<unsigned int>(NObjectCountTypes);

  // Allocate enough memory for triplet candidates that would suffice for every
  // middle spacepoint.
  auto allTriplets =
      make_device_array<Triplet>(nParallelMiddleSPs * dubletCounts.maxTriplets);
  auto filteredTriplets =
      make_device_array<Triplet>(nParallelMiddleSPs * dubletCounts.maxTriplets);
  auto filteredTripletsHost =
      make_host_array<Triplet>(nParallelMiddleSPs * dubletCounts.maxTriplets);

  // Allocate and initialise the array holding the per bottom dublet triplet
  // numbers.
  auto tripletsPerBottomDubletHost = make_host_array<unsigned int>(
      nParallelMiddleSPs * dubletCounts.maxMBDublets);
  memset(tripletsPerBottomDubletHost.get(), 0,
         nParallelMiddleSPs * dubletCounts.maxMBDublets * sizeof(unsigned int));
  auto tripletsPerBottomDublet = make_device_array<unsigned int>(
      nParallelMiddleSPs * dubletCounts.maxMBDublets);

  // Allocate the array holding the indices of the triplets found for a given
  // bottom-middle spacepoint combination.
  auto tripletIndices = make_device_array<std::size_t>(
      nParallelMiddleSPs * dubletCounts.maxMBDublets *
      dubletCounts.maxMTDublets);

  // Allocate and initialise the arrays holding the per-middle-spacepoint
  // filtered triplet counts.
  auto filteredTripletCountsHostNull =
      make_host_array<unsigned int>(nParallelMiddleSPs);
  memset(filteredTripletCountsHostNull.get(), 0,
         nParallelMiddleSPs * sizeof(unsigned int));
  auto filteredTripletCountsHost =
      make_host_array<unsigned int>(nParallelMiddleSPs);
  auto filteredTripletCounts =
      make_device_array<unsigned int>(nParallelMiddleSPs);

  // Block size used in the triplet finding.
  const std::size_t blockSize = std::sqrt(maxBlockSize);

  // Create the result object.
  std::vector<std::vector<Triplet>> result(nMiddleSPs);

  // Copy the dublet counts back to the host.
  auto middleBottomCountsHost = make_host_array<unsigned int>(nMiddleSPs);
  copyToHost(middleBottomCountsHost, middleBottomCounts, nMiddleSPs);
  auto middleTopCountsHost = make_host_array<unsigned int>(nMiddleSPs);
  copyToHost(middleTopCountsHost, middleTopCounts, nMiddleSPs);

  // Execute the triplet finding and filtering in the maximal allowed groups of
  // middle spacepoints.
  for (std::size_t middleIndex = 0; middleIndex < nMiddleSPs;
       middleIndex += nParallelMiddleSPs) {
    // Reset the device arrays.
    copyToDevice(objectCounts, objectCountsHostNull, NObjectCountTypes);
    copyToDevice(tripletsPerBottomDublet, tripletsPerBottomDubletHost,
                 nParallelMiddleSPs * dubletCounts.maxMBDublets);

    // The number of middle spacepoints to process in this iteration.
    const std::size_t nMiddleSPsProcessed =
        std::min(nParallelMiddleSPs, nMiddleSPs - middleIndex);

    // Calculate the parallelisation for the triplet finding for this collection
    // of middle spacepoints.
    const dim3 blockSizeFT(1, maxBlockSize);
    const dim3 numBlocksFT(
        (nMiddleSPsProcessed + blockSizeFT.x - 1) / blockSizeFT.x,
        (dubletCounts.maxTriplets + blockSizeFT.y - 1) / blockSizeFT.y);
    assert(dubletCounts.maxTriplets > 0);

    // Launch the triplet finding for this middle spacepoint.
    Kernels::findTriplets<<<numBlocksFT, blockSizeFT>>>(
        // Parameters needed to use all the arrays.
        middleIndex, dubletCounts.maxMBDublets, dubletCounts.maxMTDublets,
        dubletCounts.maxTriplets, nParallelMiddleSPs, nMiddleSPsProcessed,
        // Parameters of all of the spacepoints.
        nBottomSPs, bottomSPs.get(), nMiddleSPs, middleSPs.get(), nTopSPs,
        topSPs.get(),
        // Arrays describing the identified dublets.
        middleBottomCounts.get(), middleBottomDublets.get(),
        middleTopCounts.get(), middleTopDublets.get(),
        // The transformed parameters of the bottom and top spacepoints for
        // spacepoints taking part in dublets.
        bottomSPLinTransArray.get(), topSPLinTransArray.get(),
        // Configuration constants.
        maxScatteringAngle2, sigmaScattering, minHelixDiameter2, pT2perRadius,
        impactMax, seedConfig.impactWeightFactor,
        // Variables storing the results of the triplet finding.
        tripletsPerBottomDublet.get(), tripletIndices.get(),
        objectCounts.get() + MaxTripletsPerSpB,
        objectCounts.get() + AllTriplets, allTriplets.get());
    ACTS_CUDA_ERROR_CHECK(cudaGetLastError());
    ACTS_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

    // Retrieve the object counts.
    copyToHost(objectCountsHost, objectCounts, NObjectCountTypes);
    const unsigned int nAllTriplets = objectCountsHost.get()[AllTriplets];
    const unsigned int nMaxTripletsPerSpB =
        objectCountsHost.get()[MaxTripletsPerSpB];

    // If no triplet has been found, stop here for this middle spacepoint range.
    if (nAllTriplets == 0) {
      continue;
    }

    // Calculate the parallelisation for the "2SpFixed" filtering of the
    // triplets.
    const dim3 blockSizeF2SP(1, blockSize, blockSize);
    const dim3 numBlocksF2SP(
        (nMiddleSPsProcessed + blockSizeF2SP.x - 1) / blockSizeF2SP.x,
        (dubletCounts.maxMBDublets + blockSizeF2SP.y - 1) / blockSizeF2SP.y,
        (nMaxTripletsPerSpB + blockSizeF2SP.z - 1) / blockSizeF2SP.z);
    assert(dubletCounts.maxMBDublets > 0);
    assert(nMaxTripletsPerSpB > 0);

    // Launch the "2SpFixed" filtering of the triplets.
    assert(filterConfig.seedWeight != nullptr);
    assert(filterConfig.singleSeedCut != nullptr);
    Kernels::filterTriplets2Sp<<<numBlocksF2SP, blockSizeF2SP>>>(
        // Pointers to the user provided filter functions.
        filterConfig.seedWeight, filterConfig.singleSeedCut,
        // Parameters needed to use all the arrays.
        middleIndex, dubletCounts.maxMBDublets, dubletCounts.maxMTDublets,
        dubletCounts.maxTriplets, nAllTriplets, nParallelMiddleSPs,
        nMiddleSPsProcessed, middleBottomCounts.get(),
        // Parameters of all of the spacepoints.
        nBottomSPs, bottomSPs.get(), nMiddleSPs, middleSPs.get(), nTopSPs,
        topSPs.get(),
        // Variables holding the results of the triplet finding.
        tripletsPerBottomDublet.get(), tripletIndices.get(), allTriplets.get(),
        // Configuration constants.
        seedConfig.deltaInvHelixDiameter, seedConfig.deltaRMin,
        seedConfig.compatSeedWeight, seedConfig.compatSeedLimit,
        // Variables storing the results of the filtering.
        objectCounts.get() + FilteredTriplets, filteredTriplets.get());
    ACTS_CUDA_ERROR_CHECK(cudaGetLastError());
    ACTS_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

    // Retrieve the result counts of the filtering.
    copyToHost(objectCountsHost, objectCounts, NObjectCountTypes);

    // The number of triplets that survived the 2Sp filtering.
    const unsigned int nFilteredTriplets =
        objectCountsHost.get()[FilteredTriplets];
    if (nFilteredTriplets == 0) {
      continue;
    }

    // Move the filtered triplets back to the host for the final selection.
    ACTS_CUDA_ERROR_CHECK(cudaMemcpy(
        filteredTripletsHost.get(), filteredTriplets.get(),
        nFilteredTriplets * sizeof(Triplet), cudaMemcpyDeviceToHost));

    // Fill the output variable.
    for (std::size_t i = 0; i < nFilteredTriplets; ++i) {
      // Access the triplet.
      const Triplet& triplet = filteredTripletsHost.get()[i];
      // Put it into the output object.
      result[triplet.middleIndex].push_back(triplet);
    }
  }
  
  // Free host array
  /*
  Cuda::Nmm::MemoryResource::HostMemoryResource *hmr = Cuda::getCPUmmr();
  hmr->deallocate(static_cast<void*>(objectCountsHostNull.get()), NObjectCountTypes);
  hmr->deallocate(static_cast<void*>(objectCountsHost.get()), NObjectCountTypes);
  hmr->deallocate(static_cast<void*>(filteredTripletsHost.get()), nParallelMiddleSPs * dubletCounts.maxTriplets);
  hmr->deallocate(static_cast<void*>(tripletsPerBottomDubletHost.get()), nParallelMiddleSPs * dubletCounts.maxMBDublets);
  hmr->deallocate(static_cast<void*>(filteredTripletCountsHostNull.get()), nParallelMiddleSPs);
  hmr->deallocate(static_cast<void*>(filteredTripletCountsHost.get()), nParallelMiddleSPs);
  hmr->deallocate(static_cast<void*>(middleBottomCountsHost.get()), nMiddleSPs);
  hmr->deallocate(static_cast<void*>(middleTopCountsHost.get()), nMiddleSPs);
  */
  // Return the indices of all identified triplets.
  assert(result.size() == nMiddleSPs);
  return result;
}

}  // namespace Details
}  // namespace Cuda
}  // namespace Acts
