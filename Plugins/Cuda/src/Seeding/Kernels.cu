// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Cuda/Cuda.hpp"
#include "../Utilities/NewMemoryManager.cu"
#include "Acts/Plugins/Cuda/Seeding/Kernels.cuh"

#include <iostream>

#include <cuda.h>
#include <cuda_runtime.h>

__global__ void cuSearchDoublet(
    const int* nSpM, const float* spMmat, const int* nSpB, const float* spBmat,
    const int* nSpT, const float* spTmat, const float* deltaRMin,
    const float* deltaRMax, const float* cotThetaMax,
    const float* collisionRegionMin, const float* collisionRegionMax,
    int* nSpMcomp, int* nSpBcompPerSpM_Max, int* nSpTcompPerSpM_Max,
    int* nSpBcompPerSpM, int* nSpTcompPerSpM, int* McompIndex, int* BcompIndex,
    int* tmpBcompIndex, int* TcompIndex, int* tmpTcompIndex);

__global__ void cuTransformCoordinate(
    const int* nSpM, const float* spMmat, const int* McompIndex,
    const int* nSpB, const float* spBmat, const int* nSpBcompPerSpM_Max,
    const int* BcompIndex, const int* nSpT, const float* spTmat,
    const int* nSpTcompPerSpM_Max, const int* TcompIndex, float* spMcompMat,
    float* spBcompMatPerSpM, float* circBcompMatPerSpM, float* spTcompMatPerSpM,
    float* circTcompMatPerSpM);

__device__ void sp2circle(bool isBottom, const float* spM, const float* spB,
                          float* circB);

__global__ void cuSearchTriplet(
    const int* nSpTcompPerSpM, const int* nSpMcomp, const float* spMcompMat,
    const int* nSpBcompPerSpM_Max, const int* BcompIndex,
    const float* circBcompMatPerSpM, const int* nSpTcompPerSpM_Max,
    const int* TcompIndex, const float* spTcompMatPerSpM,
    const float* circTcompMatPerSpM, const float* maxScatteringAngle2,
    const float* sigmaScattering, const float* minHelixDiameter2,
    const float* pT2perRadius, const float* impactMax,
    const int* nTrplPerSpMLimit, const int* nTrplPerSpBLimit,
    const float* deltaInvHelixDiameter, const float* impactWeightFactor,
    const float* deltaRMin, const float* compatSeedWeight,
    const size_t* compatSeedLimit, int* nTrplPerSpM, Triplet* TripletsPerSpM);

namespace Acts {

void searchDoublet(const dim3 grid, const dim3 block, const int* nSpM,
                   const float* spMmat, const int* nSpB, const float* spBmat,
                   const int* nSpT, const float* spTmat, const float* deltaRMin,
                   const float* deltaRMax, const float* cotThetaMax,
                   const float* collisionRegionMin,
                   const float* collisionRegionMax, int* nSpMcomp,
                   int* nSpBcompPerSpM_Max, int* nSpTcompPerSpM_Max,
                   int* nSpBcompPerSpM, int* nSpTcompPerSpM, int* McompIndex,
                   int* BcompIndex, int* tmpBcompIndex, int* TcompIndex,
                   int* tmpTcompIndex) {
  if (grid.x == 0) {
    return;
  }
  // int sharedMemSize = (2*sizeof(int))*block.x + 2*sizeof(int);
  int sharedMemSize = 2 * sizeof(int);

  cuSearchDoublet<<<grid, block, sharedMemSize>>>(
      nSpM, spMmat, nSpB, spBmat, nSpT, spTmat, deltaRMin, deltaRMax,
      cotThetaMax, collisionRegionMin, collisionRegionMax, nSpMcomp,
      nSpBcompPerSpM_Max, nSpTcompPerSpM_Max, nSpBcompPerSpM, nSpTcompPerSpM,
      McompIndex, BcompIndex, tmpBcompIndex, TcompIndex, tmpTcompIndex);
  ACTS_CUDA_ERROR_CHECK(cudaGetLastError());
}

void transformCoordinate(const dim3 grid, const dim3 block, const int* nSpM,
                         const float* spMmat, const int* McompIndex,
                         const int* nSpB, const float* spBmat,
                         const int* nSpBcompPerSpM_Max, const int* BcompIndex,
                         const int* nSpT, const float* spTmat,
                         const int* nSpTcompPerSpM_Max, const int* TcompIndex,
                         float* spMcompMat, float* spBcompMatPerSpM,
                         float* circBcompMatPerSpM, float* spTcompMatPerSpM,
                         float* circTcompMatPerSpM) {
  if (grid.x == 0) {
    return;
  }
  int sharedMemSize = sizeof(float) * 6;
  cuTransformCoordinate<<<grid, block, sharedMemSize>>>(
      nSpM, spMmat, McompIndex, nSpB, spBmat, nSpBcompPerSpM_Max, BcompIndex,
      nSpT, spTmat, nSpTcompPerSpM_Max, TcompIndex, spMcompMat,
      spBcompMatPerSpM, circBcompMatPerSpM, spTcompMatPerSpM,
      circTcompMatPerSpM);
  ACTS_CUDA_ERROR_CHECK(cudaGetLastError());
}

void searchTriplet(
    const dim3 grid, const dim3 block, const int* nSpTcompPerSpM_cpu,
    const int* nSpTcompPerSpM_cuda, const int* nSpMcomp,
    const float* spMcompMat, const int* nSpBcompPerSpM_Max,
    const int* BcompIndex, const float* circBcompMatPerSpM,
    const int* nSpTcompPerSpM_Max, const int* TcompIndex,
    const float* spTcompMatPerSpM, const float* circTcompMatPerSpM,
    const float* maxScatteringAngle2, const float* sigmaScattering,
    const float* minHelixDiameter2, const float* pT2perRadius,
    const float* impactMax, const int* nTrplPerSpMLimit,
    const int* nTrplPerSpBLimit_cpu, const int* nTrplPerSpBLimit_cuda,
    const float* deltaInvHelixDiameter, const float* impactWeightFactor,
    const float* deltaRMin, const float* compatSeedWeight,
    const size_t* compatSeedLimit_cpu, const size_t* compatSeedLimit_cuda,
    int* nTrplPerSpM, Triplet* TripletsPerSpM, cudaStream_t* stream) {
  if (grid.x == 0) {
    return;
  }
  int sharedMemSize = sizeof(Triplet) * (*nTrplPerSpBLimit_cpu);
  sharedMemSize += sizeof(float) * (*compatSeedLimit_cpu);
  sharedMemSize += sizeof(int);

  cuSearchTriplet<<<grid, block, sharedMemSize, *stream>>>(
      nSpTcompPerSpM_cuda, nSpMcomp, spMcompMat, nSpBcompPerSpM_Max, BcompIndex,
      circBcompMatPerSpM, nSpTcompPerSpM_Max, TcompIndex, spTcompMatPerSpM,
      circTcompMatPerSpM, maxScatteringAngle2, sigmaScattering,
      minHelixDiameter2, pT2perRadius, impactMax, nTrplPerSpMLimit,
      nTrplPerSpBLimit_cuda, deltaInvHelixDiameter, impactWeightFactor,
      deltaRMin, compatSeedWeight, compatSeedLimit_cuda, nTrplPerSpM,
      TripletsPerSpM);
  ACTS_CUDA_ERROR_CHECK(cudaGetLastError());
}

}  // namespace Acts

__global__ void cuSearchDoublet(
    const int* nSpM, const float* spMmat, const int* nSpB, const float* spBmat,
    const int* nSpT, const float* spTmat, const float* deltaRMin,
    const float* deltaRMax, const float* cotThetaMax,
    const float* collisionRegionMin, const float* collisionRegionMax,
    int* nSpMcomp, int* nSpBcompPerSpM_Max, int* nSpTcompPerSpM_Max,
    int* nSpBcompPerSpM, int* nSpTcompPerSpM, int* McompIndex, int* BcompIndex,
    int* tmpBcompIndex, int* TcompIndex, int* tmpTcompIndex) {
  extern __shared__ float sharedMem[];
  int* mPos = (int*)sharedMem;
  int* isMcompat = (int*)&mPos[1];

  if (threadIdx.x == 0) {
    *isMcompat = false;
  }
  __syncthreads();

  float rM = spMmat[blockIdx.x + (*nSpM) * 3];
  float zM = spMmat[blockIdx.x + (*nSpM) * 2];

  bool isBcompat(true);
  bool isTcompat(true);

  int offset(0);

  while (offset < max(*nSpB, *nSpT)) {
    isBcompat = true;

    // Doublet search for bottom hits
    if (threadIdx.x + offset < *nSpB) {
      float rB = spBmat[threadIdx.x + offset + (*nSpB) * 3];
      float zB = spBmat[threadIdx.x + offset + (*nSpB) * 2];

      float deltaR = rM - rB;
      if (deltaR > *deltaRMax) {
        isBcompat = false;
      }

      if (deltaR < *deltaRMin) {
        isBcompat = false;
      }

      float cotTheta = (zM - zB) / deltaR;
      if (fabsf(cotTheta) > *cotThetaMax) {
        isBcompat = false;
      }

      float zOrigin = zM - rM * cotTheta;
      if (zOrigin < *collisionRegionMin || zOrigin > *collisionRegionMax) {
        isBcompat = false;
      }

      if (isBcompat == true) {
        int bPos = atomicAdd(&nSpBcompPerSpM[blockIdx.x], 1);
        tmpBcompIndex[bPos + (*nSpB) * blockIdx.x] = threadIdx.x + offset;
      }
    }

    isTcompat = true;

    // Doublet search for top hits
    if (threadIdx.x + offset < *nSpT) {
      float rT = spTmat[threadIdx.x + offset + (*nSpT) * 3];
      float zT = spTmat[threadIdx.x + offset + (*nSpT) * 2];
      float deltaR = rT - rM;
      if (deltaR < *deltaRMin) {
        isTcompat = false;
      }

      if (deltaR > *deltaRMax) {
        isTcompat = false;
      }

      if (isTcompat == true) {
        float cotTheta = (zT - zM) / deltaR;
        if (fabsf(cotTheta) > *cotThetaMax) {
          isTcompat = false;
        }

        float zOrigin = zM - rM * cotTheta;
        if (zOrigin < *collisionRegionMin || zOrigin > *collisionRegionMax) {
          isTcompat = false;
        }
      }

      if (isTcompat == true) {
        int tPos = atomicAdd(&nSpTcompPerSpM[blockIdx.x], 1);
        tmpTcompIndex[tPos + (*nSpT) * blockIdx.x] = threadIdx.x + offset;
      }
    }

    offset += blockDim.x;
  }

  __syncthreads();

  if (threadIdx.x == 0) {
    if (nSpBcompPerSpM[blockIdx.x] > 0 && nSpTcompPerSpM[blockIdx.x] > 0) {
      *mPos = atomicAdd(nSpMcomp, 1);
      *isMcompat = true;
      McompIndex[*mPos] = blockIdx.x;

      int bMax = atomicMax(nSpBcompPerSpM_Max, nSpBcompPerSpM[blockIdx.x]);
      int tMax = atomicMax(nSpTcompPerSpM_Max, nSpTcompPerSpM[blockIdx.x]);
    }
  }

  __syncthreads();

  if (*isMcompat == true) {
    offset = 0;
    while (offset <
           max(nSpBcompPerSpM[blockIdx.x], nSpTcompPerSpM[blockIdx.x])) {
      if (threadIdx.x + offset < nSpBcompPerSpM[blockIdx.x]) {
        BcompIndex[threadIdx.x + offset + (*nSpB) * (*mPos)] =
            tmpBcompIndex[threadIdx.x + offset + (*nSpB) * blockIdx.x];
      }

      if (threadIdx.x + offset < nSpTcompPerSpM[blockIdx.x]) {
        TcompIndex[threadIdx.x + offset + (*nSpT) * (*mPos)] =
            tmpTcompIndex[threadIdx.x + offset + (*nSpT) * blockIdx.x];
      }
      offset += blockDim.x;
    }
  }
}

__global__ void cuTransformCoordinate(
    const int* nSpM, const float* spMmat, const int* McompIndex,
    const int* nSpB, const float* spBmat, const int* nSpBcompPerSpM_Max,
    const int* BcompIndex, const int* nSpT, const float* spTmat,
    const int* nSpTcompPerSpM_Max, const int* TcompIndex, float* spMcompMat,
    float* spBcompMatPerSpM, float* circBcompMatPerSpM, float* spTcompMatPerSpM,
    float* circTcompMatPerSpM) {
  extern __shared__ float spM[];

  if (threadIdx.x == 0) {
    int mIndex = McompIndex[blockIdx.x];
    for (int i = 0; i < 6; i++) {
      spM[i] = spMcompMat[blockIdx.x + gridDim.x * i] =
          spMmat[mIndex + (*nSpM) * i];
    }
  }

  __syncthreads();

  int offset(0);
  while (offset < max(*nSpBcompPerSpM_Max, *nSpTcompPerSpM_Max)) {
    if (threadIdx.x + offset < *nSpBcompPerSpM_Max) {
      float spB[6];
      float circB[6];
      int bIndex = BcompIndex[threadIdx.x + offset + (*nSpB) * blockIdx.x];

      // matrix reduction
      for (int i = 0; i < 6; i++) {
        spB[i] =
            spBcompMatPerSpM[threadIdx.x + offset +
                             (*nSpBcompPerSpM_Max) * (6 * blockIdx.x + i)] =
                spBmat[bIndex + (*nSpB) * i];
      }

      // do transform coordinate (xy->uv)
      sp2circle(true, spM, spB, circB);
      for (int i = 0; i < 6; i++) {
        circBcompMatPerSpM[threadIdx.x + offset +
                           (*nSpBcompPerSpM_Max) * (6 * blockIdx.x + i)] =
            circB[i];
      }
    }

    if (threadIdx.x + offset < *nSpTcompPerSpM_Max) {
      float spT[6];
      float circT[6];
      int tIndex = TcompIndex[threadIdx.x + offset + (*nSpT) * blockIdx.x];

      // matrix reduction
      for (int i = 0; i < 6; i++) {
        spT[i] =
            spTcompMatPerSpM[threadIdx.x + offset +
                             (*nSpTcompPerSpM_Max) * (6 * blockIdx.x + i)] =
                spTmat[tIndex + (*nSpT) * i];
      }

      // do transform coordinate (xy->uv)
      sp2circle(false, spM, spT, circT);
      for (int i = 0; i < 6; i++) {
        circTcompMatPerSpM[threadIdx.x + offset +
                           (*nSpTcompPerSpM_Max) * (6 * blockIdx.x + i)] =
            circT[i];
      }
    }

    offset += blockDim.x;
  }
}

__device__ void sp2circle(bool isBottom, const float* spM, const float* spB,
                          float* circB) {
  float xM = spM[0];
  float yM = spM[1];
  float zM = spM[2];
  float rM = spM[3];
  float varianceRM = spM[4];
  float varianceZM = spM[5];

  float cosPhiM = xM / rM;
  float sinPhiM = yM / rM;

  float xB = spB[0];
  float yB = spB[1];
  float zB = spB[2];
  // float rB = spB[3];
  float varianceRB = spB[4];
  float varianceZB = spB[5];

  float deltaX = xB - xM;
  float deltaY = yB - yM;
  float deltaZ = zB - zM;

  // calculate projection fraction of spM->sp vector pointing in same
  // direction as
  // vector origin->spM (x) and projection fraction of spM->sp vector pointing
  // orthogonal to origin->spM (y)
  float x = deltaX * cosPhiM + deltaY * sinPhiM;
  float y = deltaY * cosPhiM - deltaX * sinPhiM;
  // 1/(length of M -> SP)
  float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
  float iDeltaR = sqrtf(iDeltaR2);

  int bottomFactor = 1 * (int(!(isBottom))) - 1 * (int(isBottom));

  // cot_theta = (deltaZ/deltaR)
  float cot_theta = deltaZ * iDeltaR * bottomFactor;
  // VERY frequent (SP^3) access

  // location on z-axis of this SP-duplet
  float Zo = zM - rM * cot_theta;

  // transformation of circle equation (x,y) into linear equation (u,v)
  // x^2 + y^2 - 2x_0*x - 2y_0*y = 0
  // is transformed into
  // 1 - 2x_0*u - 2y_0*v = 0
  // using the following m_U and m_V
  // (u = A + B*v); A and B are created later on
  float U = x * iDeltaR2;
  float V = y * iDeltaR2;
  // error term for sp-pair without correlation of middle space point
  float Er = ((varianceZM + varianceZB) +
              (cot_theta * cot_theta) * (varianceRM + varianceRB)) *
             iDeltaR2;

  circB[0] = Zo;
  circB[1] = cot_theta;
  circB[2] = iDeltaR;
  circB[3] = Er;
  circB[4] = U;
  circB[5] = V;
}

__global__ void cuSearchTriplet(
    const int* nSpTcompPerSpM, const int* nSpMcomp, const float* spMcompMat,
    const int* nSpBcompPerSpM_Max, const int* BcompIndex,
    const float* circBcompMatPerSpM, const int* nSpTcompPerSpM_Max,
    const int* TcompIndex, const float* spTcompMatPerSpM,
    const float* circTcompMatPerSpM, const float* maxScatteringAngle2,
    const float* sigmaScattering, const float* minHelixDiameter2,
    const float* pT2perRadius, const float* impactMax,
    const int* nTrplPerSpMLimit, const int* nTrplPerSpBLimit,
    const float* deltaInvHelixDiameter, const float* impactWeightFactor,
    const float* deltaRMin, const float* compatSeedWeight,
    const size_t* compatSeedLimit, int* nTrplPerSpM, Triplet* TripletsPerSpM) {
  extern __shared__ Triplet sh[];
  Triplet* triplets = (Triplet*)sh;
  float* compatibleSeedR = (float*)&triplets[*nTrplPerSpBLimit];
  int* nTrplPerSpB = (int*)&compatibleSeedR[*compatSeedLimit];

  if (threadIdx.x == 0) {
    *nTrplPerSpB = 0;
  }

  float rM = spMcompMat[(*nSpMcomp) * 3];
  float varianceRM = spMcompMat[(*nSpMcomp) * 4];
  float varianceZM = spMcompMat[(*nSpMcomp) * 5];
  // Zob values from CPU and CUDA are slightly different
  // float Zob        = circBcompMatPerSpM[blockId+(*nSpBcompPerSpM_Max)*0];
  float cotThetaB = circBcompMatPerSpM[blockIdx.x + (*nSpBcompPerSpM_Max) * 1];
  float iDeltaRB = circBcompMatPerSpM[blockIdx.x + (*nSpBcompPerSpM_Max) * 2];
  float ErB = circBcompMatPerSpM[blockIdx.x + (*nSpBcompPerSpM_Max) * 3];
  float Ub = circBcompMatPerSpM[blockIdx.x + (*nSpBcompPerSpM_Max) * 4];
  float Vb = circBcompMatPerSpM[blockIdx.x + (*nSpBcompPerSpM_Max) * 5];
  float iSinTheta2 = (1. + cotThetaB * cotThetaB);
  float scatteringInRegion2 = (*maxScatteringAngle2) * iSinTheta2;
  scatteringInRegion2 *= (*sigmaScattering) * (*sigmaScattering);

  int offset(0);

  while (offset < *nSpTcompPerSpM) {
    bool isPassed(1);
    float impact;
    float invHelix;
    if (threadIdx.x + offset < *nSpTcompPerSpM) {
      // float Zot        = circTcompMatPerSpM[threadId+(*nSpTcompPerSpM)*0];
      float cotThetaT =
          circTcompMatPerSpM[threadIdx.x + offset + (*nSpTcompPerSpM_Max) * 1];
      float iDeltaRT =
          circTcompMatPerSpM[threadIdx.x + offset + (*nSpTcompPerSpM_Max) * 2];
      float ErT =
          circTcompMatPerSpM[threadIdx.x + offset + (*nSpTcompPerSpM_Max) * 3];
      float Ut =
          circTcompMatPerSpM[threadIdx.x + offset + (*nSpTcompPerSpM_Max) * 4];
      float Vt =
          circTcompMatPerSpM[threadIdx.x + offset + (*nSpTcompPerSpM_Max) * 5];

      // add errors of spB-spM and spM-spT pairs and add the correlation term
      // for errors on spM
      float error2 = ErT + ErB +
                     2 * (cotThetaB * cotThetaT * varianceRM + varianceZM) *
                         iDeltaRB * iDeltaRT;

      float deltaCotTheta = cotThetaB - cotThetaT;
      float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;
      float error;
      float dCotThetaMinusError2;

      // if the error is larger than the difference in theta, no need to
      // compare with scattering
      if (deltaCotTheta2 - error2 > 0) {
        deltaCotTheta = fabsf(deltaCotTheta);
        // if deltaTheta larger than the scattering for the lower pT cut, skip
        error = sqrtf(error2);
        dCotThetaMinusError2 =
            deltaCotTheta2 + error2 - 2 * deltaCotTheta * error;
        // avoid taking root of scatteringInRegion
        // if left side of ">" is positive, both sides of unequality can be
        // squared
        // (scattering is always positive)

        if (dCotThetaMinusError2 > scatteringInRegion2) {
          isPassed = 0;
        }
      }

      // protects against division by 0
      float dU = Ut - Ub;
      if (dU == 0.) {
        isPassed = 0;
      }

      // A and B are evaluated as a function of the circumference parameters
      // x_0 and y_0
      float A = (Vt - Vb) / dU;
      float S2 = 1. + A * A;
      float B = Vb - A * Ub;
      float B2 = B * B;
      // sqrtf(S2)/B = 2 * helixradius
      // calculated radius must not be smaller than minimum radius
      if (S2 < B2 * (*minHelixDiameter2)) {
        isPassed = 0;
      }

      // 1/helixradius: (B/sqrtf(S2))/2 (we leave everything squared)
      float iHelixDiameter2 = B2 / S2;
      // calculate scattering for p(T) calculated from seed curvature
      float pT2scatter = 4 * iHelixDiameter2 * (*pT2perRadius);
      // TODO: include upper pT limit for scatter calc
      // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
      // from rad to deltaCotTheta
      float p2scatter = pT2scatter * iSinTheta2;
      // if deltaTheta larger than allowed scattering for calculated pT, skip
      if ((deltaCotTheta2 - error2 > 0) &&
          (dCotThetaMinusError2 >
           p2scatter * (*sigmaScattering) * (*sigmaScattering))) {
        isPassed = 0;
      }
      // A and B allow calculation of impact params in U/V plane with linear
      // function
      // (in contrast to having to solve a quadratic function in x/y plane)
      impact = fabsf((A - B * rM) * rM);
      invHelix = B / sqrtf(S2);
      if (impact > (*impactMax)) {
        isPassed = 0;
      }
    }

    __syncthreads();

    if (threadIdx.x + offset < *nSpTcompPerSpM) {
      // The index will be different (and not deterministic) becuase of atomic
      // operation It will be resorted after kernel call
      if (isPassed == 1) {
        int tPos = atomicAdd(nTrplPerSpB, 1);

        if (tPos < *nTrplPerSpBLimit) {
          triplets[tPos].weight = 0;
          triplets[tPos].bIndex = BcompIndex[blockIdx.x];
          triplets[tPos].tIndex = TcompIndex[threadIdx.x + offset];
          triplets[tPos].topRadius =
              spTcompMatPerSpM[threadIdx.x + offset +
                               (*nSpTcompPerSpM_Max) * 3];
          triplets[tPos].impactParameter = impact;
          triplets[tPos].invHelixDiameter = invHelix;
        }
      }
    }
    offset += blockDim.x;
  }
  __syncthreads();

  if (threadIdx.x == 0 && *nTrplPerSpB > *nTrplPerSpBLimit) {
    *nTrplPerSpB = *nTrplPerSpBLimit;
  }

  __syncthreads();
  int jj = threadIdx.x;

  // bubble sort tIndex
  for (int i = 0; i < *nTrplPerSpB / 2 + 1; i++) {
    if (threadIdx.x < *nTrplPerSpB) {
      if (jj % 2 == 0 && jj < *nTrplPerSpB - 1) {
        if (triplets[jj + 1].tIndex < triplets[jj].tIndex) {
          Triplet tempVal = triplets[jj];
          triplets[jj] = triplets[jj + 1];
          triplets[jj + 1] = tempVal;
        }
      }
    }
    __syncthreads();
    if (threadIdx.x < *nTrplPerSpB) {
      if (jj % 2 == 1 && jj < *nTrplPerSpB - 1) {
        if (triplets[jj + 1].tIndex < triplets[jj].tIndex) {
          Triplet tempVal = triplets[jj];
          triplets[jj] = triplets[jj + 1];
          triplets[jj + 1] = tempVal;
        }
      }
    }
    __syncthreads();
  }
  __syncthreads();

  // serial algorithm for seed filtering
  // Need to optimize later
  if (threadIdx.x == 0) {
    int nCompatibleSeedR;
    float lowerLimitCurv;
    float upperLimitCurv;
    float deltaR;
    bool newCompSeed;

    for (int i = 0; i < *nTrplPerSpB; i++) {
      nCompatibleSeedR = 0;
      lowerLimitCurv = triplets[i].invHelixDiameter - *deltaInvHelixDiameter;
      upperLimitCurv = triplets[i].invHelixDiameter + *deltaInvHelixDiameter;
      float& currentTop_r = triplets[i].topRadius;
      float& weight = triplets[i].weight;
      weight = -(triplets[i].impactParameter * (*impactWeightFactor));

      for (int j = 0; j < *nTrplPerSpB; j++) {
        if (i == j) {
          continue;
        }
        float& otherTop_r = triplets[j].topRadius;
        deltaR = currentTop_r - otherTop_r;

        if (fabsf(deltaR) < *deltaRMin) {
          continue;
        }
        // curvature difference within limits?
        // TODO: how much slower than sorting all vectors by curvature
        // and breaking out of loop? i.e. is vector size large (e.g. in jets?)
        if (triplets[j].invHelixDiameter < lowerLimitCurv) {
          continue;
        }
        if (triplets[j].invHelixDiameter > upperLimitCurv) {
          continue;
        }
        newCompSeed = true;

        for (int k = 0; k < nCompatibleSeedR; k++) {
          // original ATLAS code uses higher min distance for 2nd found
          // compatible seed (20mm instead of 5mm) add new compatible seed only
          // if distance larger than rmin to all other compatible seeds
          float& previousDiameter = compatibleSeedR[k];
          if (fabsf(previousDiameter - otherTop_r) < *deltaRMin) {
            newCompSeed = false;
            break;
          }
        }

        if (newCompSeed) {
          compatibleSeedR[nCompatibleSeedR] = otherTop_r;
          nCompatibleSeedR++;
          weight += *compatSeedWeight;
        }
        if (nCompatibleSeedR >= *compatSeedLimit) {
          break;
        }
      }

      int pos = atomicAdd(nTrplPerSpM, 1);

      if (pos < *nTrplPerSpMLimit) {
        TripletsPerSpM[pos] = triplets[i];
      }
    }
  }

  if (threadIdx.x == 0 && *nTrplPerSpM > *nTrplPerSpMLimit) {
    *nTrplPerSpM = *nTrplPerSpMLimit;
  }
}
