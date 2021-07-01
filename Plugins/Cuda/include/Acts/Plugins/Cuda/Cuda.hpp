// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Cuda/Utilities/CpuMatrix.hpp"
#include "Acts/Plugins/Cuda/Utilities/CpuScalar.hpp"
#include "Acts/Plugins/Cuda/Utilities/CpuVector.hpp"
#include "Acts/Plugins/Cuda/Utilities/CudaMatrix.cu"
#include "Acts/Plugins/Cuda/Utilities/CudaScalar.cu"
#include "Acts/Plugins/Cuda/Utilities/CudaUtils.cu"
#include "Acts/Plugins/Cuda/Utilities/CudaVector.cu"
#include "Acts/Plugins/Cuda/Utilities/UsmMatrix.cu"
#include "Acts/Plugins/Cuda/Utilities/UsmScalar.cu"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager.hpp"

namespace Acts {
class Cuda;
}
