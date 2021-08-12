// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/DeviceMemoryResource.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/CudaMemoryResource.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/detail/ArenaCleaner.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/CudaStreamView.hpp"

// CUDA include(s).
#include "cuda_runtime.h"

// System include(s)
#include <map>
#include <shared_mutex>
#include <thread>

namespace Acts {
namespace Nmm {
namespace MemoryResource {



} // namaspace MemoryResource
} // namaspace Nmm
} // namaspace Acts