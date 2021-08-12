// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/CudaStreamView.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/detail/Aligned.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/detail/Block.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/detail/Arena.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/detail/ArenaCleaner.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/detail/GlobalArena.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/detail/Methods.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/ArenaMemoryResource.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/CudaMemoryResource.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/PerDeviceResource.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/DeviceMemoryResource.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Host/CPUMemoryResource.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Host/HostMemoryResource.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Host/PinnedMemoryResource.hpp"

namespace Acts {

// return the memory resource for the CPU
Nmm::MemoryResource::HostMemoryResource* getCPUmmr();

// return the memory resource for the GPU
Nmm::MemoryResource::DeviceMemoryResource* getGPUmmr();

// return the memory resource for pinned memory
Nmm::MemoryResource::HostMemoryResource* getPinnedmmr();

} // namespace Acts