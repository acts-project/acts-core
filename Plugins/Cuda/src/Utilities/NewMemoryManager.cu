// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager.hpp"

namespace Acts {

Nmm::MemoryResource::HostMemoryResource* getCPUmmr() {
  // Create a simple CPU memory resource for now
  static Nmm::MemoryResource::HostMemoryResource *hmr = new Nmm::MemoryResource::CPUMemoryResource();
  return hmr;
}

Nmm::MemoryResource::HostMemoryResource* getPinnedmmr() {
  // Create a simple CPU memory resource for now
  static Nmm::MemoryResource::HostMemoryResource *hmr = new Nmm::MemoryResource::PinnedMemoryResource();
  return hmr;
}

Nmm::MemoryResource::DeviceMemoryResource* getGPUmmr() {
  // Create an upstream type, Cuda for now
  static Nmm::MemoryResource::CudaMemoryResource upstream;
  // Create an arena memory resource (in this part of the program will be setteable for the user in whenever it appears more memory resource)
  static Nmm::MemoryResource::ArenaMemoryResource< Nmm::MemoryResource::CudaMemoryResource > arenaMr{&upstream};
  // Set this new memory resource, because in the strat is by default a Cuda memory resource
  static Nmm::MemoryResource::DeviceMemoryResource *dmrOld = Nmm::MemoryResource::detail::setCurrentDeviceResource(&arenaMr);
  // Get the arena memory resource just set
  static Nmm::MemoryResource::DeviceMemoryResource *dmr = Nmm::MemoryResource::detail::getCurrentDeviceResource();
  return dmr;
}

} // namespace Acts