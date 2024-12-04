# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 3.25)

# CMake include(s).
include(FetchContent)

# Tell the user what's happening.
message(STATUS "Building VecMem as part of the Acts project")

set(VECMEM_VERSION "v${_acts_vecmem_version}")

# Declare where to get VecMem from.
FetchContent_Declare(VecMem SYSTEM ${ACTS_VECMEM_SOURCE})

# Options used in the build of VecMem.
set(VECMEM_BUILD_TESTING
    FALSE
    CACHE BOOL
    "Turn off the build of the VecMem unit tests"
)
set(VECMEM_BUILD_BENCHMARKING
    FALSE
    CACHE BOOL
    "Turn off the build of VecMem benchmarking"
)

# Now set up its build.
FetchContent_MakeAvailable(VecMem)
