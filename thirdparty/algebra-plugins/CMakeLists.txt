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
message(STATUS "Building algebra-plugins as part of the Acts project")

FetchContent_Declare(AlgebraPlugins SYSTEM ${ACTS_ALGEBRAPLUGINS_SOURCE})

# Options used in the build of Algebra Plugins.
set(ALGEBRA_PLUGINS_BUILD_TESTING
    FALSE
    CACHE BOOL
    "Turn off the build of the Algebra Plugins unit tests"
)
set(ALGEBRA_PLUGINS_INCLUDE_EIGEN
    TRUE
    CACHE BOOL
    "Turn on the build of algebra::eigen"
)
set(ALGEBRA_PLUGINS_INCLUDE_VC
    FALSE
    CACHE BOOL
    "Turn off the build of algebra::vc_array"
)
set(ALGEBRA_PLUGINS_INCLUDE_VECMEM
    TRUE
    CACHE BOOL
    "Turn on the build of algebra::vecmem_array"
)

# Build smatrix plugin only when it needs to be
# ROOT is not provided in native github compilers :/
set(ALGEBRA_PLUGINS_INCLUDE_SMATRIX
    ${DETRAY_SMATRIX_PLUGIN}
    CACHE BOOL
    "Turn on/off the build of algebra::smatrix"
)

set(ALGEBRA_PLUGINS_SETUP_EIGEN3
    OFF
    CACHE BOOL
    "Do not have Algebra Plugins set up Eigen3 for itself"
)
set(ALGEBRA_PLUGINS_SETUP_VC
    OFF
    CACHE BOOL
    "Do not have Algebra Plugins set up Vc for itself"
)
set(ALGEBRA_PLUGINS_USE_SYSTEM_VC
    OFF
    CACHE BOOL
    "Have Algebra Plugins build Vc itself"
)
set(ALGEBRA_PLUGINS_SETUP_VECMEM
    FALSE
    CACHE BOOL
    "Do not set up VecMem in Algebra Plugins"
)
set(ALGEBRA_PLUGINS_SETUP_GOOGLETEST
    FALSE
    CACHE BOOL
    "Do not set up GoogleTest in Algebra Plugins"
)

set(ALGEBRA_PLUGINS_SETUP_BENCHMARK
    FALSE
    CACHE BOOL
    "Do not setup Algebra Plugins benchmark"
)

# Now set up its build.
FetchContent_MakeAvailable(AlgebraPlugins)
