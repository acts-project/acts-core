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
message(STATUS "Building Detray as part of the Acts project")

set(DETRAY_VERSION "v${_acts_detray_version}")

# Declare where to get Detray from.
FetchContent_Declare(Detray SYSTEM ${ACTS_DETRAY_SOURCE})

# Options used in the build of Detray.
set(ACTS_DETRAY_SCALARTYPE "double")

set(DETRAY_CUSTOM_SCALARTYPE
    "${ACTS_DETRAY_SCALARTYPE}"
    CACHE STRING
    "Scalar type to use in the Detray code"
)

set(DETRAY_BUILD_TEST_UTILS
    ON
    CACHE BOOL
    "Turn on the build of the Detray test utilities"
)

set(DETRAY_BUILD_UNITTESTS
    OFF
    CACHE BOOL
    "Turn off the build of the Detray unit tests"
)
set(DETRAY_BUILD_INTEGRATIONTESTS
    OFF
    CACHE BOOL
    "Turn off the build of the Detray integration tests"
)
set(DETRAY_BUILD_TUTORIALS
    OFF
    CACHE BOOL
    "Turn off the build of the Detray tutorials"
)
set(DETRAY_BUILD_BENCHMARKS
    OFF
    CACHE BOOL
    "Turn off the build of the Detray benchmarks"
)
set(DETRAY_EIGEN_PLUGIN
    ON
    CACHE BOOL
    "Turn on the build of the Detray Eigen code"
)
set(DETRAY_VC_PLUGIN OFF CACHE BOOL "Turn off the build of the Detray Vc code")

set(DETRAY_SETUP_VECMEM OFF CACHE BOOL "Do not set up VecMem as part of Detray")
set(DETRAY_SETUP_ALGEBRA_PLUGINS
    OFF
    CACHE BOOL
    "Do not set up Algebra Plugins as part of Detray"
)
set(DETRAY_SETUP_GOOGLETEST
    OFF
    CACHE BOOL
    "Do not set up GoogleTest as part of Detray"
)
set(DETRAY_SETUP_BENCHMARK
    OFF
    CACHE BOOL
    "Do not set up Google Benchmark as part of Detray"
)
set(DETRAY_SETUP_THRUST ON CACHE BOOL "Set up Thrust as part of Detray")
set(DETRAY_SETUP_COVFIE OFF CACHE BOOL "Do not set up covfie as part of Detray")
set(DETRAY_SETUP_NLOHMANN
    OFF
    CACHE BOOL
    "Do not set up Nlohmann as part of Detray"
)
set(DETRAY_SETUP_ACTSVG OFF CACHE BOOL "Do not set up Actsvg as part of Detray")
set(DETRAY_SETUP_DFELIBS
    ON
    CACHE BOOL
    "Do not set up Dfelibs as part of Detray"
)
set(DETRAY_SVG_DISPLAY OFF CACHE BOOL "No not build the ActSVG display module")

#Now set up its build.
FetchContent_MakeAvailable(Detray)
