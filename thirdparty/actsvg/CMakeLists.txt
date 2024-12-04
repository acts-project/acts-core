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
message(STATUS "Building actsvg as part of the ACTS project")

FetchContent_Declare(actsvg SYSTEM ${ACTS_ACTSVG_SOURCE})

set(ACTSVG_BUILD_PYTHON_BINDINGS
    OFF
    CACHE BOOL
    "Do not build the Python bindings for actsvg"
)

# Now set up its build.
FetchContent_MakeAvailable(actsvg)
