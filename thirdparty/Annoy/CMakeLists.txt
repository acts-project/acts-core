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

# Information message.
message(STATUS "Building Annoy as part of the ACTS project")

# Version of the Annoy library to fetch.
set(ANNOY_VERSION "v${_acts_annoy_version}")

# Declare the fetching of the Annoy library.
FetchContent_Declare(
    Annoy
    SYSTEM
    ${ACTS_ANNOY_SOURCE}
    PATCH_COMMAND
        patch -p1 < ${CMAKE_CURRENT_SOURCE_DIR}/0001-Modify-annoy.patch
)

# Make the fetched content available.
# Annoy contains a CMakeLists.txt file, so it will be added as a subdirectory.
FetchContent_MakeAvailable(Annoy)
