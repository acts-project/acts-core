# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 3.25)

include(FetchContent)

message(STATUS "Building FRNN as part of the ACTS project")

# Declare where to get frnncontent from
FetchContent_Declare(frnncontent SYSTEM ${ACTS_FRNN_SOURCE})

# FRNN does not provide a CMakeLists.txt, so we use a custom one. Because of this,
# we have to implement the populate step manually
FetchContent_GetProperties(frnncontent)
if(NOT frnncontent_POPULATED)
    FetchContent_Populate(frnncontent)
    configure_file(
        CMakeLists.txt.in
        "${frnncontent_SOURCE_DIR}/CMakeLists.txt"
        COPYONLY
    )
    add_subdirectory(${frnncontent_SOURCE_DIR} ${frnncontent_BINARY_DIR})
endif()
