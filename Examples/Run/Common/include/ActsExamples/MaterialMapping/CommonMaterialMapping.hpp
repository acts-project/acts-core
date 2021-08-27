// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace ActsExamples {
class IBaseDetector;
}

/// @brief The material validation example, it runs a propagation
/// and then writes out the material information
///
/// @param argc the number of argumetns of the call
/// @param atgv the argument list
/// @param detector the detector instance
///
int materialMappingExample(int argc, char* argv[],
                           ActsExamples::IBaseDetector& detector);
