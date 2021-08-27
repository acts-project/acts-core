// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <memory>

namespace ActsExamples {
class IBaseDetector;
}

using AlignedDetElementGetter =
    std::function<std::vector<Acts::DetectorElementBase*>(
        const std::shared_ptr<ActsExamples::IBaseDetector>&,
        const std::vector<Acts::GeometryIdentifier>&)>;

/// The options for running alignment or not
///
/// @param desc The options description to add options to
void addAlignmentOptions(ActsExamples::Options::Description& desc);

/// Main function for running alignment for specific detector.
///
/// @param argc number of command line arguments
/// @param argv command line arguments
/// @param detector is the detector to be aligned
int runDetectorAlignment(
    int argc, char* argv[],
    std::shared_ptr<ActsExamples::IBaseDetector> detector,
    ActsAlignment::AlignedTransformUpdater alignedTransformUpdater,
    AlignedDetElementGetter alignedDetElementsGetter);
