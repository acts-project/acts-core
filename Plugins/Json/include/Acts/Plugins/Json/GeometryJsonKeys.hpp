// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <string_view>

namespace Acts::JsonKey {

/// The name identification
static constexpr std::string_view kName = "NAME";
/// The bin key
static constexpr std::string_view kBin = "binUtility";
/// The material key
static constexpr std::string_view kMaterial = "material";
/// The local to global transformation key
static constexpr std::string_view kTransform = "transformation";
/// The type key -> proto, else
static constexpr std::string_view kType = "type";
/// The data key
static constexpr std::string_view kData = "data";
/// The mapping key, add surface to mapping procedure if true
static constexpr std::string_view kMaterialMap = "mapMaterial";
/// The mapping type key, used to select the type od material mapping for the
/// surface
static constexpr std::string_view kMappingType = "mappingType";

}  // namespace Acts::JsonKey
