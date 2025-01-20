// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

#include <nlohmann/json.hpp>

/// Custom Json encoder/decoders. Naming is mandated by nlohmann::json and thus
/// can not match our naming guidelines.
///
/// This uses a custom API and nemclature as it would
/// otherwise require the ProtoAxis to have a default
/// constructor which is disallowed to to the unique_ptr
namespace Acts::ProtoAxisJsonConverter {
/// @brief Write the ProtoAxis to a json object
///
/// @param pa the proto axis to be written out
nlohmann::json toJson(const ProtoAxis& pa);

/// Create a ProtoAxis from a json object
///
/// @param j the json object to be read from
Acts::ProtoAxis fromJson(const nlohmann::json& j);

}  // namespace Acts::ProtoAxisJsonConverter
