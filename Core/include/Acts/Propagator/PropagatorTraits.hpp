// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <type_traits>

namespace Acts {

template <typename stepper_t>
struct SupportsBoundParameters : public std::false_type {};

template <typename stepper_t>
constexpr bool SupportsBoundParameters_v =
    SupportsBoundParameters<stepper_t>::value;

}  // namespace Acts
