// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <string>        // for string printing
#include <system_error>  // bring in std::error_code et al

namespace ActsAlignment {
// This is the custom error code enum
enum class AlignmentError {
  NoAlignmentDofOnTrack = 1,
  AlignmentParametersUpdateFailure = 2,
  ConvergeFailure = 3
};

namespace detail {
// Define a custom error code category derived from std::error_category
class AlignmentErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category
  const char* name() const noexcept final { return "AlignmentError"; }
  // Return what each enum means in text
  std::string message(int c) const final {
    switch (static_cast<AlignmentError>(c)) {
      case AlignmentError::NoAlignmentDofOnTrack:
        return "No alignment parameters on the track";
      case AlignmentError::AlignmentParametersUpdateFailure:
        return "Update to alignment parameters failure";
      case AlignmentError::ConvergeFailure:
        return "The alignment is not converged";
      default:
        return "unknown";
    }
  }
};
}  // namespace detail

// Declare a global function returning a static instance of the custom category
extern inline const detail::AlignmentErrorCategory& AlignmentErrorCategory() {
  static detail::AlignmentErrorCategory c;
  return c;
}

inline std::error_code make_error_code(ActsAlignment::AlignmentError e) {
  return {static_cast<int>(e), ActsAlignment::AlignmentErrorCategory()};
}
}  // namespace ActsAlignment

namespace std {
// register with STL
template <>
struct is_error_code_enum<ActsAlignment::AlignmentError> : std::true_type {};
}  // namespace std
