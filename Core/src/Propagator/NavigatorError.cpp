// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/NavigatorError.hpp"

#include <string>

namespace {

class NavigatorErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category.
  const char* name() const noexcept final { return "NavigatorError"; }

  // Return what each enum means in text.
  std::string message(int c) const final {
    using Acts::NavigatorError;

    switch (static_cast<NavigatorError>(c)) {
      case NavigatorError::NotInsideExpectedVolume:
        return "We did not end up inside the volume.";
      case NavigatorError::NotOnExpectedSurface:
        return "Stepper not on surface";
      default:
        return "unknown";
    }
  }
};

}  // namespace

std::error_code Acts::make_error_code(Acts::NavigatorError e) {
  static NavigatorErrorCategory c;
  return {static_cast<int>(e), c};
}
