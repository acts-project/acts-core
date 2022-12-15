// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <ostream>
#include <iostream>

std::ostream& Acts::operator<<(std::ostream& os, Acts::GeometryIdentifier id) {
  // zero represents an invalid/undefined identifier
  if (id.value() == 0u) {
    return (os << "undefined");
  }

  static const char* const names[] = {
      "vol=", "bnd=", "lay=", "apr=", "sen=", "ext=",
  };
  const GeometryIdentifier::Value levels[] = {id.volume(),    id.boundary(),
                                              id.layer(),     id.approach(),
                                              id.sensitive(), id.extra()};

  bool writeSeparator = false;
  for (auto i = 0u; i < (sizeof(levels) / sizeof(levels[0])); ++i) {
    if (levels[i] != 0u) {
      if (writeSeparator) {
        os << '|';
      }
      os << names[i] << levels[i];
      writeSeparator = true;
    }
  }
  return os;
}

Acts::GeometryIdentifier Acts::GeometryIdentifierHook::decorateIdentifier(
    Acts::GeometryIdentifier identifier, const Acts::Surface&) const {
  std::cout << "NO-HOOK\n";
  return identifier;
}
