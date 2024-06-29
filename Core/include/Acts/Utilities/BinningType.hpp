// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <string>
#include <type_traits>
#include <vector>

namespace Acts {

/// @enum BinningType, BinningOption & BinningAccess
///
///- BinningType:
///
///  Enumeration to qualify the binning type for the use of the
///  LayerArrayCreator and the TrackingVolumeArrayCreator
///
/// - BinningOption:
///   open:   [0,max]
///   closed:  0 -> nextbin -> max -> 0
///
/// - BinningValue
///   necessary access to global positions
///
enum BinningType { equidistant, arbitrary };

/// @brief flag for open/closed bins
enum BinningOption { open, closed };

/// @enum BinningValue how to take the global / local position
enum class BinningValue : int {
  binX = 0,
  binY = 1,
  binZ = 2,
  binR = 3,
  binPhi = 4,
  binRPhi = 5,
  binH = 6,
  binEta = 7,
  binMag = 8,
  binValues = 9
};

/// @brief static list of all binning values
static const std::vector<BinningValue> s_binningValues = {
    BinningValue::binX, BinningValue::binY,   BinningValue::binZ,
    BinningValue::binR, BinningValue::binPhi, BinningValue::binRPhi,
    BinningValue::binH, BinningValue::binEta, BinningValue::binMag};

inline const std::vector<std::string>& binningValueNames() {
  static const std::vector<std::string> _binningValueNames = {
      "binX",    "binY", "binZ",   "binR",  "binPhi",
      "binRPhi", "binH", "binEta", "binMag"};
  return _binningValueNames;
}

/// @brief screen output option
inline const std::string& binningValueName(BinningValue bValue) {
  return binningValueNames()[static_cast<std::underlying_type_t<BinningValue>>(
      bValue)];
}

inline std::ostream& operator<<(std::ostream& os, BinningValue bValue) {
  os << binningValueName(bValue);
  return os;
}

}  // namespace Acts
