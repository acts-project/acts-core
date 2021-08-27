// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsAlignment/Kernel/detail/AlignmentEngine.hpp"

namespace ActsAlignment {
namespace detail {

void resetAlignmentDerivative(Acts::AlignmentToBoundMatrix& alignToBound,
                              AlignmentMask mask) {
  if (not ACTS_CHECK_BIT(mask, AlignmentMask::Center0)) {
    alignToBound.col(Acts::eAlignmentCenter0) = Acts::AlignmentVector::Zero();
  }
  if (not ACTS_CHECK_BIT(mask, AlignmentMask::Center1)) {
    alignToBound.col(Acts::eAlignmentCenter1) = Acts::AlignmentVector::Zero();
  }
  if (not ACTS_CHECK_BIT(mask, AlignmentMask::Center2)) {
    alignToBound.col(Acts::eAlignmentCenter2) = Acts::AlignmentVector::Zero();
  }
  if (not ACTS_CHECK_BIT(mask, AlignmentMask::Rotation0)) {
    alignToBound.col(Acts::eAlignmentRotation0) = Acts::AlignmentVector::Zero();
  }
  if (not ACTS_CHECK_BIT(mask, AlignmentMask::Rotation1)) {
    alignToBound.col(Acts::eAlignmentRotation1) = Acts::AlignmentVector::Zero();
  }
  if (not ACTS_CHECK_BIT(mask, AlignmentMask::Rotation2)) {
    alignToBound.col(Acts::eAlignmentRotation2) = Acts::AlignmentVector::Zero();
  }
}

}  // namespace detail
}  // namespace ActsAlignment
