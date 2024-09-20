// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/NavigationDelegate.hpp"

#include "Acts/Geometry/TrackingVolume.hpp"

namespace Acts {

void TryAllPortalNavigationDelegate::updateState(
    Experimental::Gen3Geometry::NavigationState& state) const {
  assert(m_volume != nullptr);
  assert(state.currentVolume == m_volume);

  for (const auto& portal : state.currentVolume->portals()) {
    state.main.addPortalCandidate(portal);
  };
}
}  // namespace Acts
