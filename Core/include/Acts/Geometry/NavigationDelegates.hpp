// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/NavigationState.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <memory>
#include <tuple>

/// @note this is foreseen for the 'Geometry' module

namespace Acts {

class Surface;

namespace Experimental {

class Portal;
class DetectorVolume;
class Detector;

/// Base class for navigation delegates
class INavigationDelegate {
    public:
        INavigationDelegate() = default;

};

/// Declare an updator for the local navigation, i.e. the
/// navigation inside a detector volume. This can be called
/// either directly after a volume switch or in order to update
/// within a volume after some progression
///
/// This delegate dispatches the local navigation action
/// to a dedicated struct or function that is optimised for
/// the given environment.
///
/// @param gctx is the current geometry context
/// @param nState [in,out] is the navigation state to be updated
///
/// @note it relies on the detector volume to be set to the state
/// Memory  managed navigation state updator
using SurfaceCandidatesUpdator =
    OwningDelegate<void(const GeometryContext& gctx, NavigationState& nState),
                   INavigationDelegate>;

/// Declare a Detctor Volume finding or switching delegate
///
/// @param gctx is the current geometry context
/// @param nState [in, out] is the navigation state to be updated
///
/// @return the new DetectorVolume into which one changes at this switch
using DetectorVolumeUpdator =
    OwningDelegate<void(const GeometryContext& gctx, NavigationState& nState),
                   INavigationDelegate>;

// Construct and return a void delegate for initialization
inline static OwningDelegate<void(const GeometryContext& gctx,
                           NavigationState& nState),
                      INavigationDelegate>
voidDelegate() {
  OwningDelegate<void(const GeometryContext& gctx, NavigationState& nState),
                 INavigationDelegate>
      d;
  return d;
}

}  // namespace Experimental
}  // namespace Acts
