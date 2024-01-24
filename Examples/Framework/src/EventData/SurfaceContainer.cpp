// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/SurfaceContainer.hpp"

Acts::SurfacePtrsContainer Acts::SurfaceContainer::getPtrs(
    DetectorPtr detector) const {
  Acts::SurfaceContainer::SurfaceVisitor visitor;
  for (auto& vol : detector->rootVolumePtrs()) {
    for (auto& surf : vol->surfacePtrs()) {
      visitor(&(*surf));
    }
  }
  return visitor.surfacePtrs;
}

Acts::SurfacePtrsContainer Acts::SurfaceContainer::getPtrs(
    TrackingGeometryPtr tGeometryPtr) const {
  Acts::SurfaceContainer::SurfaceVisitor visitor;
  tGeometryPtr->visitSurfaces(visitor);
  return visitor.surfacePtrs;
}