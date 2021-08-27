// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/VolumeJsonConverter.hpp"

#include "Acts/Plugins/Json/GeometryJsonKeys.hpp"
#include "Acts/Plugins/Json/MaterialJsonConverter.hpp"

void Acts::to_json(nlohmann::json& j,
                   const Acts::TrackingVolumeAndMaterial& volume) {
  j[Acts::jsonKey().namekey] = volume.first->volumeName();
  to_json(j, volume.second.get());
}

void Acts::to_json(nlohmann::json& j,
                   const Acts::TrackingVolumePointer& volume) {
  to_json(j, *volume);
}

void Acts::to_json(nlohmann::json& j, const Acts::TrackingVolume& volume) {
  j[Acts::jsonKey().namekey] = volume.volumeName();
  if (volume.volumeMaterial() != nullptr) {
    to_json(j, volume.volumeMaterial());
  }
  return;
}
