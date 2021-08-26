// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Plugins/Json/MaterialMapJsonConverter.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <fstream>
#include <map>
#include <mutex>

// Convenience shorthand

namespace Acts {

/// @brief Material decorator from Json format
///
/// This reads in material maps for surfaces and volumes
/// from a json file
class JsonMaterialDecorator : public IMaterialDecorator {
 public:
  using SurfaceMaterialMap =
      std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>;

  using VolumeMaterialMap =
      std::map<GeometryIdentifier, std::shared_ptr<const IVolumeMaterial>>;

  JsonMaterialDecorator(const MaterialMapJsonConverter::Config& rConfig,
                        const std::string& jFileName,
                        Acts::Logging::Level level,
                        bool clearSurfaceMaterial = true,
                        bool clearVolumeMaterial = true)
      : m_readerConfig(rConfig),
        m_clearSurfaceMaterial(clearSurfaceMaterial),
        m_clearVolumeMaterial(clearVolumeMaterial),
        m_logger{getDefaultLogger("JsonMaterialDecorator", level)} {
    // the material reader
    Acts::MaterialMapJsonConverter jmConverter(rConfig, level);

    ACTS_VERBOSE("Reading JSON material description from: " << jFileName);
    std::ifstream ifj(jFileName.c_str());
    if (!ifj.good()) {
      throw std::runtime_error{"Unable to open input JSON material file"};
    }
    nlohmann::json jin;

    if (jFileName.find(".cbor") != std::string::npos) {
      std::vector<std::uint8_t> iCbor((std::istreambuf_iterator<char>(ifj)),
                                      std::istreambuf_iterator<char>());
      jin = nlohmann::json::from_cbor(iCbor);
    } else {
      ifj >> jin;
    }

    auto maps = jmConverter.jsonToMaterialMaps(jin);
    m_surfaceMaterialMap = maps.first;
    m_volumeMaterialMap = maps.second;
    ACTS_VERBOSE("JSON material description read complete");
  }

  /// Decorate a surface
  ///
  /// @param surface the non-cost surface that is decorated
  void decorate(Surface& surface) const final {
    ACTS_VERBOSE("Processing surface: " << surface.geometryId());
    // Clear the material if registered to do so
    if (m_clearSurfaceMaterial) {
      ACTS_VERBOSE("-> Clearing surface material");
      surface.assignSurfaceMaterial(nullptr);
    }
    // Try to find the surface in the map
    auto sMaterial = m_surfaceMaterialMap.find(surface.geometryId());
    if (sMaterial != m_surfaceMaterialMap.end()) {
      ACTS_VERBOSE("-> Found material for surface, assigning");
      surface.assignSurfaceMaterial(sMaterial->second);
    }
  }

  /// Decorate a TrackingVolume
  ///
  /// @param volume the non-cost volume that is decorated
  void decorate(TrackingVolume& volume) const final {
    ACTS_VERBOSE("Processing volume: " << volume.geometryId());
    // Clear the material if registered to do so
    if (m_clearVolumeMaterial) {
      ACTS_VERBOSE("-> Clearing volume material");
      volume.assignVolumeMaterial(nullptr);
    }
    // Try to find the volume in the map
    auto vMaterial = m_volumeMaterialMap.find(volume.geometryId());
    if (vMaterial != m_volumeMaterialMap.end()) {
      ACTS_VERBOSE("-> Found material for volume, assigning");
      volume.assignVolumeMaterial(vMaterial->second);
    }
  }

 private:
  MaterialMapJsonConverter::Config m_readerConfig;
  SurfaceMaterialMap m_surfaceMaterialMap;
  VolumeMaterialMap m_volumeMaterialMap;

  bool m_clearSurfaceMaterial{true};
  bool m_clearVolumeMaterial{true};

  std::unique_ptr<const Logger> m_logger;

  const Logger& logger() const { return *m_logger; }
};
}  // namespace Acts
