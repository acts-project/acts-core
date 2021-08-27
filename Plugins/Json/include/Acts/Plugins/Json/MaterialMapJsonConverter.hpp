// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Plugins/Json/GeometryHierarchyMapJsonConverter.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <map>
#include <memory>

namespace Acts {

using SurfaceAndMaterialWithContext =
    std::tuple<std::shared_ptr<const Acts::Surface>,
               std::shared_ptr<const Acts::ISurfaceMaterial>,
               Acts::GeometryContext>;
using TrackingVolumeAndMaterial =
    std::pair<const Acts::TrackingVolume*,
              std::shared_ptr<const Acts::IVolumeMaterial>>;

/// @class MaterialMapJsonConverter
///
/// @brief read the material from Json
class MaterialMapJsonConverter {
 public:
  using SurfaceMaterialMap =
      std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>;
  using VolumeMaterialMap =
      std::map<GeometryIdentifier, std::shared_ptr<const IVolumeMaterial>>;
  using DetectorMaterialMaps = std::pair<SurfaceMaterialMap, VolumeMaterialMap>;

  /// @class Config
  /// Configuration of the Converter
  class Config {
   public:
    /// Default geometry context to extract surface tranforms
    GeometryContext context = GeometryContext();

    /// Steering to handle sensitive data
    bool processSensitives = true;
    /// Steering to handle approach data
    bool processApproaches = true;
    /// Steering to handle representing data
    bool processRepresenting = true;
    /// Steering to handle boundary data
    bool processBoundaries = true;
    /// Steering to handle volume data
    bool processVolumes = true;
    /// Steering to handle volume data
    bool processDenseVolumes = false;
    /// Add proto material to all surfaces
    bool processNonMaterial = false;
  };

  /// Constructor
  ///
  /// @param config configuration struct for the reader
  /// @param level The log level
  MaterialMapJsonConverter(const Config& config, Acts::Logging::Level level);

  /// Destructor
  ~MaterialMapJsonConverter() = default;

  /// Convert a json material map to a DetectorMaterialMaps
  ///
  /// @param materialmaps The json material
  DetectorMaterialMaps jsonToMaterialMaps(const nlohmann::json& materialmaps);

  /// Convert a DetectorMaterialMaps to json
  ///
  /// @param maps The material map collection
  nlohmann::json materialMapsToJson(const DetectorMaterialMaps& maps);

  /// Convert a tracking geometry to json.
  /// Can be used to initialise the material mapping process.
  ///
  /// @param tGeometry is the tracking geometry
  nlohmann::json trackingGeometryToJson(const TrackingGeometry& tGeometry);

  /// Go through a volume to find subvolume, layers and surfaces.
  /// Store volumes and surfaces in two vector used to initialised the geometry
  /// hierachy.
  ///
  /// @param volumeHierarchy is a vector of volume to be filled
  /// @param surfaceHierarchy is a vector of surfaces to be filled
  /// @param tVolume is a volume
  void convertToHierarchy(
      std::vector<std::pair<GeometryIdentifier,
                            Acts::TrackingVolumeAndMaterial>>& volumeHierarchy,
      std::vector<
          std::pair<GeometryIdentifier, Acts::SurfaceAndMaterialWithContext>>&
          surfaceHierarchy,
      const Acts::TrackingVolume* tVolume);

 private:
  /// The config class
  Config m_cfg;

  /// The logger instance
  std::unique_ptr<const Logger> m_logger{nullptr};

  /// Name of the volume hierarchy
  std::string m_volumeName = "Material Volume Map";
  /// Geometry hierarchy writer for volume material.
  Acts::GeometryHierarchyMapJsonConverter<const IVolumeMaterial*>
      m_volumeMaterialConverter;
  /// Geometry hierarchy writer for tracking volume.
  Acts::GeometryHierarchyMapJsonConverter<Acts::TrackingVolumeAndMaterial>
      m_volumeConverter;

  /// Name of the surface hierarchy
  std::string m_surfaceName = "Material Surface Map";
  /// Geometry hierarchy writer for surface material.
  Acts::GeometryHierarchyMapJsonConverter<const ISurfaceMaterial*>
      m_surfaceMaterialConverter;
  /// Geometry hierarchy writer for surface.
  Acts::GeometryHierarchyMapJsonConverter<Acts::SurfaceAndMaterialWithContext>
      m_surfaceConverter;

  /// Private access to the logging instance
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts
