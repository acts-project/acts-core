// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/Detector.hpp"
#include "Acts/Plugins/Detray/DetrayConversionUtils.hpp"
#include "Acts/Plugins/Detray/DetrayGeometryConverter.hpp"
#include "Acts/Plugins/Detray/DetrayMaterialConverter.hpp"
#include "Acts/Plugins/Detray/DetraySurfaceGridsConverter.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

#include <detray/io/common/geometry_reader.hpp>
#include <detray/io/common/material_map_reader.hpp>
#include <detray/io/common/surface_grid_reader.hpp>
#include <detray/io/frontend/detector_writer_config.hpp>
#include <detray/utils/consistency_checker.hpp>

namespace Acts {

class DetrayConverter {
 public:
  /// Detray conversion options
  struct Options {
    /// Option to switch on/off the material conversion
    bool convertMaterial = true;
    /// Option to switch on/off the surface grid conversin
    bool convertSurfaceGrids = true;
    /// Option to switch on/off the export to json
    bool writeToJson = false;
  };

  /// Constructor with logger
  explicit DetrayConverter(
      std::unique_ptr<const Logger> logger = getDefaultLogger("DetrayConverter",
                                                              Logging::INFO));

  /// Convert an Acts::Experimental::Detector to a detray::detector object
  ///
  /// @param gctx the geometry context
  /// @param detector the detector to be converted
  /// @param mr the memory resource to be used
  /// @param options the conversion options
  ///
  /// @returns a detector of requested return type
  template <typename detector_t = DetrayHostDetector>
  detector_t convert(const GeometryContext& gctx,
                     const Experimental::Detector& detector,
                     vecmem::memory_resource& mr, const Options& options) {
    // The building cache object
    DetrayConversionUtils::Cache cCache(detector.volumes());

    typename detector_t::name_map names = {{0u, detector.name()}};

    // build detector
    detray::detector_builder<typename detector_t::metadata> detectorBuilder{};
    // (1) geometry
    detray::io::detector_payload detectorPayload =
        DetrayGeometryConverter::convertDetector(cCache, gctx, detector,
                                                 logger());
    detray::io::geometry_reader::convert<detector_t>(detectorBuilder, names,
                                                     detectorPayload);

    // (2a) homogeneous material
    if constexpr (detray::concepts::has_homogeneous_material<detector_t>) {
      if (options.convertMaterial) {
        detray::io::detector_homogeneous_material_payload materialSlabsPayload =
            DetrayMaterialConverter::convertHomogeneousSurfaceMaterial(
                cCache, detector, logger());
        detray::io::homogeneous_material_reader::convert<detector_t>(
            detectorBuilder, names, std::move(materialSlabsPayload));
      }
    }

    // (2b) material grids
    if constexpr (detray::concepts::has_material_rods<detector_t>) {
      if (options.convertMaterial) {
        detray::io::detector_grids_payload<detray::io::material_slab_payload,
                                           detray::io::material_id>
            materialGridsPayload =
                DetrayMaterialConverter::convertGridSurfaceMaterial(
                    cCache, detector, logger());
        detray::io::material_map_reader<std::integral_constant<
            std::size_t, 2>>::convert<detector_t>(detectorBuilder, names,
                                                  std::move(
                                                      materialGridsPayload));
      }
    }

    // (3) surface grids
    if (options.convertSurfaceGrids) {
      detray::io::detector_grids_payload<std::size_t, detray::io::accel_id>
          surfaceGridsPayload =
              DetraySurfaceGridsConverter::convertSurfaceGrids(detector);

      // Capacity 0 (dynamic bin size) and dimension 2 (2D grids)
      detray::io::surface_grid_reader<typename detector_t::surface_type,
                                      std::integral_constant<std::size_t, 0>,
                                      std::integral_constant<std::size_t, 2>>::
          template convert<detector_t>(detectorBuilder, names,
                                       std::move(surfaceGridsPayload));
    }

    detector_t detrayDetector(detectorBuilder.build(mr));

    // Checks and print
    detray::detail::check_consistency(detrayDetector);

    // If configured, write the detector to json
    if (options.writeToJson) {
      // Create a writer configuration and write it out
      detray::io::detector_writer_config writerConfig{};
      writerConfig.m_write_material = options.convertMaterial;
      writerConfig.m_write_grids = options.convertSurfaceGrids;
      writeToJson(detrayDetector, names, writerConfig);
    }

    return detrayDetector;
  }

  /// Write the detector to json output
  ///
  /// @param dDetector is the detray detector (converted)
  /// @param names a name map for the detector volumes
  /// @param writer_cfg the writer configuration
  static void writeToJson(
      const DetrayHostDetector& dDetector,
      const typename DetrayHostDetector::name_map& names = {},
      detray::io::detector_writer_config writer_cfg = {});

 private:
  /// The logger instance
  std::unique_ptr<const Logger> m_logger = nullptr;

  // Return the logging instance
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace Acts
