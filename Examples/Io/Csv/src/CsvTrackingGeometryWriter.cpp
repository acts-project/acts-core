// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvTrackingGeometryWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Plugins/Digitization/CartesianSegmentation.hpp"
#include "Acts/Plugins/Digitization/DigitizationModule.hpp"
#include "Acts/Plugins/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <dfe/dfe_io_dsv.hpp>

#include "CsvOutputData.hpp"

using namespace ActsExamples;

CsvTrackingGeometryWriter::CsvTrackingGeometryWriter(
    const CsvTrackingGeometryWriter::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_world(nullptr),
      m_logger(Acts::getDefaultLogger("CsvTrackingGeometryWriter", lvl))

{
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  m_world = m_cfg.trackingGeometry->highestTrackingVolume();
  if (not m_world) {
    throw std::invalid_argument("Could not identify the world volume");
  }
}

std::string CsvTrackingGeometryWriter::name() const {
  return "CsvTrackingGeometryWriter";
}

namespace {

using SurfaceWriter = dfe::NamedTupleCsvWriter<SurfaceData>;
using BoundarySurface = Acts::BoundarySurfaceT<Acts::TrackingVolume>;

/// Write a single surface.
void fillSurfaceData(SurfaceData& data, const Acts::Surface& surface,
                     const Acts::GeometryContext& geoCtx) {
  // encoded and partially decoded geometry identifier
  data.geometry_id = surface.geometryId().value();
  data.volume_id = surface.geometryId().volume();
  data.boundary_id = surface.geometryId().boundary();
  data.layer_id = surface.geometryId().layer();
  data.module_id = surface.geometryId().sensitive();
  // center position
  auto center = surface.center(geoCtx);
  data.cx = center.x() / Acts::UnitConstants::mm;
  data.cy = center.y() / Acts::UnitConstants::mm;
  data.cz = center.z() / Acts::UnitConstants::mm;
  // rotation matrix components are unit-less
  auto transform = surface.transform(geoCtx);
  data.rot_xu = transform(0, 0);
  data.rot_xv = transform(0, 1);
  data.rot_xw = transform(0, 2);
  data.rot_yu = transform(1, 0);
  data.rot_yv = transform(1, 1);
  data.rot_yw = transform(1, 2);
  data.rot_zu = transform(2, 0);
  data.rot_zv = transform(2, 1);
  data.rot_zw = transform(2, 2);

  std::array<float*, 7> dataBoundParameters = {
      &data.p0, &data.p1, &data.p2, &data.p3, &data.p4, &data.p5, &data.p6};

  const auto& bounds = surface.bounds();
  data.bounds_type = static_cast<int>(bounds.type());
  auto boundValues = bounds.values();
  for (size_t ipar = 0; ipar < boundValues.size(); ++ipar) {
    (*dataBoundParameters[ipar]) = boundValues[ipar];
  }
}

/// Write a single surface.
void writeSurface(SurfaceWriter& writer, const Acts::Surface& surface,
                  const Acts::GeometryContext& geoCtx) {
  SurfaceData data;
  fillSurfaceData(data, surface, geoCtx);
  writer.append(data);
}

/// Write a single surface.
void writeBoundarySurface(SurfaceWriter& writer,
                          const BoundarySurface& bsurface,
                          const Acts::GeometryContext& geoCtx) {
  SurfaceData data;
  fillSurfaceData(data, bsurface.surfaceRepresentation(), geoCtx);
  writer.append(data);
}

/// Write all child surfaces and descend into confined volumes.
void writeVolume(SurfaceWriter& writer, const Acts::TrackingVolume& volume,
                 bool writeSensitive, bool writeBoundary,
                 const Acts::GeometryContext& geoCtx) {
  // process all layers that are directly stored within this volume
  if (volume.confinedLayers()) {
    for (auto layer : volume.confinedLayers()->arrayObjects()) {
      // We jump navigation layers
      if (layer->layerType() == Acts::navigation) {
        continue;
      }
      // check for sensitive surfaces
      if (layer->surfaceArray() and writeSensitive) {
        for (auto surface : layer->surfaceArray()->surfaces()) {
          if (surface) {
            writeSurface(writer, *surface, geoCtx);
          }
        }
      }
    }
    // This is a navigation volume, write the boundaries
    if (writeBoundary) {
      for (auto bsurface : volume.boundarySurfaces()) {
        writeBoundarySurface(writer, *bsurface, geoCtx);
      }
    }
  }
  // step down into hierarchy to process all child volumnes
  if (volume.confinedVolumes()) {
    for (auto confined : volume.confinedVolumes()->arrayObjects()) {
      writeVolume(writer, *confined.get(), writeSensitive, writeBoundary,
                  geoCtx);
    }
  }
}
}  // namespace

ProcessCode CsvTrackingGeometryWriter::write(const AlgorithmContext& ctx) {
  if (not m_cfg.writePerEvent) {
    return ProcessCode::SUCCESS;
  }
  SurfaceWriter writer(
      perEventFilepath(m_cfg.outputDir, "detectors.csv", ctx.eventNumber),
      m_cfg.outputPrecision);
  writeVolume(writer, *m_world, m_cfg.writeSensitive, m_cfg.writeSensitive,
              ctx.geoContext);
  return ProcessCode::SUCCESS;
}

ProcessCode CsvTrackingGeometryWriter::endRun() {
  SurfaceWriter writer(joinPaths(m_cfg.outputDir, "detectors.csv"),
                       m_cfg.outputPrecision);
  writeVolume(writer, *m_world, m_cfg.writeSensitive, m_cfg.writeSensitive,
              Acts::GeometryContext());
  return ProcessCode::SUCCESS;
}
