// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Plugins/Obj/ObjTrackingGeometryWriter.hpp"

#include <Acts/Geometry/Layer.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Visualization/GeometryView3D.hpp>
#include <Acts/Visualization/ObjVisualization3D.hpp>

#include <iostream>

ActsExamples::ObjTrackingGeometryWriter::ObjTrackingGeometryWriter(
    const ActsExamples::ObjTrackingGeometryWriter::Config& config,
    Acts::Logging::Level level)
    : m_logger{Acts::getDefaultLogger(name(), level)}, m_cfg(config) {}

std::string ActsExamples::ObjTrackingGeometryWriter::name() const {
  return "ObjTrackingGeometryWriter";
}

ActsExamples::ProcessCode ActsExamples::ObjTrackingGeometryWriter::write(
    const AlgorithmContext& context, const Acts::TrackingGeometry& tGeometry) {
  ACTS_DEBUG(">>Obj: Writer for TrackingGeometry object called.");

  auto world = tGeometry.highestTrackingVolume();
  if (world) {
    write(context, *world);
  }
  return ActsExamples::ProcessCode::SUCCESS;
}

void ActsExamples::ObjTrackingGeometryWriter::write(
    const AlgorithmContext& context, const Acts::TrackingVolume& tVolume) {
  ACTS_DEBUG(">>Obj: Writer for TrackingVolume object called.");

  Acts::ObjVisualization3D objVis(m_cfg.outputPrecision, m_cfg.outputScalor);

  Acts::GeometryView3D::drawTrackingVolume(
      objVis, tVolume, context.geoContext, m_cfg.containerView,
      m_cfg.volumeView, m_cfg.passiveView, m_cfg.sensitiveView, m_cfg.gridView,
      true, "", m_cfg.outputDir);
}
