// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Visualization/ViewConfig.hpp>

#include <fstream>
#include <iostream>
#include <mutex>

namespace Acts {
class TrackingVolume;
class TrackingGeometry;
}  // namespace Acts

namespace ActsExamples {

/// @class ObjTrackingGeometryWriter
///
/// An Obj writer for the geometry: TrackingGeometry master
/// It delegates the writing of surfaces to the surface writers
class ObjTrackingGeometryWriter {
 public:
  // @class Config
  //
  // The nested config class
  class Config {
   public:
    double outputScalor = 1.0;   ///< scale output values
    size_t outputPrecision = 6;  ///< floating point precision
    std::string outputDir = ".";

    Acts::ViewConfig containerView = Acts::ViewConfig({220, 220, 220});
    Acts::ViewConfig volumeView = Acts::ViewConfig({220, 220, 0});
    Acts::ViewConfig sensitiveView = Acts::ViewConfig({0, 180, 240});
    Acts::ViewConfig passiveView = Acts::ViewConfig({240, 280, 0});
    Acts::ViewConfig gridView = Acts::ViewConfig({220, 0, 0});
  };

  /// Constructor
  /// @param config is the configuration class
  /// @param level the log level
  ObjTrackingGeometryWriter(const Config& config, Acts::Logging::Level level);

  /// Framework name() method
  /// @return the name of the tool
  std::string name() const;

  /// The write interface
  /// @param context the Algorithm/Event context of this call
  /// @param tGeometry is the geometry to be written out
  /// @return ProcessCode to indicate success/failure
  ActsExamples::ProcessCode write(const AlgorithmContext& context,
                                  const Acts::TrackingGeometry& tGeometry);

 private:
  std::unique_ptr<const Acts::Logger> m_logger;  ///< the logger instance

  Config m_cfg;  ///< the config class

  /// process this volume
  /// @param context the Algorithm/Event context for this call
  /// @param tVolume the volume to be processed
  void write(const AlgorithmContext& context,
             const Acts::TrackingVolume& tVolume);

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
