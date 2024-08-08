// This file is part of the Acts project.
//
// Copyright (C) 2019-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DetectorCommons/Detector.hpp"

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

namespace ActsExamples::Generic {
class GenericDetectorElement;

class GenericDetector : public ActsExamples::DetectorCommons::Detector {
 public:
  using TrackingGeometryPtr = std::shared_ptr<const Acts::TrackingGeometry>;
  using ContextDecorators =
      std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>;

  using DetectorElement = ActsExamples::Generic::GenericDetectorElement;

  struct Config {
    std::size_t buildLevel{3};
    Acts::Logging::Level logLevel{Acts::Logging::INFO};
    Acts::Logging::Level surfaceLogLevel{Acts::Logging::INFO};
    Acts::Logging::Level layerLogLevel{Acts::Logging::INFO};
    Acts::Logging::Level volumeLogLevel{Acts::Logging::INFO};
    bool buildProto{false};
    std::shared_ptr<const Acts::IMaterialDecorator> materialDecorator;
  };

  explicit GenericDetector(const Config& cfg);

 private:
  Config m_cfg;

  void buildTrackingGeometry() final;
};

}  // namespace ActsExamples::Generic
