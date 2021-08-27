// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <memory>
#include <vector>

namespace Acts {
class TGeoDetectorElement;
}

namespace ActsExamples {

struct TGeoDetector : public ActsExamples::IBaseDetector {
  using DetectorElementPtr = std::shared_ptr<const Acts::TGeoDetectorElement>;
  using DetectorStore = std::vector<DetectorElementPtr>;

  /// The Store of the detector elements (lifetime: job)
  DetectorStore detectorStore;

  struct Config {
    Acts::Logging::Level surfaceLogLevel = Acts::Logging::WARNING;
    Acts::Logging::Level layerLogLevel = Acts::Logging::WARNING;
    Acts::Logging::Level volumeLogLevel = Acts::Logging::WARNING;

    std::string fileName;
    bool buildBeamPipe = false;
    double beamPipeRadius;
    double beamPipeHalflengthZ;
    double beamPipeLayerThickness;

    double unitScalor = 1.0;

    enum SubVolume : size_t { Negative = 0, Central, Positive };

    template <typename T>
    struct LayerTriplet {
      LayerTriplet() = default;
      LayerTriplet(T value)
          : negative{value}, central{value}, positive{value} {}

      LayerTriplet(T _negative, T _central, T _positive)
          : negative{_negative}, central{_central}, positive{_positive} {}

      T negative;
      T central;
      T positive;

      T& at(SubVolume i) {
        switch (i) {
          case Negative:
            return negative;
          case Central:
            return central;
          case Positive:
            return positive;
          default:
            throw std::invalid_argument{"Unknown index"};
        }
      }

      const T& at(SubVolume i) const {
        switch (i) {
          case Negative:
            return negative;
          case Central:
            return central;
          case Positive:
            return positive;
          default:
            throw std::invalid_argument{"Unknown index"};
        }
      }
    };

    struct Volume {
      std::string name;
      LayerTriplet<bool> layers{false};
      LayerTriplet<std::string> subVolumeName;
      LayerTriplet<std::vector<std::string>> sensitiveNames;
      LayerTriplet<std::string> sensitiveAxes;
      LayerTriplet<Options::Interval> rRange;
      LayerTriplet<Options::Interval> zRange;
      LayerTriplet<double> splitTolR{0};
      LayerTriplet<double> splitTolZ{0};

      Options::Interval binToleranceR;
      Options::Interval binTolerancePhi;
      Options::Interval binToleranceZ;

      bool cylinderDiscSplit = false;
      unsigned int cylinderNZSegments;
      unsigned int cylinderNPhiSegments;
      unsigned int discNRSegments;
      unsigned int discNPhiSegments;
    };

    std::vector<Volume> volumes;
  };

  void addOptions(
      boost::program_options::options_description& opt) const override;

  std::pair<ActsExamples::IBaseDetector::TrackingGeometryPtr, ContextDecorators>
  finalize(const boost::program_options::variables_map& vm,
           std::shared_ptr<const Acts::IMaterialDecorator> mdecorator) override;

  std::pair<ActsExamples::IBaseDetector::TrackingGeometryPtr, ContextDecorators>
  finalize(const Config& cfg,
           std::shared_ptr<const Acts::IMaterialDecorator> mdecorator);
};

}  // namespace ActsExamples