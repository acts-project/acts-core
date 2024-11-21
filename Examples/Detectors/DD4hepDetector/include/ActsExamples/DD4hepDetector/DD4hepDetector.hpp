// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DetectorCommons/DetectorBase.hpp"
#include "ActsExamples/Geant4/RegionCreator.hpp"

#include <functional>
#include <memory>
#include <string>
#include <vector>

class TGeoNode;

namespace dd4hep {
class Detector;
class DetElement;
}  // namespace dd4hep

namespace ActsExamples {

void sortFCChhDetElements(std::vector<dd4hep::DetElement>& det);

/// @class DD4hepDetector
///
/// @brief geometries from dd4hep input
///
/// The DD4hepDetector creates the DD4hep, the TGeo and the ACTS
/// TrackingGeometry from DD4hep xml input. The geometries are created only on
/// demand.
class DD4hepDetector : public DetectorBase,
                       public std::enable_shared_from_this<DD4hepDetector> {
 public:
  /// @brief The context decorators
  using ContextDecorators =
      std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>;

  struct Config {
    /// Log level for the geometry service.
    Acts::Logging::Level logLevel = Acts::Logging::Level::INFO;
    /// Log level for DD4hep itself
    Acts::Logging::Level dd4hepLogLevel = Acts::Logging::Level::WARNING;
    /// XML-file with the detector description
    std::vector<std::string> xmlFileNames;
    /// The name of the service
    std::string name = "default";
    /// Binningtype in phi
    Acts::BinningType bTypePhi = Acts::equidistant;
    /// Binningtype in r
    Acts::BinningType bTypeR = Acts::arbitrary;
    /// Binningtype in z
    Acts::BinningType bTypeZ = Acts::equidistant;
    /// The tolerance added to the geometrical extension in r
    /// of the layers contained to build the volume envelope around
    /// @note this parameter only needs to be set if the volumes containing
    /// the
    /// layers (e.g. barrel, endcap volumes) have no specific shape
    /// (assemblies)
    double envelopeR = 1 * Acts::UnitConstants::mm;
    /// The tolerance added to the geometrical extension in z
    /// of the layers contained to build the volume envelope around
    /// @note this parameter only needs to be set if the volumes containing
    /// the layers (e.g. barrel, endcap volumes) have no specific shape
    /// (assemblies)
    double envelopeZ = 1 * Acts::UnitConstants::mm;
    double defaultLayerThickness = 1e-10;
    std::function<void(std::vector<dd4hep::DetElement>& detectors)>
        sortDetectors = sortFCChhDetElements;
    /// Material decorator
    std::shared_ptr<const Acts::IMaterialDecorator> materialDecorator;

    /// Optional geometry identifier hook to be used during closure
    std::shared_ptr<const Acts::GeometryIdentifierHook> geometryIdentifierHook =
        std::make_shared<const Acts::GeometryIdentifierHook>();
  };

  explicit DD4hepDetector(const Config& cfg);
  DD4hepDetector(const DD4hepDetector&) = delete;
  DD4hepDetector(DD4hepDetector&&);
  ~DD4hepDetector() override;
  DD4hepDetector& operator=(const DD4hepDetector&) = delete;
  DD4hepDetector& operator=(DD4hepDetector&&);

  /// Interface method to access to the DD4hep geometry
  dd4hep::Detector& dd4hepDetector();

  /// Interface method to access the DD4hep geometry
  /// @return The world DD4hep DetElement
  dd4hep::DetElement dd4hepGeometry();

  /// Interface method to Access the TGeo geometry
  /// @return The world TGeoNode (physical volume)
  TGeoNode& tgeoGeometry();

  Gen1GeometryHolder buildGen1Geometry() override;

  std::unique_ptr<G4VUserDetectorConstruction> buildGeant4DetectorConstruction(
      std::vector<std::shared_ptr<RegionCreator>> regionCreators) override;

 private:
  /// Private method to initiate building of the DD4hep geometry
  void buildDD4hepGeometry();

  /// The config class
  Config m_cfg;
  /// Pointer to the interface to the DD4hep geometry
  std::unique_ptr<dd4hep::Detector> m_detector;
};

}  // namespace ActsExamples
