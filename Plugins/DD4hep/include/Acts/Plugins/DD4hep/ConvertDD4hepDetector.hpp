// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <DD4hep/DetElement.h>

namespace Acts {

/// Sort function which sorts dd4hep::DetElement by their ID
/// @param [in,out] det the dd4hep::DetElements to be sorted
inline void sortDetElementsByID(std::vector<dd4hep::DetElement>& det) {
  sort(det.begin(), det.end(),
       [](const dd4hep::DetElement& a, const dd4hep::DetElement& b) {
         return (a.id() < b.id());
       });
}

/// @brief Global method which creates the TrackingGeometry from DD4hep input
///
/// This method returns a std::unique_ptr of the TrackingGeometry from the
/// World DD4hep \a DetElement.
/// @pre Before using this method make sure, that the preconditions described in
/// DD4hepPlugins are met.
///
///
/// @param [in] worldDetElement the DD4hep DetElement of the world
/// @param [in] loggingLevel is the debug logging level of the conversion and
/// geometry building
/// @param [in] bTypePhi is how the sensitive surfaces (modules) should be
/// binned in a layer in phi direction.
/// @note Possible binningtypes:
/// 	- arbitrary   - of the sizes if the surfaces and the distance inbetween
/// 		vary. This mode finds out the bin boundaries by scanning through
/// the 		surfaces.
/// 	- equidistant - if the sensitive surfaces are placed equidistantly
/// @note equidistant binningtype is recommended because it is faster not only
/// while building the geometry but also for look up during the extrapolation
/// @param [in] bTypeR is how the sensitive surfaces (modules) should be binned
/// in a layer in r direction
/// @param [in] bTypeZ is how the sensitive surfaces (modules) should be binned
/// in a layer in z direction
/// @param [in] layerEnvelopeR the tolerance added to the geometrical extension
/// in r of the layers contained to build the volume envelope around
/// @param [in] layerEnvelopeZ the tolerance added to the geometrical extension
/// in z of the layers contained to build the volume envelope around
/// @param [in] defaultLayerThickness In case no surfaces (to be contained by
/// the layer) are handed over, the layer thickness will be set to this value
/// @note Layers containing surfaces per default are not allowed to be
///       attached to each other (navigation will fail at this point).
///       However, to allow material layers (not containing surfaces) to be
///       attached to each other, this default thickness is needed. In this
///       way, the layer will be thin (with space to the next layer), but
///       the material will have the'real' thickness.
/// @attention The default thickness should be set thin enough that no
///            touching or overlapping with the next layer can happen.
/// @param [in] sortSubDetectors @c std::function which should be used in order to
///                              sort all sub detectors (=all Detelements
///                              collected by the method
///                              @c collectSubDetectors() ) from bottom to top to
///                              ensure correct wrapping of the volumes, which
///                              is needed for navigation. Therefore the
///                              different hierachies need to be sorted
///                              ascending. The default is sorting by ID.
/// @param gctx The geometry context to use
/// @param matDecorator is the material decorator that loads material maps
///
/// @exception std::logic_error if an error in the translation occurs
/// @return std::unique_ptr to the full TrackingGeometry

///	* The Tracking geometry needs to be built from bottom to top to ensure
/// Navigation. Therefore the different hierachies need to be sorted ascending.
/// Per default the sub detectors are sorted by the id of their
/// dd4hep::DetElement. In case another sorting needs to be applied, the users
/// can provide their own function

std::unique_ptr<const TrackingGeometry> convertDD4hepDetector(
    dd4hep::DetElement worldDetElement,
    Logging::Level loggingLevel = Logging::Level::INFO,
    BinningType bTypePhi = equidistant, BinningType bTypeR = equidistant,
    BinningType bTypeZ = equidistant, double layerEnvelopeR = UnitConstants::mm,
    double layerEnvelopeZ = UnitConstants::mm,
    double defaultLayerThickness = UnitConstants::fm,
    const std::function<void(std::vector<dd4hep::DetElement>& detectors)>&
        sortSubDetectors = sortDetElementsByID,
    const GeometryContext& gctx = GeometryContext(),
    std::shared_ptr<const IMaterialDecorator> matDecorator = nullptr);

/// @brief Method internally used to create an Acts::CylinderVolumeBuilder
///
/// This method creates an Acts::CylinderVolumeBuilder from a sub detector
/// (= 'compound' DetElement or DetElements which are declared as 'isBarrel' or
/// 'isBeampipe' by their extension.
///
///
/// @param [in] subDetector the DD4hep DetElement of the subdetector
/// @param [in] loggingLevel is the debug logging level of the conversion and
/// geometry building
/// @param [in] bTypePhi is how the sensitive surfaces (modules) should be
/// binned in a layer in phi direction.
/// @note Possible binningtypes:
/// 	- arbitrary   - of the sizes if the surfaces and the distance inbetween
/// 		vary. This mode finds out the bin boundaries by scanning through
/// the 		surfaces.
/// 	- equidistant - if the sensitive surfaces are placed equidistantly
/// @note equidistant binningtype is recommended because it is faster not only
/// while building the geometry but also for look up during the extrapolation
/// @param [in] bTypeR is how the sensitive surfaces (modules) should be binned
/// in a layer in r direction
/// @param [in] bTypeZ is how the sensitive surfaces (modules) should be binned
/// in a layer in z direction
/// @param [in] layerEnvelopeR the tolerance added to the geometrical extension
/// in r of the layers contained to build the volume envelope around
/// @param [in] layerEnvelopeZ the tolerance added to the geometrical extension
/// in z of the layers contained to build the volume envelope around
/// @param [in] defaultLayerThickness In case no surfaces (to be contained by
/// the layer) are handed over, the layer thickness will be set to this value
/// @note Layers containing surfaces per default are not allowed to be
///       attached to each other (navigation will fail at this point).
///       However, to allow material layers (not containing surfaces) to be
///       attached to each other, this default thickness is needed. In this
///       way, the layer will be thin (with space to the next layer), but
///       the material will have the'real' thickness.
/// @attention The default thickness should be set thin enough that no
///            touching or overlapping with the next layer can happen.
/// @return std::shared_ptr the Acts::CylinderVolumeBuilder which can be used to
/// build the full tracking geometry
std::shared_ptr<const CylinderVolumeBuilder> volumeBuilder_dd4hep(
    dd4hep::DetElement subDetector,
    Logging::Level loggingLevel = Logging::Level::INFO,
    BinningType bTypePhi = equidistant, BinningType bTypeR = equidistant,
    BinningType bTypeZ = equidistant, double layerEnvelopeR = UnitConstants::mm,
    double layerEnvelopeZ = UnitConstants::mm,
    double defaultLayerThickness = UnitConstants::fm);

/// Helper method internally used to create a default
/// Acts::CylinderVolumeBuilder
std::shared_ptr<const CylinderVolumeHelper> cylinderVolumeHelper_dd4hep(
    Logging::Level loggingLevel = Logging::Level::INFO);

/// Method internally used by convertDD4hepDetector to collect all sub detectors
/// Sub detector means each 'compound' DetElement or DetElements which are
/// declared as 'isBarrel' or 'isBeampipe' by their extension.
/// @param [in] detElement the dd4hep::DetElement of the volume of which the sub
/// detectors should be collected
/// @param [out] subdetectors the DD4hep::DetElements of the sub detectors
/// contained by detElement
void collectSubDetectors_dd4hep(dd4hep::DetElement& detElement,
                                std::vector<dd4hep::DetElement>& subdetectors);

/// Method internally used by convertDD4hepDetector to collect all volumes of a
/// compound detector
/// @param [in] detElement the dd4hep::DetElement of the volume of which the
/// compounds should be collected
/// @param [out] compounds the DD4hep::DetElements of the compounds contained by
/// detElement
void collectCompounds_dd4hep(dd4hep::DetElement& detElement,
                             std::vector<dd4hep::DetElement>& compounds);

/// Method internally used by convertDD4hepDetector
/// @param [in] detElement the dd4hep::DetElement of the volume of which the
/// layers should be collected
/// @param [out] layers the DD4hep::DetElements of the layers contained by
/// detElement
void collectLayers_dd4hep(dd4hep::DetElement& detElement,
                          std::vector<dd4hep::DetElement>& layers);

/// Method internally used by convertDD4hepDetector
/// @param [in] detElement the dd4hep::DetElement of the volume of which the
/// volumes should be collected
/// @param [out] volumes the DD4hep::DetElements of the volumes contained by
/// detElement
void collectVolumes_dd4hep(dd4hep::DetElement& detElement,
                           std::vector<dd4hep::DetElement>& volumes);
}  // namespace Acts
