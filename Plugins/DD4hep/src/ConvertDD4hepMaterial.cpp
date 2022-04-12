// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/ConvertDD4hepMaterial.hpp"

#include "Acts/Geometry/ApproachDescriptor.hpp"
#include "Acts/Geometry/CylinderLayer.hpp"
#include "Acts/Geometry/DiscLayer.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinUtility.hpp"

#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

#include "XML/XMLElements.h"

std::shared_ptr<Acts::ProtoSurfaceMaterial> Acts::createProtoMaterial(
    const ActsExtension& actsExtension, const std::string& valueTag,
    const std::vector<std::pair<const std::string, Acts::BinningOption> >&
        binning) {
  // Create the bin utility
  Acts::BinUtility bu;
  // Loop over the bins
  for (auto& bin : binning) {
    // finding the iterator position to determine the binning value
    auto bit = std::find(Acts::binningValueNames().begin(),
                         Acts::binningValueNames().end(), bin.first);
    size_t indx = std::distance(Acts::binningValueNames().begin(), bit);
    Acts::BinningValue bval = Acts::BinningValue(indx);
    Acts::BinningOption bopt = bin.second;
    double min = 0.;
    double max = 0.;
    if (bopt == Acts::closed) {
      min = -M_PI;
      max = M_PI;
    }
    int bins = actsExtension.getValue(bin.first, valueTag);
    if (bins >= 1) {
      bu += Acts::BinUtility(bins, min, max, bopt, bval);
    }
  }
  return std::make_shared<Acts::ProtoSurfaceMaterial>(bu);
}

void Acts::addLayerProtoMaterial(
    const ActsExtension& actsExtension, Layer& layer,
    const std::vector<std::pair<const std::string, Acts::BinningOption> >&
        binning) {
  // Start with the representing surface
  std::vector<std::string> materialOptions = {"layer_material_representing"};
  std::vector<const Surface*> materialSurfaces = {
      &(layer.surfaceRepresentation())};
  // Now fill (optionally) with the approach surfaces
  auto aDescriptor = layer.approachDescriptor();
  if (aDescriptor != nullptr and aDescriptor->containedSurfaces().size() >= 2) {
    // Add the inner and outer approach surface
    const std::vector<const Surface*>& aSurfaces =
        aDescriptor->containedSurfaces();
    materialOptions.push_back("layer_material_inner");
    materialSurfaces.push_back(aSurfaces[0]);
    materialOptions.push_back("layer_material_outer");
    materialSurfaces.push_back(aSurfaces[1]);
  }

  // Now loop over it and create the ProtoMaterial
  for (unsigned int is = 0; is < materialOptions.size(); ++is) {
    if (actsExtension.hasValue(materialOptions[is])) {
      // Create the material and assign it
      auto psMaterial =
          createProtoMaterial(actsExtension, materialOptions[is], binning);
      // const_cast (ugly - to be changed after internal geometry stored
      // non-const)
      Surface* surface = const_cast<Surface*>(materialSurfaces[is]);
      surface->assignSurfaceMaterial(psMaterial);
    }
  }
}

void Acts::addCylinderLayerProtoMaterial(dd4hep::DetElement detElement,
                                         Layer& cylinderLayer,
                                         Logging::Level loggingLevel) {
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("DD4hepConversion", loggingLevel));
  ACTS_VERBOSE(
      "Translating DD4hep material into Acts material for CylinderLayer : "
      << detElement.name());
  // Get the Acts extension in order to prepare the ProtoMaterial
  auto actsExtension = detElement.extension<ActsExtension>();
  if (actsExtension != nullptr and actsExtension->hasType("layer_material")) {
    addLayerProtoMaterial(*actsExtension, cylinderLayer,
                          {{"binPhi", Acts::closed}, {"binZ", Acts::open}});
  }
}

void Acts::addDiscLayerProtoMaterial(dd4hep::DetElement detElement,
                                     Layer& discLayer,
                                     Logging::Level loggingLevel) {
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("DD4hepConversion", loggingLevel));
  ACTS_VERBOSE("Translating DD4hep material into Acts material for DiscLayer : "
               << detElement.name());

  // Get the Acts extension
  auto actsExtension = detElement.extension<ActsExtension>();
  if (actsExtension != nullptr and actsExtension->hasType("layer_material")) {
    addLayerProtoMaterial(*actsExtension, discLayer,
                          {{"binPhi", Acts::closed}, {"binR", Acts::open}});
  }
}
