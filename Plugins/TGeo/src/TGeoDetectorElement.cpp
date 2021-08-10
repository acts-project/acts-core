// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Plugins/TGeo/TGeoSurfaceConverter.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <fstream>
#include <iostream>
#include <utility>

#include <boost/algorithm/string.hpp>

#include "TGeoArb8.h"
#include "TGeoBBox.h"
#include "TGeoBoolNode.h"
#include "TGeoCompositeShape.h"
#include "TGeoTrd2.h"
#include "TGeoTube.h"

using Line2D = Eigen::Hyperplane<double, 2>;

Acts::TGeoDetectorElement::TGeoDetectorElement(
    const Identifier& identifier, const TGeoNode& tGeoNode,
    const TGeoMatrix& tGeoMatrix, const std::string& axes, double scalor,
    std::shared_ptr<const Acts::ISurfaceMaterial> material)
    : Acts::IdentifiedDetectorElement(),
      m_detElement(&tGeoNode),
      m_identifier(identifier) {
  // Create temporary local non const surface (to allow setting the
  // material)
  const Double_t* translation = tGeoMatrix.GetTranslation();
  const Double_t* rotation = tGeoMatrix.GetRotationMatrix();

  auto sensor = m_detElement->GetVolume();
  auto tgShape = sensor->GetShape();

  auto cylinderComps = TGeoSurfaceConverter::cylinderComponents(
      *tgShape, rotation, translation, axes, scalor);
  auto cylinderBounds = std::get<0>(cylinderComps);
  if (cylinderBounds != nullptr) {
    m_transform = std::get<1>(cylinderComps);
    m_bounds = cylinderBounds;
    m_thickness = std::get<2>(cylinderComps);
    m_surface = Surface::makeShared<CylinderSurface>(cylinderBounds, *this);
  }

  // Check next if you do not have a surface
  if (m_surface == nullptr) {
    auto discComps = TGeoSurfaceConverter::discComponents(
        *tgShape, rotation, translation, axes, scalor);
    auto discBounds = std::get<0>(discComps);
    if (discBounds != nullptr) {
      m_bounds = discBounds;
      m_transform = std::get<1>(discComps);
      m_thickness = std::get<2>(discComps);
      m_surface = Surface::makeShared<DiscSurface>(discBounds, *this);
    }
  }

  // Check next if you do not have a surface
  if (m_surface == nullptr) {
    auto planeComps = TGeoSurfaceConverter::planeComponents(
        *tgShape, rotation, translation, axes, scalor);
    auto planeBounds = std::get<0>(planeComps);
    if (planeBounds != nullptr) {
      m_bounds = planeBounds;
      m_transform = std::get<1>(planeComps);
      m_thickness = std::get<2>(planeComps);
      m_surface = Surface::makeShared<PlaneSurface>(planeBounds, *this);
    }
  }

  // set the asscoiated material (non const method)
  if (m_surface != nullptr) {
    m_surface->assignSurfaceMaterial(std::move(material));
  }
}

Acts::TGeoDetectorElement::TGeoDetectorElement(
    const Identifier& identifier, const TGeoNode& tGeoNode,
    const Transform3& tgTransform, std::shared_ptr<const PlanarBounds> tgBounds,
    double tgThickness)
    : Acts::IdentifiedDetectorElement(),
      m_detElement(&tGeoNode),
      m_transform(tgTransform),
      m_identifier(identifier),
      m_bounds(tgBounds),
      m_thickness(tgThickness) {
  m_surface = Surface::makeShared<PlaneSurface>(tgBounds, *this);
}

Acts::TGeoDetectorElement::TGeoDetectorElement(
    const Identifier& identifier, const TGeoNode& tGeoNode,
    const Transform3& tgTransform, std::shared_ptr<const DiscBounds> tgBounds,
    double tgThickness)
    : Acts::IdentifiedDetectorElement(),
      m_detElement(&tGeoNode),
      m_transform(tgTransform),
      m_identifier(identifier),
      m_bounds(tgBounds),
      m_thickness(tgThickness) {
  m_surface = Surface::makeShared<DiscSurface>(tgBounds, *this);
}

Acts::TGeoDetectorElement::~TGeoDetectorElement() = default;
