// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/MaterialJsonConverter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/InterpolatedMaterialMap.hpp"
#include "Acts/Material/MaterialGridHelper.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "Acts/Plugins/Json/GeometryJsonKeys.hpp"
#include "Acts/Plugins/Json/GridJsonConverter.hpp"
#include "Acts/Plugins/Json/UtilitiesJsonConverter.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridAccessHelpers.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"
#include "Acts/Utilities/TypeList.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iosfwd>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

// Grid definition : eq bound
template <typename value_type>
using GridEqBound =
    Acts::Grid<value_type, Acts::Axis<Acts::AxisType::Equidistant,
                                      Acts::AxisBoundaryType::Bound>>;
// Grid definition : eq closed
template <typename value_type>
using GridEqClosed =
    Acts::Grid<value_type, Acts::Axis<Acts::AxisType::Equidistant,
                                      Acts::AxisBoundaryType::Closed>>;

// Grid definition : eq bound eq bound
template <typename value_type>
using GridEqBoundEqBound = Acts::Grid<
    value_type,
    Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Bound>,
    Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Bound>>;

// Grid definition : eq bound eq closed
template <typename value_type>
using GridEqBoundEqClosed = Acts::Grid<
    value_type,
    Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Bound>,
    Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Closed>>;

// Grid definition : eq closed eq bound
template <typename value_type>
using GridEqClosedEqBound = Acts::Grid<
    value_type,
    Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Closed>,
    Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Bound>>;

/// @brief Helper function to convert a grid surface material to json
///
/// @tparam indexed_grid_materital_t
/// @param jMaterial the json object to written into
/// @param indexedMaterialCandidate the actual indexed material
template <typename indexed_grid_materital_t>
void convertIndexedGridMaterial(
    nlohmann::json& jMaterial,
    const Acts::ISurfaceMaterial& indexedMaterialCandidate) {
  // Check if the material is of the right type
  const indexed_grid_materital_t* indexedMaterial =
      dynamic_cast<const indexed_grid_materital_t*>(&indexedMaterialCandidate);

  if (indexedMaterial != nullptr) {
    // It is a grid type material
    jMaterial[Acts::JsonKey::kType] = "grid";
    nlohmann::json jMaterialAccessor;
    // Assume globally indexed first
    jMaterialAccessor["type"] = "globally_indexed";

    // If we have a globally indexed map, the material data is loaded elsewhere,
    // locally indexed material vectors are written though
    const auto& materialAccessor = indexedMaterial->materialAccessor();

    if constexpr (std::is_same_v<decltype(materialAccessor),
                                 const Acts::IndexedMaterialAccessor&>) {
      // It's actually locally indexed
      jMaterialAccessor["type"] = "indexed";

      nlohmann::json jMaterialData;
      for (const auto& msl : materialAccessor.material) {
        jMaterialData.push_back(msl);
      }
      jMaterialAccessor["storage_vector"] = jMaterialData;
    }
    // Write the index grid
    jMaterialAccessor["grid"] =
        Acts::GridJsonConverter::toJson(indexedMaterial->grid());
    jMaterial["accessor"] = jMaterialAccessor;

    // Global and bound -> grid local
    jMaterial["global_to_grid_local"] = Acts::GridAccessJsonConverter::toJson(
        *(indexedMaterial->globalToGridLocal().instance()));

    jMaterial["bound_to_grid_local"] = Acts::GridAccessJsonConverter::toJson(
        *(indexedMaterial->boundToGridLocal().instance()));
  }
}

/// @brief Unrolling function for catching the right instance
///
/// @param jMaterial is the json object to be written into
/// @param indexedMaterial is the indexed material
template <typename... Args>
void unrollIndexedGridConversion(nlohmann::json& jMaterial,
                                 const Acts::ISurfaceMaterial& indexedMaterial,
                                 Acts::TypeList<Args...> /*unused*/) {
  (convertIndexedGridMaterial<Args>(jMaterial, indexedMaterial), ...);
}

template <typename IndexedAccessorType>
Acts::ISurfaceMaterial* indexedMaterialFromJson(nlohmann::json& jMaterial) {
  // Load accessor and grid
  nlohmann::json jMaterialAccessor = jMaterial["accessor"];

  // Prepare the material and its accessor
  IndexedAccessorType materialAccessor{};

  // If it's locally indexed, we need to load the material vector
  if constexpr (std::is_same_v<IndexedAccessorType,
                               Acts::IndexedMaterialAccessor>) {
    // It's actually locally indexed
    for (const auto& msl : jMaterialAccessor["storage_vector"]) {
      materialAccessor.material.push_back(msl);
    }
  }

  // Now make the grid and the axes
  nlohmann::json jGrid = jMaterialAccessor["grid"];
  nlohmann::json jGridAxes = jGrid["axes"];

  Acts::AxisBoundaryType boundaryType0 = jGridAxes[0]["boundary_type"];

  // 1-dimensional case
  if (jGridAxes.size() == 1u) {
    // Bound case
    if (boundaryType0 == Acts::AxisBoundaryType::Bound) {
      Acts::GridAxisGenerators::EqBound eqBound{jGridAxes[0]["range"],
                                                jGridAxes[0]["bins"]};
      auto grid =
          Acts::GridJsonConverter::fromJson<decltype(eqBound), std::size_t>(
              jGrid, eqBound);

      auto boundToGridLocal =
          Acts::GridAccessJsonConverter::boundToGridLocal1DimDelegateFromJson(
              jMaterial["bound_to_grid_local"]);

      auto globalToGridLocal =
          Acts::GridAccessJsonConverter::globalToGridLocal1DimDelegateFromJson(
              jMaterial["global_to_grid_local"]);

      return new Acts::IndexedSurfaceMaterial<decltype(grid)>(
          std::move(grid), std::move(materialAccessor),
          std::move(boundToGridLocal), std::move(globalToGridLocal));
    }
    // Closed case
    if (boundaryType0 == Acts::AxisBoundaryType::Closed) {
      Acts::GridAxisGenerators::EqClosed eqClosed{jGridAxes[0]["range"],
                                                  jGridAxes[0]["bins"]};
      auto grid =
          Acts::GridJsonConverter::fromJson<decltype(eqClosed), std::size_t>(
              jGrid, eqClosed);

      auto boundToGridLocal =
          Acts::GridAccessJsonConverter::boundToGridLocal1DimDelegateFromJson(
              jMaterial["bound_to_grid_local"]);

      auto globalToGridLocal =
          Acts::GridAccessJsonConverter::globalToGridLocal1DimDelegateFromJson(
              jMaterial["global_to_grid_local"]);

      return new Acts::IndexedSurfaceMaterial<decltype(grid)>(
          std::move(grid), std::move(materialAccessor),
          std::move(boundToGridLocal), std::move(globalToGridLocal));
    }
  }

  // 2-dimensional case
  if (jGridAxes.size() == 2u) {
    // Second boundary type
    Acts::AxisBoundaryType boundaryType1 = jGridAxes[1]["boundary_type"];

    // Bound-bound setup
    if (boundaryType0 == Acts::AxisBoundaryType::Bound &&
        boundaryType1 == Acts::AxisBoundaryType::Bound) {
      Acts::GridAxisGenerators::EqBoundEqBound eqBoundEqBound{
          jGridAxes[0]["range"], jGridAxes[0]["bins"], jGridAxes[1]["range"],
          jGridAxes[1]["bins"]};
      auto grid =
          Acts::GridJsonConverter::fromJson<decltype(eqBoundEqBound),
                                            std::size_t>(jGrid, eqBoundEqBound);

      auto boundToGridLocal =
          Acts::GridAccessJsonConverter::boundToGridLocal2DimDelegateFromJson(
              jMaterial["bound_to_grid_local"]);

      auto globalToGridLocal =
          Acts::GridAccessJsonConverter::globalToGridLocal2DimDelegateFromJson(
              jMaterial["global_to_grid_local"]);

      return new Acts::IndexedSurfaceMaterial<decltype(grid)>(
          std::move(grid), std::move(materialAccessor),
          std::move(boundToGridLocal), std::move(globalToGridLocal));
    }

    // Bound-closed setup
    if (boundaryType0 == Acts::AxisBoundaryType::Bound &&
        boundaryType1 == Acts::AxisBoundaryType::Closed) {
      Acts::GridAxisGenerators::EqBoundEqClosed eqBoundEqClosed{
          jGridAxes[0]["range"], jGridAxes[0]["bins"], jGridAxes[1]["range"],
          jGridAxes[1]["bins"]};
      auto grid = Acts::GridJsonConverter::fromJson<decltype(eqBoundEqClosed),
                                                    std::size_t>(
          jGrid, eqBoundEqClosed);

      auto boundToGridLocal =
          Acts::GridAccessJsonConverter::boundToGridLocal2DimDelegateFromJson(
              jMaterial["bound_to_grid_local"]);

      auto globalToGridLocal =
          Acts::GridAccessJsonConverter::globalToGridLocal2DimDelegateFromJson(
              jMaterial["global_to_grid_local"]);

      return new Acts::IndexedSurfaceMaterial<decltype(grid)>(
          std::move(grid), std::move(materialAccessor),
          std::move(boundToGridLocal), std::move(globalToGridLocal));
    }

    // Closed-bound setup
    if (boundaryType0 == Acts::AxisBoundaryType::Closed &&
        boundaryType1 == Acts::AxisBoundaryType::Bound) {
      Acts::GridAxisGenerators::EqClosedEqBound eqClosedEqBound{
          jGridAxes[0]["range"], jGridAxes[0]["bins"], jGridAxes[1]["range"],
          jGridAxes[1]["bins"]};
      auto grid = Acts::GridJsonConverter::fromJson<decltype(eqClosedEqBound),
                                                    std::size_t>(
          jGrid, eqClosedEqBound);

      auto boundToGridLocal =
          Acts::GridAccessJsonConverter::boundToGridLocal2DimDelegateFromJson(
              jMaterial["bound_to_grid_local"]);

      auto globalToGridLocal =
          Acts::GridAccessJsonConverter::globalToGridLocal2DimDelegateFromJson(
              jMaterial["global_to_grid_local"]);

      return new Acts::IndexedSurfaceMaterial<decltype(grid)>(
          std::move(grid), std::move(materialAccessor),
          std::move(boundToGridLocal), std::move(globalToGridLocal));
    }
  }

  return nullptr;
}

}  // namespace

void Acts::to_json(nlohmann::json& j, const Material& t) {
  if (!t.isValid()) {
    return;
  }
  for (unsigned i = 0; i < t.parameters().size(); ++i) {
    j.push_back(t.parameters()[i]);
  }
}

void Acts::from_json(const nlohmann::json& j, Material& t) {
  if (j.is_null()) {
    return;
  }
  Acts::Material::ParametersVector params =
      Acts::Material::ParametersVector::Zero();
  for (auto i = params.size(); 0 < i--;) {
    // .at(...) ensures bound checks
    params[i] = j.at(i);
  }
  t = Acts::Material(params);
  return;
}

void Acts::to_json(nlohmann::json& j, const MaterialSlab& t) {
  nlohmann::json jmat(t.material());
  j["material"] = jmat;
  j["thickness"] = t.thickness();
}

void Acts::from_json(const nlohmann::json& j, MaterialSlab& t) {
  Material mat(j["material"].get<Material>());
  t = Acts::MaterialSlab(mat, j.at("thickness").get<float>());
}

void Acts::from_json(const nlohmann::json& j, MaterialSlabMatrix& t) {
  // the input data must be array[array[object]]
  for (auto& outer : j) {
    Acts::MaterialSlabVector mpVector;
    for (auto& inner : outer) {
      MaterialSlab mat = inner.get<MaterialSlab>();
      mpVector.emplace_back(mat);
    }
    t.push_back(std::move(mpVector));
  }
}

void Acts::to_json(nlohmann::json& j, const surfaceMaterialPointer& material) {
  nlohmann::json jMaterial;
  // A bin utility needs to be written
  const Acts::BinUtility* bUtility = nullptr;

  // First: Check if we have a proto material
  auto psMaterial = dynamic_cast<const Acts::ProtoSurfaceMaterial*>(material);
  if (psMaterial != nullptr) {
    // Type is proto material
    jMaterial[Acts::JsonKey::kType] = "proto";
    // Set mapping type
    nlohmann::json mapType(material->mappingType());
    jMaterial[Acts::JsonKey::kMappingType] = mapType;
    // by default the protoMaterial is not used for mapping
    jMaterial[Acts::JsonKey::kMaterialMap] = false;
    // write the bin utility
    bUtility = &(psMaterial->binning());
    // Check in the number of bin is different from 1
    auto& binningData = bUtility->binningData();
    for (std::size_t ibin = 0; ibin < binningData.size(); ++ibin) {
      if (binningData[ibin].bins() > 1) {
        jMaterial[Acts::JsonKey::kMaterialMap] = true;
        break;
      }
    }
    nlohmann::json jBin(*bUtility);
    jMaterial[Acts::JsonKey::kBin] = jBin;
    j[Acts::JsonKey::kMaterial] = jMaterial;
    return;
  }

  // Second: check if we have a homogeneous material
  auto hsMaterial =
      dynamic_cast<const Acts::HomogeneousSurfaceMaterial*>(material);
  if (hsMaterial != nullptr) {
    // type is homogeneous
    jMaterial[Acts::JsonKey::kType] = "homogeneous";
    // Set mapping type
    nlohmann::json mapType(material->mappingType());
    jMaterial[Acts::JsonKey::kMappingType] = mapType;
    // Material has been mapped
    jMaterial[Acts::JsonKey::kMaterialMap] = true;
    nlohmann::json jmat(hsMaterial->materialSlab(Acts::Vector3(0., 0., 0.)));
    jMaterial[Acts::JsonKey::kData] = nlohmann::json::array({
        nlohmann::json::array({
            jmat,
        }),
    });
    j[Acts::JsonKey::kMaterial] = jMaterial;
    return;
  }

  // Next option remaining: BinnedSurface material
  auto bsMaterial = dynamic_cast<const Acts::BinnedSurfaceMaterial*>(material);
  if (bsMaterial != nullptr) {
    // type is binned
    jMaterial[Acts::JsonKey::kType] = "binned";
    // Set mapping type
    nlohmann::json mapType(material->mappingType());
    jMaterial[Acts::JsonKey::kMappingType] = mapType;
    // Material has been mapped
    jMaterial[Acts::JsonKey::kMaterialMap] = true;
    bUtility = &(bsMaterial->binUtility());
    // convert the data
    // get the material matrix
    nlohmann::json mmat = nlohmann::json::array();
    for (const auto& mpVector : bsMaterial->fullMaterial()) {
      nlohmann::json mvec = nlohmann::json::array();
      for (const auto& mp : mpVector) {
        nlohmann::json jmat(mp);
        mvec.push_back(jmat);
      }
      mmat.push_back(std::move(mvec));
    }
    jMaterial[Acts::JsonKey::kData] = std::move(mmat);
    // write the bin utility
    nlohmann::json jBin(*bUtility);
    jMaterial[Acts::JsonKey::kBin] = jBin;
    j[Acts::JsonKey::kMaterial] = jMaterial;
    return;
  }

  // Possible indexed grid types
  using IndexedSurfaceGrids = Acts::TypeList<
      Acts::IndexedSurfaceMaterial<GridEqBound<std::size_t>>,
      Acts::IndexedSurfaceMaterial<GridEqClosed<std::size_t>>,
      Acts::IndexedSurfaceMaterial<GridEqBoundEqBound<std::size_t>>,
      Acts::IndexedSurfaceMaterial<GridEqBoundEqClosed<std::size_t>>,
      Acts::IndexedSurfaceMaterial<GridEqClosedEqBound<std::size_t>>>;

  unrollIndexedGridConversion(jMaterial, *material, IndexedSurfaceGrids{});
  if (!jMaterial.empty()) {
    j[Acts::JsonKey::kMaterial] = jMaterial;
    return;
  }

  // Possible: globally indexed grid types
  using GloballyIndexedSurfaceGrids = Acts::TypeList<
      Acts::GloballyIndexedSurfaceMaterial<GridEqBound<std::size_t>>,
      Acts::GloballyIndexedSurfaceMaterial<GridEqClosed<std::size_t>>,
      Acts::GloballyIndexedSurfaceMaterial<GridEqBoundEqBound<std::size_t>>,
      Acts::GloballyIndexedSurfaceMaterial<GridEqBoundEqClosed<std::size_t>>,
      Acts::GloballyIndexedSurfaceMaterial<GridEqClosedEqBound<std::size_t>>>;

  unrollIndexedGridConversion(jMaterial, *material,
                              GloballyIndexedSurfaceGrids{});
  if (!jMaterial.empty()) {
    j[Acts::JsonKey::kMaterial] = jMaterial;
    return;
  }

  // Possible: material grid types
  // using MaterialSurfaceGrids = Acts::TypeList<
  //    Acts::GridSurfaceMaterial<GridEqBound<std::size_t>>,
  //    Acts::GridSurfaceMaterial<GridEqClosed<std::size_t>>,
  //    Acts::GridSurfaceMaterial<GridEqBoundEqBound<std::size_t>>,
  //    Acts::GridSurfaceMaterial<GridEqBoundEqClosed<std::size_t>>,
  //    Acts::GridSurfaceMaterial<GridEqClosedEqBound<std::size_t>>>;

  // No material the json object is left empty.
  return;
}

void Acts::from_json(const nlohmann::json& j,
                     surfaceMaterialPointer& material) {
  if (j.find(Acts::JsonKey::kMaterial) == j.end()) {
    return;
  }
  nlohmann::json jMaterial = j[Acts::JsonKey::kMaterial];
  // By default no material is return.
  material = nullptr;
  if (jMaterial[Acts::JsonKey::kMaterialMap] == false) {
    return;
  }

  // Grid based material maps
  if (jMaterial[Acts::JsonKey::kType] == "grid") {
    material =
        indexedMaterialFromJson<Acts::IndexedMaterialAccessor>(jMaterial);
    return;
  }

  // The bin utility and material
  Acts::BinUtility bUtility;
  Acts::MaterialSlabMatrix mpMatrix;
  Acts::MappingType mapType = Acts::MappingType::Default;
  for (auto& [key, value] : jMaterial.items()) {
    if (key == Acts::JsonKey::kBin && !value.empty()) {
      from_json(value, bUtility);
    }
    if (key == Acts::JsonKey::kData && !value.empty()) {
      from_json(value, mpMatrix);
    }
    if (key == Acts::JsonKey::kMappingType && !value.empty()) {
      from_json(value, mapType);
    }
  }
  // Return the appropriate typr of material
  if (mpMatrix.empty()) {
    material = new Acts::ProtoSurfaceMaterial(bUtility, mapType);
  } else if (bUtility.bins() == 1) {
    material = new Acts::HomogeneousSurfaceMaterial(mpMatrix[0][0], mapType);
  } else {
    material = new Acts::BinnedSurfaceMaterial(bUtility, mpMatrix, mapType);
  }
}

void Acts::to_json(nlohmann::json& j, const volumeMaterialPointer& material) {
  nlohmann::json jMaterial;
  // A bin utility needs to be written
  const Acts::BinUtility* bUtility = nullptr;
  // Check if we have a proto material
  auto pvMaterial = dynamic_cast<const Acts::ProtoVolumeMaterial*>(material);
  if (pvMaterial != nullptr) {
    // Type is proto material
    jMaterial[Acts::JsonKey::kType] = "proto";
    // By default the protoMaterial is not used for mapping
    jMaterial[Acts::JsonKey::kMaterialMap] = false;
    bUtility = &(pvMaterial->binUtility());
    // Check in the number of bin is different from 1
    auto& binningData = bUtility->binningData();
    for (std::size_t ibin = 0; ibin < binningData.size(); ++ibin) {
      if (binningData[ibin].bins() > 1) {
        jMaterial[Acts::JsonKey::kMaterialMap] = true;
        break;
      }
    }
    // Write the bin utility
    nlohmann::json jBin(*bUtility);
    jMaterial[Acts::JsonKey::kBin] = jBin;
    j[Acts::JsonKey::kMaterial] = jMaterial;
    return;
  }
  // Now check if we have a homogeneous material
  auto hvMaterial =
      dynamic_cast<const Acts::HomogeneousVolumeMaterial*>(material);
  if (hvMaterial != nullptr) {
    // type is homogeneous
    jMaterial[Acts::JsonKey::kType] = "homogeneous";
    jMaterial[Acts::JsonKey::kMaterialMap] = true;
    // array of encoded materials w/ one entry
    nlohmann::json jmat(hvMaterial->material({0, 0, 0}));
    jMaterial[Acts::JsonKey::kData] = nlohmann::json::array({
        jmat,
    });
    j[Acts::JsonKey::kMaterial] = jMaterial;
    return;
  }
  // Only option remaining: material map
  auto bvMaterial2D = dynamic_cast<const Acts::InterpolatedMaterialMap<
      Acts::MaterialMapper<Acts::MaterialGrid2D>>*>(material);
  // Now check if we have a 2D map
  if (bvMaterial2D != nullptr) {
    // type is binned
    jMaterial[Acts::JsonKey::kType] = "interpolated2D";
    jMaterial[Acts::JsonKey::kMaterialMap] = true;
    bUtility = &(bvMaterial2D->binUtility());
    // convert the data
    nlohmann::json mmat = nlohmann::json::array();
    Acts::MaterialGrid2D grid = bvMaterial2D->getMapper().getGrid();
    for (std::size_t bin = 0; bin < grid.size(); bin++) {
      nlohmann::json jmat(Material(grid.at(bin)));
      mmat.push_back(jmat);
    }
    jMaterial[Acts::JsonKey::kData] = std::move(mmat);
    // Write the bin utility
    nlohmann::json jBin(*bUtility);
    jMaterial[Acts::JsonKey::kBin] = jBin;
    j[Acts::JsonKey::kMaterial] = jMaterial;
    return;
  }
  // Only option remaining: material map
  auto bvMaterial3D = dynamic_cast<const Acts::InterpolatedMaterialMap<
      Acts::MaterialMapper<Acts::MaterialGrid3D>>*>(material);
  // Now check if we have a 3D map
  if (bvMaterial3D != nullptr) {
    // type is binned
    jMaterial[Acts::JsonKey::kType] = "interpolated3D";
    jMaterial[Acts::JsonKey::kMaterialMap] = true;
    bUtility = &(bvMaterial3D->binUtility());
    // convert the data
    nlohmann::json mmat = nlohmann::json::array();
    Acts::MaterialGrid3D grid = bvMaterial3D->getMapper().getGrid();
    for (std::size_t bin = 0; bin < grid.size(); bin++) {
      nlohmann::json jmat(Material(grid.at(bin)));
      mmat.push_back(jmat);
    }
    jMaterial[Acts::JsonKey::kData] = std::move(mmat);
    // Write the bin utility
    nlohmann::json jBin(*bUtility);
    jMaterial[Acts::JsonKey::kBin] = jBin;
    j[Acts::JsonKey::kMaterial] = jMaterial;
    return;
  }
}

void Acts::from_json(const nlohmann::json& j, volumeMaterialPointer& material) {
  if (j.find(Acts::JsonKey::kMaterial) == j.end()) {
    return;
  }
  nlohmann::json jMaterial = j[Acts::JsonKey::kMaterial];
  // By default no material is return.
  material = nullptr;
  if (jMaterial[Acts::JsonKey::kMaterialMap] == false) {
    return;
  }
  // The bin utility and material
  Acts::BinUtility bUtility;
  std::vector<Acts::Material> mmat;
  for (auto& [key, value] : jMaterial.items()) {
    if (key == Acts::JsonKey::kBin && !value.empty()) {
      from_json(value, bUtility);
    }
    if (key == Acts::JsonKey::kData && !value.empty()) {
      for (const auto& bin : value) {
        Acts::Material mat(bin.get<Acts::Material>());
        mmat.push_back(mat);
      }
    }
  }
  // We have protoMaterial
  if (mmat.empty()) {
    material = new Acts::ProtoVolumeMaterial(bUtility);
    return;
  }
  if (mmat.size() == 1) {
    material = new Acts::HomogeneousVolumeMaterial(mmat[0]);
    return;
  }
  if (bUtility.dimensions() == 2) {
    std::function<Acts::Vector2(Acts::Vector3)> transfoGlobalToLocal;
    Acts::Grid2D grid = createGrid2D(bUtility, transfoGlobalToLocal);

    Acts::Grid2D::point_t min = grid.minPosition();
    Acts::Grid2D::point_t max = grid.maxPosition();
    Acts::Grid2D::index_t nBins = grid.numLocalBins();

    Acts::EAxis axis1(min[0], max[0], nBins[0]);
    Acts::EAxis axis2(min[1], max[1], nBins[1]);

    // Build the grid and fill it with data
    Acts::MaterialGrid2D mGrid(std::make_tuple(axis1, axis2));

    for (std::size_t bin = 0; bin < mmat.size(); bin++) {
      mGrid.at(bin) = mmat[bin].parameters();
    }
    Acts::MaterialMapper<Acts::MaterialGrid2D> matMap(transfoGlobalToLocal,
                                                      mGrid);
    material = new Acts::InterpolatedMaterialMap<
        Acts::MaterialMapper<Acts::MaterialGrid2D>>(std::move(matMap),
                                                    bUtility);
    return;
  }
  if (bUtility.dimensions() == 3) {
    std::function<Acts::Vector3(Acts::Vector3)> transfoGlobalToLocal;
    Acts::Grid3D grid = createGrid3D(bUtility, transfoGlobalToLocal);

    Acts::Grid3D::point_t min = grid.minPosition();
    Acts::Grid3D::point_t max = grid.maxPosition();
    Acts::Grid3D::index_t nBins = grid.numLocalBins();

    Acts::EAxis axis1(min[0], max[0], nBins[0]);
    Acts::EAxis axis2(min[1], max[1], nBins[1]);
    Acts::EAxis axis3(min[2], max[2], nBins[2]);

    // Build the grid and fill it with data
    Acts::MaterialGrid3D mGrid(std::make_tuple(axis1, axis2, axis3));

    for (std::size_t bin = 0; bin < mmat.size(); bin++) {
      mGrid.at(bin) = mmat[bin].parameters();
    }
    Acts::MaterialMapper<Acts::MaterialGrid3D> matMap(transfoGlobalToLocal,
                                                      mGrid);
    material = new Acts::InterpolatedMaterialMap<
        Acts::MaterialMapper<Acts::MaterialGrid3D>>(std::move(matMap),
                                                    bUtility);
    return;
  }
}

nlohmann::json Acts::MaterialJsonConverter::toJsonDetray(
    const Acts::ISurfaceMaterial& surfaceMaterial, const Acts::Surface& surface,
    std::size_t surfaceIndex, std::map<std::size_t, std::size_t>& gridLink) {
  nlohmann::json jSurfaceMaterial;

  // Binned material conversion
  if (auto binnedMaterial =
          dynamic_cast<const BinnedSurfaceMaterial*>(&surfaceMaterial);
      binnedMaterial != nullptr) {
    // BinUtility modifications
    bool swapped = false;
    // Get the bin utility (make a copy as we may modify it)
    // Detray expects 2-dimensional grid, currently supported are
    // x-y, r-phi, phi-z
    BinUtility bUtility = binnedMaterial->binUtility();
    // Turn the bin value into a 2D grid
    if (bUtility.dimensions() == 1u) {
      if (bUtility.binningData()[0u].binvalue == AxisDirection::AxisR) {
        // Turn to R-Phi
        bUtility += BinUtility(1u, -std::numbers::pi, std::numbers::pi, closed,
                               AxisDirection::AxisPhi);
      } else if (bUtility.binningData()[0u].binvalue == AxisDirection::AxisZ) {
        // Turn to Phi-Z - swap needed
        BinUtility nbUtility(1u, -std::numbers::pi, std::numbers::pi, closed,
                             AxisDirection::AxisPhi);
        nbUtility += bUtility;
        bUtility = std::move(nbUtility);
        swapped = true;
      } else {
        std::runtime_error("Unsupported binning for Detray");
      }
    } else if (bUtility.dimensions() == 2u &&
               bUtility.binningData()[0u].binvalue == AxisDirection::AxisZ &&
               bUtility.binningData()[1u].binvalue == AxisDirection::AxisPhi) {
      BinUtility nbUtility(bUtility.binningData()[1u]);
      nbUtility += BinUtility{bUtility.binningData()[0u]};
      bUtility = std::move(nbUtility);
      swapped = true;
    }

    AxisDirection bVal0 = bUtility.binningData()[0u].binvalue;
    AxisDirection bVal1 = bUtility.binningData()[1u].binvalue;

    // Translate into grid index type
    int gridIndexType = 0;
    if (bVal0 == AxisDirection::AxisR && bVal1 == AxisDirection::AxisPhi) {
      gridIndexType = 0;
    } else if (bVal0 == AxisDirection::AxisPhi &&
               bVal1 == AxisDirection::AxisZ) {
      gridIndexType = 3;
    } else if (bVal0 == AxisDirection::AxisX && bVal1 == AxisDirection::AxisY) {
      gridIndexType = 2;
    } else {
      std::runtime_error("Unsupported binning for Detray");
    }
    // Convert the axes
    nlohmann::json jAxes = toJsonDetray(bUtility, surface);
    // Create  a grid index, i.e. type, index tuple
    nlohmann::json jGridLink;
    jGridLink["type"] = gridIndexType;
    std::size_t gridIndex = 0;
    if (gridLink.contains(gridIndexType)) {
      std::size_t& fGridIndex = gridLink[gridIndex];
      gridIndex = fGridIndex;
      fGridIndex++;
    } else {
      gridLink[gridIndexType] = 1;
    }
    jGridLink["index"] = gridIndex;

    // The grid data
    jSurfaceMaterial["axes"] = jAxes;
    jSurfaceMaterial["grid_link"] = jGridLink;
    jSurfaceMaterial["owner_link"] = surfaceIndex;

    // The bins to be filled
    nlohmann::json jBins;
    auto materialMatrix = binnedMaterial->fullMaterial();
    for (std::size_t ib1 = 0; ib1 < materialMatrix.size(); ++ib1) {
      for (std::size_t ib0 = 0; ib0 < materialMatrix[0u].size(); ++ib0) {
        nlohmann::json jBin;
        // Look up the material slab
        MaterialSlab slab = materialMatrix[ib1][ib0];
        // Translate into a local bin
        std::size_t lb0 = swapped ? ib1 : ib0;
        std::size_t lb1 = swapped ? ib0 : ib1;
        jBin["loc_index"] = std::array<std::size_t, 2u>{lb0, lb1};

        const Material& material = slab.material();
        // The content
        nlohmann::json jContent;
        jContent["thickness"] = slab.thickness();
        // The actual material
        nlohmann::json jMaterialParams;
        if (slab.thickness() > 0.) {
          jMaterialParams["params"] =
              std::vector<double>{material.X0(),
                                  material.L0(),
                                  material.Ar(),
                                  material.Z(),
                                  material.massDensity(),
                                  material.molarDensity(),
                                  0.};

        } else {
          jMaterialParams["params"] =
              std::vector<double>{0., 0., 0., 0., 0., 0., 0.};
        }
        jContent["material"] = jMaterialParams;
        jContent["type"] = 6;
        jContent["surface_idx"] = surfaceIndex;

        nlohmann::json jContentVector;
        jContentVector.push_back(jContent);
        jBin["content"] = jContentVector;
        jBins.push_back(jBin);
      }
    }
    jSurfaceMaterial["bins"] = jBins;
  }
  return jSurfaceMaterial;
}

nlohmann::json Acts::MaterialJsonConverter::toJsonDetray(
    const Acts::BinUtility& binUtility, const Surface& surface) {
  nlohmann::json jAxes;
  for (const auto [ib, bData] : enumerate(binUtility.binningData())) {
    nlohmann::json jAxis;
    jAxis["bounds"] = bData.option == closed ? 2 : 1;
    jAxis["binning"] = 0u;
    jAxis["label"] = ib;
    jAxis["bins"] = bData.bins();
    double offset = 0;
    if (bData.binvalue == AxisDirection::AxisZ) {
      offset = surface.center(Acts::GeometryContext{}).z();
    }
    jAxis["edges"] =
        std::array<double, 2>{bData.min + offset, bData.max + offset};
    jAxes.push_back(jAxis);
  }
  return jAxes;
}
