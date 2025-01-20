// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/IAxis.hpp"

namespace Acts {
/// @brief Description of a ProtoAxis which holds an IAxis
/// and lets the user to define a certain axis type, boundary type
/// and associated binning.
///
/// The IAxis allows via the visitor pattern to access the actual axis type
/// which helps to create grid creation code by the compiler as done
/// in the makeGrid helper functions.
///
/// In addition to a simple axis definitions, it holds also a description
/// of the axis direction.
class ProtoAxis {
 public:
  /// Convenience constructors - for variable binning
  ///
  /// @param aDir the value/cast in which this is binned
  /// @param abType the axis boundary type
  /// @param edges the bin edges (variable binning)
  ProtoAxis(AxisDirection aDir, Acts::AxisBoundaryType abType,
            const std::vector<double>& edges);

  /// Convenience constructors - for equidistant binning
  ///
  /// @param aDir the value/cast in which this is binned
  /// @param abType the axis boundary type
  /// @param minE the lowest edge of the binning
  /// @param maxE the highest edge of the binning
  /// @param nbins the number of bins
  ProtoAxis(AxisDirection aDir, AxisBoundaryType abType, double minE,
            double maxE, std::size_t nbins);

  /// Placeholder constructor for auto-range binning
  ///
  /// @param aDir the value/cast in which this is binned
  /// @param abType the axis boundary type
  /// @param nbins the number of bins
  ///
  /// @note that auto-range is only supported for equidistant binning
  ProtoAxis(AxisDirection aDir, AxisBoundaryType abType, std::size_t nbins);

  ProtoAxis() = delete;
  ProtoAxis(const ProtoAxis&) = delete;
  ProtoAxis& operator=(const ProtoAxis&) = delete;
  ProtoAxis(ProtoAxis&&) = default;
  ProtoAxis& operator=(ProtoAxis&&) = default;

  /// @brief returns the axis direction
  ///
  /// @return @c AxisDirection of this axis
  AxisDirection getAxisDirection() const;

  /// @brief Return the IAxis representation
  ///
  /// @return @c AxisType of this axis
  const IAxis& getAxis() const;

  /// @brief check if this is an auto-range binning
  bool isAutorange() const;

  /// Dump into a string
  /// @return the string representation
  std::string toString() const;

 protected:
  /// Dispatch to the correct stream operator
  /// @param os output stream
  void toStream(std::ostream& os) const;

  /// The axis direction
  AxisDirection m_axisDir = AxisDirection::AxisX;

  /// The axis representation
  std::unique_ptr<IAxis> m_axis = nullptr;

  /// Indicate if this is a place holder auto-range binning
  bool m_autorange = false;
};

/// @brief Helper method to create a 1D grid from a single proto axis
///
/// @tparam payload_t the grid payloat type
///
/// @param a the proto axis
///
/// @return an IGrid unique ptr and hence transfers ownership
template <typename payload_t>
std::unique_ptr<IGrid> makeGrid(const ProtoAxis& a) {
  return a.getAxis().visit([&](const auto& axis) -> std::unique_ptr<IGrid> {
    using AxisTypeA = std::decay_t<decltype(axis)>;
    using GridType = Grid<payload_t, AxisTypeA>;
    return std::make_unique<GridType>(axis);
  });
}

/// @brief Helper method to create a 2D grid from a two proto axes
///
/// @tparam payload_t the grid payloat type
///
/// @param a the first proto axis
/// @param b the second proto axis
///
/// @return an IGrid unique ptr and hence transfers ownership
template <typename payload_t>
std::unique_ptr<IGrid> makeGrid(const ProtoAxis& a, const ProtoAxis& b) {
  return a.getAxis().visit([&](const auto& axisA) -> std::unique_ptr<IGrid> {
    return b.getAxis().visit([&](const auto& axisB) -> std::unique_ptr<IGrid> {
      using AxisTypeA = std::decay_t<decltype(axisA)>;
      using AxisTypeB = std::decay_t<decltype(axisB)>;
      using GridType = Grid<payload_t, AxisTypeA, AxisTypeB>;
      return std::make_unique<GridType>(axisA, axisB);
    });
  });
}

}  // namespace Acts
