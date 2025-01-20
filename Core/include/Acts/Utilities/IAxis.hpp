// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

#include <iosfwd>
#include <memory>
#include <vector>

namespace Acts {

/// Common base class for all Axis instance. This allows generice handling
/// such as for inspection.
class IAxis {
 public:
  /// Virtual destructor
  virtual ~IAxis() = default;

  /// @brief returns whether the axis is equidistant
  ///
  /// @return bool is equidistant
  virtual bool isEquidistant() const = 0;

  /// @brief returns whether the axis is variable
  ///
  /// @return bool is variable
  virtual bool isVariable() const = 0;

  /// @brief returns the type of the axis
  /// @return @c AxisType of this axis
  virtual AxisType getType() const = 0;

  /// @brief returns the boundary type set in the template param
  ///
  /// @return @c AxisBoundaryType of this axis
  virtual AxisBoundaryType getBoundaryType() const = 0;

  /// @brief Return a vector of bin edges
  /// @return Vector which contains the bin edges
  virtual std::vector<double> getBinEdges() const = 0;

  /// @brief get minimum of binning range
  ///
  /// @return minimum of binning range
  virtual double getMin() const = 0;

  /// @brief get maximum of binning range
  ///
  /// @return maximum of binning range
  virtual double getMax() const = 0;

  /// @brief get total number of bins
  ///
  /// @return total number of bins (excluding under-/overflow bins)
  virtual std::size_t getNBins() const = 0;

  /// Centralized axis factory for equidistant binning
  ///
  /// @param aBoundaryType the axis boundary type
  /// @param min the minimum edge of the axis
  /// @param max the maximum edge of the axis
  /// @param nbins the number of bins
  ///
  /// @throws std::invalid_argument if min >= max or nbins == 0
  ///
  /// @return a unique pointer to the axis
  static std::unique_ptr<IAxis> create(AxisBoundaryType aBoundaryType,
                                       double min, double max,
                                       std::size_t nbins);

  /// Centralized axis factory for variable binning
  ///
  /// @param aBoundaryType the axis boundary type
  /// @param edges are the bin edges
  ///
  /// @throws std::invalid_argument if edges is empty or not strictly increasing
  ///
  /// @return a unique pointer to the axis
  static std::unique_ptr<IAxis> create(AxisBoundaryType aBoundaryType,
                                       const std::vector<double>& edges);

  /// Helper function that dispatches from the @c IAxis base class
  /// to a concrete axis type. It will call the provided @p callable
  /// with a const reference to the concrete axis type.
  /// @tparam callable_t the callable type
  /// @param callable the callable object
  /// @return the value returned by the callable
  template <typename callable_t>
  decltype(auto) visit(const callable_t& callable) const {
    auto switchOnType =
        [this, &callable]<AxisBoundaryType bdt>(AxisBoundaryTypeTag<bdt>) {
          switch (getType()) {
            using enum AxisType;
            case Equidistant:
              return callable(
                  dynamic_cast<const Axis<AxisType::Equidistant, bdt>&>(*this));
            case Variable:
              return callable(
                  dynamic_cast<const Axis<AxisType::Variable, bdt>&>(*this));
            default:
              throw std::logic_error("Unknown axis type");
          }
        };

    switch (getBoundaryType()) {
      using enum AxisBoundaryType;
      case Open:
        return switchOnType(AxisOpen);
      case Bound:
        return switchOnType(AxisBound);
      case Closed:
        return switchOnType(AxisClosed);
      default:
        throw std::logic_error("Unknown axis type");
    }
  }

  /// Check if two axes are equal
  /// @param lhs first axis
  /// @param rhs second axis
  /// @return true if the axes are equal
  friend bool operator==(const IAxis& lhs, const IAxis& rhs) {
    return lhs.getType() == rhs.getType() &&
           lhs.getBoundaryType() == rhs.getBoundaryType() &&
           lhs.getMin() == rhs.getMin() && lhs.getMax() == rhs.getMax() &&
           lhs.getNBins() == rhs.getNBins() &&
           lhs.getBinEdges() == rhs.getBinEdges();
  }

  /// Output stream operator
  /// @param os output stream
  /// @param axis the axis to be printed
  /// @return the output stream
  friend std::ostream& operator<<(std::ostream& os, const IAxis& axis) {
    axis.toStream(os);
    return os;
  }

 protected:
  /// Dispatch to the correct stream operator
  /// @param os output stream
  virtual void toStream(std::ostream& os) const = 0;
};

template <typename T>
concept AxisConcept = std::derived_from<T, IAxis>;

}  // namespace Acts
