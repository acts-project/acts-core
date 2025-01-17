// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/IAxis.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace Acts {

// This object can be iterated to produce up to two sequences of integer
// indices, corresponding to the half-open integer ranges [begin1, end1[ and
// [begin2, end2[.
//
// The goal is to emulate the effect of enumerating a range of neighbor
// indices on an axis (which may go out of bounds and wrap around since we
// have AxisBoundaryType::Closed), inserting them into an std::vector, and
// discarding duplicates, without paying the price of duplicate removal
// and dynamic memory allocation in hot magnetic field interpolation code.
//
class NeighborHoodIndices {
 public:
  NeighborHoodIndices() = default;

  NeighborHoodIndices(std::size_t begin, std::size_t end)
      : m_begin1(begin), m_end1(end), m_begin2(end), m_end2(end) {}

  NeighborHoodIndices(std::size_t begin1, std::size_t end1, std::size_t begin2,
                      std::size_t end2)
      : m_begin1(begin1), m_end1(end1), m_begin2(begin2), m_end2(end2) {}

  class iterator {
   public:
    iterator() = default;

    // Specialized constructor for end() iterator
    explicit iterator(std::size_t current)
        : m_current(current), m_wrapped(true) {}

    iterator(std::size_t begin1, std::size_t end1, std::size_t begin2)
        : m_current(begin1),
          m_end1(end1),
          m_begin2(begin2),
          m_wrapped(begin1 == begin2) {}

    std::size_t operator*() const { return m_current; }

    iterator& operator++() {
      ++m_current;
      if (m_current == m_end1) {
        m_current = m_begin2;
        m_wrapped = true;
      }
      return *this;
    }

    bool operator==(const iterator& it) const {
      return (m_current == it.m_current) && (m_wrapped == it.m_wrapped);
    }

   private:
    std::size_t m_current = 0, m_end1 = 0, m_begin2 = 0;
    bool m_wrapped = false;
  };

  iterator begin() const { return iterator(m_begin1, m_end1, m_begin2); }

  iterator end() const { return iterator(m_end2); }

  // Number of indices that will be produced if this sequence is iterated
  std::size_t size() const { return (m_end1 - m_begin1) + (m_end2 - m_begin2); }

  // Collect the sequence of indices into an std::vector
  std::vector<std::size_t> collect() const {
    std::vector<std::size_t> result;
    result.reserve(this->size());
    for (std::size_t idx : *this) {
      result.push_back(idx);
    }
    return result;
  }

 private:
  std::size_t m_begin1 = 0, m_end1 = 0, m_begin2 = 0, m_end2 = 0;
};

/// @brief calculate bin indices for an equidistant binning
///
/// This class provides some basic functionality for calculating bin indices
/// for a given equidistant binning.
template <AxisBoundaryType bdt>
class Axis<AxisType::Equidistant, bdt> : public IAxis {
 public:
  static constexpr AxisType type = AxisType::Equidistant;

  /// @brief default constructor
  ///
  /// @param [in] xmin lower boundary of axis range
  /// @param [in] xmax upper boundary of axis range
  /// @param [in] nBins number of bins to divide the axis range into
  ///
  /// Divide the range \f$[\text{xmin},\text{xmax})\f$ into \f$\text{nBins}\f$
  /// equidistant bins.
  Axis(double xmin, double xmax, std::size_t nBins)
      : m_min(xmin),
        m_max(xmax),
        m_width((xmax - xmin) / nBins),
        m_bins(nBins) {}

  /// Constructor with a tag for the boundary type
  ///
  /// @param [in] typeTag boundary type tag
  /// @param [in] xmin lower boundary of axis range
  /// @param [in] xmax upper boundary of axis range
  /// @param [in] nBins number of bins to divide the axis range into
  ///
  /// Divide the range \f$[\text{xmin},\text{xmax})\f$ into \f$\text{nBins}\f$
  /// equidistant bins.
  Axis(AxisBoundaryTypeTag<bdt> typeTag, double xmin, double xmax,
       std::size_t nBins)
      : Axis(xmin, xmax, nBins) {
    static_cast<void>(typeTag);
  }

  /// @brief returns whether the axis is equidistant
  ///
  /// @return bool is equidistant
  bool isEquidistant() const override { return true; }

  /// @brief returns whether the axis is variable
  ///
  /// @return bool is variable
  bool isVariable() const override { return false; }

  /// @brief returns the type of the axis
  /// @return @c AxisType of this axis
  AxisType getType() const override { return type; }

  /// @brief returns the boundary type set in the template param
  ///
  /// @return @c AxisBoundaryType of this axis
  AxisBoundaryType getBoundaryType() const override { return bdt; }

  /// @brief Get #size bins which neighbor the one given
  ///
  /// Generic overload with symmetric size
  ///
  /// @param [in] idx requested bin index
  /// @param [in] size how many neighboring bins (up/down)
  /// @return Set of neighboring bin indices (global)
  NeighborHoodIndices neighborHoodIndices(std::size_t idx,
                                          std::size_t size = 1) const {
    return neighborHoodIndices(idx,
                               std::make_pair(-static_cast<int>(size), size));
  }

  /// @brief Get #size bins which neighbor the one given
  ///
  /// This is the version for Open
  ///
  /// @param [in] idx requested bin index
  /// @param [in] sizes how many neighboring bins (up/down)
  /// @return Set of neighboring bin indices (global)
  /// @note Open varies given bin and allows 0 and NBins+1 (underflow,
  /// overflow)
  ///       as neighbors
  NeighborHoodIndices neighborHoodIndices(std::size_t idx,
                                          std::pair<int, int> sizes = {-1,
                                                                       1}) const
    requires(bdt == AxisBoundaryType::Open)
  {
    constexpr int min = 0;
    const int max = getNBins() + 1;
    const int itmin = std::clamp(static_cast<int>(idx + sizes.first), min, max);
    const int itmax =
        std::clamp(static_cast<int>(idx + sizes.second), min, max);
    return NeighborHoodIndices(itmin, itmax + 1);
  }

  /// @brief Get #size bins which neighbor the one given
  ///
  /// This is the version for Bound
  ///
  /// @param [in] idx requested bin index
  /// @param [in] sizes how many neighboring bins (up/down)
  /// @return Set of neighboring bin indices (global)
  /// @note Bound varies given bin and allows 1 and NBins (regular bins)
  ///       as neighbors
  NeighborHoodIndices neighborHoodIndices(std::size_t idx,
                                          std::pair<int, int> sizes = {-1,
                                                                       1}) const
    requires(bdt == AxisBoundaryType::Bound)
  {
    if (idx <= 0 || idx >= (getNBins() + 1)) {
      return NeighborHoodIndices();
    }
    constexpr int min = 1;
    const int max = getNBins();
    const int itmin = std::clamp(static_cast<int>(idx) + sizes.first, min, max);
    const int itmax =
        std::clamp(static_cast<int>(idx) + sizes.second, min, max);
    return NeighborHoodIndices(itmin, itmax + 1);
  }

  /// @brief Get #size bins which neighbor the one given
  ///
  /// This is the version for Closed (i.e. Wrapping)
  ///
  /// @param [in] idx requested bin index
  /// @param [in] sizes how many neighboring bins (up/down)
  /// @return Set of neighboring bin indices (global)
  /// @note Closed varies given bin and allows bins on the opposite
  ///       side of the axis as neighbors. (excludes underflow / overflow)
  NeighborHoodIndices neighborHoodIndices(std::size_t idx,
                                          std::pair<int, int> sizes = {-1,
                                                                       1}) const
    requires(bdt == AxisBoundaryType::Closed)
  {
    // Handle invalid indices
    if (idx <= 0 || idx >= (getNBins() + 1)) {
      return NeighborHoodIndices();
    }

    // Handle corner case where user requests more neighbours than the number
    // of bins on the axis. All bins are returned in this case.

    const int max = getNBins();
    sizes.first = std::clamp(sizes.first, -max, max);
    sizes.second = std::clamp(sizes.second, -max, max);
    if (std::abs(sizes.first - sizes.second) >= max) {
      sizes.first = 1 - idx;
      sizes.second = max - idx;
    }

    // If the entire index range is not covered, we must wrap the range of
    // targeted neighbor indices into the range of valid bin indices. This may
    // split the range of neighbor indices in two parts:
    //
    // Before wraparound - [        XXXXX]XXX
    // After wraparound  - [ XXXX   XXXX ]
    //
    const int itmin = idx + sizes.first;
    const int itmax = idx + sizes.second;
    const std::size_t itfirst = wrapBin(itmin);
    const std::size_t itlast = wrapBin(itmax);
    if (itfirst <= itlast) {
      return NeighborHoodIndices(itfirst, itlast + 1);
    } else {
      return NeighborHoodIndices(itfirst, max + 1, 1, itlast + 1);
    }
  }

  /// @brief Converts bin index into a valid one for this axis.
  ///
  /// @note Open: bin index is clamped to [0, nBins+1]
  ///
  /// @param [in] bin The bin to wrap
  /// @return valid bin index
  std::size_t wrapBin(int bin) const
    requires(bdt == AxisBoundaryType::Open)
  {
    return std::max(std::min(bin, static_cast<int>(getNBins()) + 1), 0);
  }

  /// @brief Converts bin index into a valid one for this axis.
  ///
  /// @note Bound: bin index is clamped to [1, nBins]
  ///
  /// @param [in] bin The bin to wrap
  /// @return valid bin index
  std::size_t wrapBin(int bin) const
    requires(bdt == AxisBoundaryType::Bound)
  {
    return std::max(std::min(bin, static_cast<int>(getNBins())), 1);
  }

  /// @brief Converts bin index into a valid one for this axis.
  ///
  /// @note Closed: bin index wraps around to other side
  ///
  /// @param [in] bin The bin to wrap
  /// @return valid bin index
  std::size_t wrapBin(int bin) const
    requires(bdt == AxisBoundaryType::Closed)
  {
    const int w = getNBins();
    return 1 + (w + ((bin - 1) % w)) % w;
    // return int(bin<1)*w - int(bin>w)*w + bin;
  }

  /// @brief get corresponding bin index for given coordinate
  ///
  /// @param  [in] x input coordinate
  /// @return index of bin containing the given value
  ///
  /// @note Bin intervals are defined with closed lower bounds and open upper
  ///       bounds, that is \f$l <= x < u\f$ if the value @c x lies within a
  ///       bin with lower bound @c l and upper bound @c u.
  /// @note Bin indices start at @c 1. The underflow bin has the index @c 0
  ///       while the index <tt>nBins + 1</tt> indicates the overflow bin .
  std::size_t getBin(double x) const {
    return wrapBin(
        static_cast<int>(std::floor((x - getMin()) / getBinWidth()) + 1));
  }

  /// @brief get bin width
  ///
  /// @return constant width for all bins
  double getBinWidth(std::size_t /*bin*/ = 0) const { return m_width; }

  /// @brief get lower bound of bin
  ///
  /// @param  [in] bin index of bin
  /// @return lower bin boundary
  ///
  /// @pre @c bin must be a valid bin index (excluding the underflow bin),
  ///      i.e. \f$1 \le \text{bin} \le \text{nBins} + 1\f$
  ///
  /// @note Bin intervals have a closed lower bound, i.e. the lower boundary
  ///       belongs to the bin with the given bin index.
  double getBinLowerBound(std::size_t bin) const {
    return getMin() + (bin - 1) * getBinWidth();
  }

  /// @brief get upper bound of bin
  ///
  /// @param  [in] bin index of bin
  /// @return upper bin boundary
  ///
  /// @pre @c bin must be a valid bin index (excluding the overflow bin),
  ///      i.e. \f$0 \le \text{bin} \le \text{nBins}\f$
  ///
  /// @note Bin intervals have an open upper bound, i.e. the upper boundary
  ///       does @b not belong to the bin with the given bin index.
  double getBinUpperBound(std::size_t bin) const {
    return getMin() + bin * getBinWidth();
  }

  /// @brief get bin center
  ///
  /// @param  [in] bin index of bin
  /// @return bin center position
  ///
  /// @pre @c bin must be a valid bin index (excluding under-/overflow bins),
  ///      i.e. \f$1 \le \text{bin} \le \text{nBins}\f$
  double getBinCenter(std::size_t bin) const {
    return getMin() + (bin - 0.5) * getBinWidth();
  }

  /// @brief get maximum of binning range
  ///
  /// @return maximum of binning range
  double getMax() const override { return m_max; }

  /// @brief get minimum of binning range
  ///
  /// @return minimum of binning range
  double getMin() const override { return m_min; }

  /// @brief get total number of bins
  ///
  /// @return total number of bins (excluding under-/overflow bins)
  std::size_t getNBins() const override { return m_bins; }

  /// @brief check whether value is inside axis limits
  ///
  /// @return @c true if \f$\text{xmin} \le x < \text{xmax}\f$, otherwise
  ///         @c false
  ///
  /// @post If @c true is returned, the bin containing the given value is a
  ///       valid bin, i.e. it is neither the underflow nor the overflow bin.
  bool isInside(double x) const { return (m_min <= x) && (x < m_max); }

  /// @brief Return a vector of bin edges
  /// @return Vector which contains the bin edges
  std::vector<double> getBinEdges() const override {
    std::vector<double> binEdges;
    for (std::size_t i = 1; i <= m_bins; i++) {
      binEdges.push_back(getBinLowerBound(i));
    }
    binEdges.push_back(getBinUpperBound(m_bins));
    return binEdges;
  }

  friend std::ostream& operator<<(std::ostream& os, const Axis& axis) {
    os << "Axis<Equidistant, " << bdt << ">(";
    os << axis.m_min << ", ";
    os << axis.m_max << ", ";
    os << axis.m_bins << ")";
    return os;
  }

 protected:
  void toStream(std::ostream& os) const override { os << *this; }

 private:
  /// minimum of binning range
  double m_min{};
  /// maximum of binning range
  double m_max{};
  /// constant bin width
  double m_width{};
  /// number of bins (excluding under-/overflow bins)
  std::size_t m_bins{};
};

/// @brief calculate bin indices for a variable binning
///
/// This class provides some basic functionality for calculating bin indices
/// for a given binning with variable bin sizes.
template <AxisBoundaryType bdt>
class Axis<AxisType::Variable, bdt> : public IAxis {
 public:
  static constexpr AxisType type = AxisType::Variable;

  /// @param [in] binEdges vector of bin edges
  /// @pre @c binEdges must be strictly sorted in ascending order.
  /// @pre @c binEdges must contain at least two entries.
  ///
  /// Create a binning structure with @c nBins variable-sized bins from the
  /// given bin boundaries. @c nBins is given by the number of bin edges
  /// reduced by one.
  explicit Axis(std::vector<double> binEdges)
      : m_binEdges(std::move(binEdges)) {}

  /// @param [in] typeTag boundary type tag
  /// @param [in] binEdges vector of bin edges
  /// @pre @c binEdges must be strictly sorted in ascending order.
  /// @pre @c binEdges must contain at least two entries.
  ///
  /// Create a binning structure with @c nBins variable-sized bins from the
  /// given bin boundaries. @c nBins is given by the number of bin edges
  /// reduced by one.
  Axis(AxisBoundaryTypeTag<bdt> typeTag, std::vector<double> binEdges)
      : Axis(std::move(binEdges)) {
    static_cast<void>(typeTag);
  }

  /// @brief returns whether the axis is equidistante
  ///
  /// @return bool is equidistant
  bool isEquidistant() const override { return false; }

  /// @brief returns whether the axis is variable
  ///
  /// @return bool is variable
  bool isVariable() const override { return true; }

  /// @brief returns the type of the axis
  /// @return @c AxisType of this axis
  AxisType getType() const override { return type; }

  /// @brief returns the boundary type set in the template param
  ///
  /// @return @c AxisBoundaryType of this axis
  AxisBoundaryType getBoundaryType() const override { return bdt; }

  /// @brief Get #size bins which neighbor the one given
  ///
  /// Generic overload with symmetric size
  ///
  /// @param [in] idx requested bin index
  /// @param [in] size how many neighboring bins
  /// @return Set of neighboring bin indices (global)
  NeighborHoodIndices neighborHoodIndices(std::size_t idx,
                                          std::size_t size = 1) const {
    return neighborHoodIndices(idx,
                               std::make_pair(-static_cast<int>(size), size));
  }

  /// @brief Get #size bins which neighbor the one given
  ///
  /// This is the version for Open
  ///
  /// @param [in] idx requested bin index
  /// @param [in] sizes how many neighboring bins (up/down)
  /// @return Set of neighboring bin indices (global)
  /// @note Open varies given bin and allows 0 and NBins+1 (underflow,
  /// overflow)
  ///       as neighbors
  NeighborHoodIndices neighborHoodIndices(std::size_t idx,
                                          std::pair<int, int> sizes = {-1,
                                                                       1}) const
    requires(bdt == AxisBoundaryType::Open)
  {
    constexpr int min = 0;
    const int max = getNBins() + 1;
    const int itmin = std::max(min, static_cast<int>(idx) + sizes.first);
    const int itmax = std::min(max, static_cast<int>(idx) + sizes.second);
    return NeighborHoodIndices(itmin, itmax + 1);
  }

  /// @brief Get #size bins which neighbor the one given
  ///
  /// This is the version for Bound
  ///
  /// @param [in] idx requested bin index
  /// @param [in] sizes how many neighboring bins (up/down)
  /// @return Set of neighboring bin indices (global)
  /// @note Bound varies given bin and allows 1 and NBins (regular bins)
  ///       as neighbors
  NeighborHoodIndices neighborHoodIndices(std::size_t idx,
                                          std::pair<int, int> sizes = {-1,
                                                                       1}) const
    requires(bdt == AxisBoundaryType::Bound)
  {
    if (idx <= 0 || idx >= (getNBins() + 1)) {
      return NeighborHoodIndices();
    }
    constexpr int min = 1;
    const int max = getNBins();
    const int itmin = std::max(min, static_cast<int>(idx) + sizes.first);
    const int itmax = std::min(max, static_cast<int>(idx) + sizes.second);
    return NeighborHoodIndices(itmin, itmax + 1);
  }

  /// @brief Get #size bins which neighbor the one given
  ///
  /// This is the version for Closed
  ///
  /// @param [in] idx requested bin index
  /// @param [in] sizes how many neighboring bins (up/down)
  /// @return Set of neighboring bin indices (global)
  /// @note Closed varies given bin and allows bins on the opposite
  ///       side of the axis as neighbors. (excludes underflow / overflow)
  NeighborHoodIndices neighborHoodIndices(std::size_t idx,
                                          std::pair<int, int> sizes = {-1,
                                                                       1}) const
    requires(bdt == AxisBoundaryType::Closed)
  {
    // Handle invalid indices
    if (idx <= 0 || idx >= (getNBins() + 1)) {
      return NeighborHoodIndices();
    }

    // Handle corner case where user requests more neighbours than the number
    // of bins on the axis. All bins are returned in this case

    const int max = getNBins();
    sizes.first = std::clamp(sizes.first, -max, max);
    sizes.second = std::clamp(sizes.second, -max, max);
    if (std::abs(sizes.first - sizes.second) >= max) {
      sizes.first = 1 - idx;
      sizes.second = max - idx;
    }

    // If the entire index range is not covered, we must wrap the range of
    // targeted neighbor indices into the range of valid bin indices. This may
    // split the range of neighbor indices in two parts:
    //
    // Before wraparound - [        XXXXX]XXX
    // After wraparound  - [ XXXX   XXXX ]
    //
    const int itmin = idx + sizes.first;
    const int itmax = idx + sizes.second;
    const std::size_t itfirst = wrapBin(itmin);
    const std::size_t itlast = wrapBin(itmax);
    if (itfirst <= itlast) {
      return NeighborHoodIndices(itfirst, itlast + 1);
    } else {
      return NeighborHoodIndices(itfirst, max + 1, 1, itlast + 1);
    }
  }

  /// @brief Converts bin index into a valid one for this axis.
  ///
  /// @note Open: bin index is clamped to [0, nBins+1]
  ///
  /// @param [in] bin The bin to wrap
  /// @return valid bin index
  std::size_t wrapBin(int bin) const
    requires(bdt == AxisBoundaryType::Open)
  {
    return std::max(std::min(bin, static_cast<int>(getNBins()) + 1), 0);
  }

  /// @brief Converts bin index into a valid one for this axis.
  ///
  /// @note Bound: bin index is clamped to [1, nBins]
  ///
  /// @param [in] bin The bin to wrap
  /// @return valid bin index
  std::size_t wrapBin(int bin) const
    requires(bdt == AxisBoundaryType::Bound)
  {
    return std::max(std::min(bin, static_cast<int>(getNBins())), 1);
  }

  /// @brief Converts bin index into a valid one for this axis.
  ///
  /// @note Closed: bin index wraps around to other side
  ///
  /// @param [in] bin The bin to wrap
  /// @return valid bin index
  std::size_t wrapBin(int bin) const
    requires(bdt == AxisBoundaryType::Closed)
  {
    const int w = getNBins();
    return 1 + (w + ((bin - 1) % w)) % w;
    // return int(bin<1)*w - int(bin>w)*w + bin;
  }

  /// @brief get corresponding bin index for given coordinate
  ///
  /// @param  [in] x input coordinate
  /// @return index of bin containing the given value
  ///
  /// @note Bin intervals are defined with closed lower bounds and open upper
  ///       bounds, that is \f$l <= x < u\f$ if the value @c x lies within a
  ///       bin with lower bound @c l and upper bound @c u.
  /// @note Bin indices start at @c 1. The underflow bin has the index @c 0
  ///       while the index <tt>nBins + 1</tt> indicates the overflow bin .
  std::size_t getBin(double x) const {
    const auto it =
        std::upper_bound(std::begin(m_binEdges), std::end(m_binEdges), x);
    return wrapBin(std::distance(std::begin(m_binEdges), it));
  }

  /// @brief get bin width
  ///
  /// @param  [in] bin index of bin
  /// @return width of given bin
  ///
  /// @pre @c bin must be a valid bin index (excluding under-/overflow bins),
  ///      i.e. \f$1 \le \text{bin} \le \text{nBins}\f$
  double getBinWidth(std::size_t bin) const {
    return m_binEdges.at(bin) - m_binEdges.at(bin - 1);
  }

  /// @brief get lower bound of bin
  ///
  /// @param  [in] bin index of bin
  /// @return lower bin boundary
  ///
  /// @pre @c bin must be a valid bin index (excluding the underflow bin),
  ///      i.e. \f$1 \le \text{bin} \le \text{nBins} + 1\f$
  ///
  /// @note Bin intervals have a closed lower bound, i.e. the lower boundary
  ///       belongs to the bin with the given bin index.
  double getBinLowerBound(std::size_t bin) const {
    return m_binEdges.at(bin - 1);
  }

  /// @brief get upper bound of bin
  ///
  /// @param  [in] bin index of bin
  /// @return upper bin boundary
  ///
  /// @pre @c bin must be a valid bin index (excluding the overflow bin),
  ///      i.e. \f$0 \le \text{bin} \le \text{nBins}\f$
  ///
  /// @note Bin intervals have an open upper bound, i.e. the upper boundary
  ///       does @b not belong to the bin with the given bin index.
  double getBinUpperBound(std::size_t bin) const { return m_binEdges.at(bin); }

  /// @brief get bin center
  ///
  /// @param  [in] bin index of bin
  /// @return bin center position
  ///
  /// @pre @c bin must be a valid bin index (excluding under-/overflow bins),
  ///      i.e. \f$1 \le \text{bin} \le \text{nBins}\f$
  double getBinCenter(std::size_t bin) const {
    return 0.5 * (getBinLowerBound(bin) + getBinUpperBound(bin));
  }

  /// @brief get maximum of binning range
  ///
  /// @return maximum of binning range
  double getMax() const override { return m_binEdges.back(); }

  /// @brief get minimum of binning range
  ///
  /// @return minimum of binning range
  double getMin() const override { return m_binEdges.front(); }

  /// @brief get total number of bins
  ///
  /// @return total number of bins (excluding under-/overflow bins)
  std::size_t getNBins() const override { return m_binEdges.size() - 1; }

  /// @brief check whether value is inside axis limits
  ///
  /// @return @c true if \f$\text{xmin} \le x < \text{xmax}\f$, otherwise
  ///         @c false
  ///
  /// @post If @c true is returned, the bin containing the given value is a
  ///       valid bin, i.e. it is neither the underflow nor the overflow bin.
  bool isInside(double x) const {
    return (m_binEdges.front() <= x) && (x < m_binEdges.back());
  }

  /// @brief Return a vector of bin edges
  /// @return Vector which contains the bin edges
  std::vector<double> getBinEdges() const override { return m_binEdges; }

  friend std::ostream& operator<<(std::ostream& os, const Axis& axis) {
    os << "Axis<Variable, " << bdt << ">(";
    os << axis.m_binEdges.front();
    for (std::size_t i = 1; i < axis.m_binEdges.size(); i++) {
      os << ", " << axis.m_binEdges[i];
    }
    os << ")";
    return os;
  }

 protected:
  void toStream(std::ostream& os) const override { os << *this; }

 private:
  /// vector of bin edges (sorted in ascending order)
  std::vector<double> m_binEdges;
};
}  // namespace Acts
