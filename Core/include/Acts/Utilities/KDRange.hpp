// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <limits>

namespace Acts {

/// @brief A one-dimensional range between two points.
///
/// This type describes a one-demensional range of values, designed to be used
/// in the construction of more complex multi-dimensional types. The class
/// provides functionality for growing and shrinking the range, as well as some
/// other utilities.
///
/// @tparam Type The scalar type of the values contained in this range.
template <typename Type>
class OneDRange {
 public:
  /// @brief Construct a new degenerate range object
  ///
  /// This constructor coonstructs a degenerate range object with a maximum
  /// lower than the minimum. In other words, this range is empty.
  OneDRange()
      : m_min(std::numeric_limits<Type>::lowest()),
        m_max(std::numeric_limits<Type>::max()) {}

  /// @brief Construct a new range object from a lower and upper bound
  ///
  /// Construct a new range object given the values for the minimum and
  /// maximum. Note that it is perfectly possible to construct a degenerate
  /// range in this way.
  ///
  /// @param min The minimum value in the range (inclusive)
  /// @param max The maximum value in the range (inclusive)
  OneDRange(Type min, Type max) : m_min(min), m_max(max) {}

  /// @brief Construct a new range object from a pair of bounds
  ///
  /// Construct a new range object from a pair of values, the first of which
  /// is taken to be the minimum, and the second of which is taken to be the
  /// maximum.
  ///
  /// @param p The pair of values to use as the minimum and maximum
  OneDRange(const std::pair<Type, Type>& p) : m_min(p.first), m_max(p.second) {}

  /// @brief Set the minimum value
  ///
  /// Override the minimum value of the range, regardless of what was already
  /// set.
  ///
  /// @note If you want to shrink or expand the range, use the shrink and
  /// expand methods.
  ///
  /// @param v The value to use as the new minimum
  void set_min(const Type& v) { m_min = v; }

  /// @brief Set the maximum value
  ///
  /// Override the maximum value of the range, regardless of what was already
  /// set.
  ///
  /// @note If you want to shrink or expand the range, use the shrink and
  /// expand methods.
  ///
  /// @param v The value to use as the new maximum
  void set_max(const Type& v) { m_max = v; }

  /// @brief Shrink a range by increasing the minimum value
  ///
  /// Shrink the range by increasing the minimum value. If the given value is
  /// smaller than the current minimum (in other words, if the proposed new
  /// range would be larger than the current range), this is a no-op.
  ///
  /// @param v The proposed new minimum for the range
  void shrink_min(const Type& v) { m_min = std::max(m_min, v); }

  /// @brief Shrink a range by decreasing the maximum value
  ///
  /// Shrink the range by decreasing the maximum value. If the given value is
  /// larger than the current maximum (in other words, if the proposed new
  /// range would be larger than the current range), this is a no-op.
  ///
  /// @param v The proposed new maximum for the range
  void shrink_max(const Type& v) { m_max = std::min(m_max, v); }

  /// @brief Expand a range by decreasing the minimum value
  ///
  /// Expand the range by decreasing the minimum value. If the given value is
  /// larger than the current minimum (in other words, if the proposed new
  /// range would be smaller than the current range), this is a no-op.
  ///
  /// @param v The proposed new minimum for the range
  void expand_min(const Type& v) { m_min = std::min(m_min, v); }

  /// @brief Expand a range by increasing the maximum value
  ///
  /// Expand the range by increasing the maximum value. If the given value is
  /// smaller than the current maximum (in other words, if the proposed new
  /// range would be smaller than the current range), this is a no-op.
  ///
  /// @param v The proposed new maximum for the range
  void expand_max(const Type& v) { m_max = std::max(m_max, v); }

  /// @brief Return the minimum value of the range (inclusive)
  Type min(void) const { return m_min; }

  /// @brief Return the maximum value of the range (inclusive)
  Type max(void) const { return m_max; }

  /// @brief Compute the size of the range
  ///
  /// The size of a range is defined as the difference between the minimum
  /// and the maximum. For degenerate ranges, this is zero.
  ///
  /// @return The size of the range
  Type size(void) const {
    return std::min(static_cast<Type>(0), m_max - m_min);
  }

  /// @brief Determine if this range is degenerate or not
  ///
  /// A degenerate range has a minimum higher than the maximum, and thus
  /// cannot contain any values.
  ///
  /// @return true The range is degenerate and has size zero
  /// @return false The range is not degenerate
  bool degenerate(void) const { return m_min >= m_max; }

  /// @brief Determine if the range contains a given value
  ///
  /// A value is inside a range if and only if it is greater than the minimum
  /// and smaller than the maximum.
  ///
  /// @param v The value to check
  ///
  /// @return true The value is inside the range
  /// @return false The value is not inside the range
  bool contains(const Type& v) const { return m_min <= v && v <= m_max; }

  /// @brief Determine whether the range intersects another range
  ///
  /// The intersection of a range is the space where both ranges overlap. If
  /// the ranges overlap at all, they are said to intersect. This operation
  /// is commutative.
  ///
  /// @param o The other range to check
  ///
  /// @return true The ranges intersect
  /// @return false The ranges do not intersect
  bool operator&&(const OneDRange<Type>& o) const {
    return m_min <= o.max() && o.min() <= m_max;
  }

  /// @brief Determine whether the range is equal to another range
  ///
  /// Two ranges are equal if and only if their minima and maxima are the
  /// same.
  ///
  /// @warning This method relies on the existence of a well-defined notion
  /// of equality for the underlying types. Using this method on floating
  /// ranges may have unintended effecrs.
  ///
  /// @param o The other range to check
  ///
  /// @return true The ranges are equal
  /// @return false The ranges are not equal
  bool operator==(const OneDRange<Type>& o) const {
    return min() == o.min() && max() == o.max();
  }

  /// @brief Determine whether the left-hand range is a subset of the
  /// right-hand range
  ///
  /// A range is a subset of another range if and only if all values
  /// contained in the first range are also contained in the second range.
  ///
  /// @param o The other range to check
  ///
  /// @return true The left-hand range is a subset of the right-hand range
  /// @return false The left-hand range is not a subset of the right-hand
  /// range
  bool operator<=(const OneDRange<Type>& o) const {
    return min() >= o.min() && max() <= o.max();
  }

  /// @brief Determine whether the left-hand range is a superset of the
  /// right-hand range
  ///
  /// A range is a superset of another range if and only if all values
  /// contained in the second range are also contained in the first range.
  ///
  /// @param o The other range to check
  ///
  /// @return true The left-hand range is a superset of thr right-hand range
  /// @return false The left-hand range is not a superset of the right-hand
  /// range
  bool operator>=(const OneDRange<Type>& o) const {
    return min() <= o.min() && max() >= o.max();
  }

  /// @brief Assignment operator
  ///
  /// Copy the right-hand range into the left-hand range, which means setting
  /// the minimum and maximum to equal the minimum and maximum of the
  /// right-hand side.
  ///
  /// @param o The range of values to copy
  ///
  /// @return This range
  OneDRange<Type>& operator=(const OneDRange<Type>& o) {
    m_min = o.min();
    m_max = o.max();

    return *this;
  }

  /// @brief Compute the intersection of two ranges
  ///
  /// The intersection of two ranges is the range containing all values
  /// contained in both ranges. If the two ranges do not intersect, the
  /// intersection is a degenerate range. This operation is commutative.
  ///
  /// @param o The range to compute the intersection with
  ///
  /// @return The intersection range between the two ranges
  OneDRange<Type> operator&(const OneDRange<Type>& o) const {
    return OneDRange<Type>(std::max(m_min, o.min()), std::min(m_max, o.max()));
  }

 private:
  Type m_min, m_max;
};

/// @brief An orthogonal range in an arbitrary number of dimensions
///
/// By combining a number one-dimensional ranges we can (under the assumption
/// that our axes are orthogonal) construct an orthogonal range of values. In
/// other words, a hyperrectangular volume in space.
///
/// @tparam Dims The number of dimensions in our range
/// @tparam Type The scalar type of our ranges
/// @tparam Vector The vector type used to define coordinates
template <std::size_t Dims, typename Type,
          template <typename, std::size_t> typename Vector = std::array>
class KDRange {
 public:
  /// @brief The type used to describe coordinates in our range
  using coordinate_t = Vector<Type, Dims>;

  /// @brief Determine whether this range is degenerate
  ///
  /// A degenerate multi-dimensional range has no volume and cannot contain any
  /// values. This is the case if any of its dimensions are degenerate.
  ///
  /// @return true The range is degenerate
  /// @return false The range is not degenerate
  bool degenerate(void) const {
    for (std::size_t i = 0; i < Dims; ++i) {
      if (m_dims[i].degenerate()) {
        return true;
      }
    }
    return false;
  }

  /// @brief Determine whether the range contains a certain point
  ///
  /// This is true if and only if the range contains the point in all of its
  /// dimensions.
  ///
  /// @param v The coordinate to check for membership in the range
  ///
  /// @return true The coordinate is inside the range
  /// @return false The coordinate is outside the range
  bool contains(const coordinate_t& v) const {
    for (std::size_t i = 0; i < Dims; ++i) {
      if (!m_dims[i].contains(v[i])) {
        return false;
      }
    }

    return true;
  }

  /// @brief Access one of the dimensional ranges of the volume
  ///
  /// @param i The index of the dimension to access
  /// @return A reference to the dimension contained in this range
  OneDRange<Type>& operator[](const std::size_t& i) {
    return m_dims[static_cast<std::size_t>(i)];
  }

  /// @brief Access one of the dimensional ranges of the volume
  ///
  /// @param i The index of the dimension to access
  /// @return A reference to the dimension contained in this range
  const OneDRange<Type>& operator[](const std::size_t& i) const {
    return m_dims[static_cast<std::size_t>(i)];
  }

  /// @brief Determine whether two ranges are equal
  ///
  /// Two n-dimensional ranges are equal if and only if they are equal in each
  /// of their n dimensions.
  ///
  /// @param o The other range to check for equality
  ///
  /// @return true The ranges are equal
  /// @return false The ranges are not equal
  bool operator==(const KDRange<Dims, Type, Vector>& o) const {
    for (std::size_t i = 0; i < Dims; ++i) {
      if (!(this->m_dims[i] == o[i])) {
        return false;
      }
    }

    return true;
  }

  /// @brief Determine whether one range is a subset of another range
  ///
  /// One range is a subset of another range if and only if all points
  /// contained within the first set are also contained within the second set.
  /// Alternatively, this is equivalent to each of the first range's
  /// one-dimensional ranges being a subset of the second range's equivalent
  /// one-dimensional range.
  ///
  /// @param o The other range to compare to
  ///
  /// @return true The first range is a subset of the second range
  /// @return false The first range is not a subset of the second range
  bool operator<=(const KDRange<Dims, Type, Vector>& o) const {
    for (std::size_t i = 0; i < Dims; ++i) {
      if (!(this->m_dims[i] <= o[i])) {
        return false;
      }
    }

    return true;
  }

  /// @brief Determine whether one range is a superset of another range
  ///
  /// One range is a superset of another range if and only if all points
  /// contained within the second range are also contained within the first
  /// range. Alternatively, this is equivalent to each of the one-dimensional
  /// ranges in the first range being a superset of the corresponding
  /// one-dimensional range in the second range.
  ///
  /// @param o The other range to compare to
  ///
  /// @return true The left-hand range is a superset of the right-hand range
  /// @return false The left-hand range is not a superset of the right-hand
  /// range
  bool operator>=(const KDRange<Dims, Type, Vector>& o) const {
    for (std::size_t i = 0; i < Dims; ++i) {
      if (!(this->m_dims[i] >= o[i])) {
        return false;
      }
    }

    return true;
  }

  /// @brief Compute the intersection of this range with another range
  ///
  /// The intersection of one orthogonal range with another orthogonal range is
  /// in itself an orthogonal range. This operation is commutative. This
  /// intersection between two n-dimensional ranges is defined simply as the
  /// intersection in each dimension of the two ranges.
  ///
  /// @param o The orthogonal range to compute the intersection with
  ///
  /// @return The intersection between the ranges
  KDRange<Dims, Type, Vector> operator&(
      const KDRange<Dims, Type, Vector>& o) const {
    KDRange<Dims, Type> res;

    for (std::size_t i = 0; i < Dims; ++i) {
      res[i] = m_dims[i] & o[i];
    }

    return res;
  }

  /// @brief Determine whether this range intersects another
  ///
  /// Two n-dimensional ranges intersect if and only if they intersect in every
  /// one of their n dimensions. Otherwise, they are disjoint.
  ///
  /// @param r The other range to check
  ///
  /// @return true The ranges intersect
  /// @return false The ranges do not intersect
  bool operator&&(const KDRange<Dims, Type, Vector>& r) const {
    for (std::size_t i = 0; i < Dims; ++i) {
      if (!(m_dims[i] && r[i])) {
        return false;
      }
    }

    return true;
  }

  std::string str(void) const {
    std::stringstream s;

    for (std::size_t i = 0; i < Dims; ++i) {
      s << m_dims[i].min() << " <= v[" << i << "] <= " << m_dims[i].max();
      if (i != Dims - 1) {
        s << ", ";
      }
    }
    return s.str();
  }

  // KDRange<Dims, Type

 private:
  std::array<OneDRange<Type>, Dims> m_dims;
};
}  // namespace Acts
