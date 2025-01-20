// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <limits>
#include <sstream>
#include <string>

namespace Acts {
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
class RangeXD {
 private:
  // @TODO: Replace with std::span or boost::span once available
  template <typename, std::size_t>
  struct SingleElementContainer {
    Type* element;

    Type& operator[](std::size_t i) {
      (void)i;
      assert(i == 0);

      return *element;
    }
  };

 public:
  RangeXD() {
    for (std::size_t i = 0; i < Dims; ++i) {
      min(i) = std::numeric_limits<Type>::lowest();
      max(i) = std::numeric_limits<Type>::max();
    }
  }

  /// @brief Construct a range from a pair of minimum and maximum values
  /// @param minima The minimum values of the range
  /// @param maxima The maximum values of the range
  RangeXD(Vector<Type, Dims> minima, Vector<Type, Dims> maxima)
      : m_minima(minima), m_maxima(maxima) {}

  /// @brief Construct a range from a pair of single minimum and maximum values
  /// @note Only available for one-dimensional ranges
  /// @param minimum The minimum value of the range
  /// @param maximum The maximum value of the range
  RangeXD(Type minimum, Type maximum)
    requires(Dims == 1)
      : m_minima({minimum}), m_maxima({maximum}) {}

  /// @brief Construct a range from a pair of minimum and maximum values
  /// @note Only available for one-dimensional ranges
  /// @param p The pair of minimum and maximum values
  explicit RangeXD(const std::pair<Type, Type>& p)
    requires(Dims == 1)
      : m_minima({p.first}), m_maxima({p.second}) {}

  /// @brief Determine whether this range is degenerate
  ///
  /// A degenerate multi-dimensional range has no volume and cannot contain
  /// any values. This is the case if any of its dimensions are degenerate.
  ///
  /// @return true The range is degenerate
  /// @return false The range is not degenerate
  bool degenerate() const {
    for (std::size_t i = 0; i < Dims; ++i) {
      if (min(i) >= max(i)) {
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
  template <template <typename, std::size_t> typename coordinate_t = std::array>
  bool contains(const coordinate_t<Type, Dims>& v) const {
    for (std::size_t i = 0; i < Dims; ++i) {
      if (!(min(i) <= v[i] && v[i] < max(i))) {
        return false;
      }
    }

    return true;
  }

  /// @brief Access one of the dimensional ranges of the volume
  ///
  /// @param i The index of the dimension to access
  /// @return A reference to the dimension contained in this range
  RangeXD<1, Type, SingleElementContainer> operator[](const std::size_t& i) {
    return RangeXD<1, Type, SingleElementContainer>{{&min(i)}, {&max(i)}};
  }

  /// @brief Access one of the dimensional ranges of the volume
  ///
  /// @param i The index of the dimension to access
  /// @return A reference to the dimension contained in this range
  RangeXD<1, Type> operator[](const std::size_t& i) const {
    return RangeXD<1, Type>{{min(i)}, {max(i)}};
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
  template <template <typename, std::size_t> typename V>
  RangeXD& operator=(const RangeXD<Dims, Type, V>& o) {
    for (std::size_t i = 0; i < Dims; ++i) {
      min(i) = o.min(i);
      max(i) = o.max(i);
    }

    return *this;
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
  bool operator==(const RangeXD<Dims, Type, Vector>& o) const {
    for (std::size_t i = 0; i < Dims; ++i) {
      if (!(min(i) == o.min(i) && max(i) == o.max(i))) {
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
  bool operator<=(const RangeXD<Dims, Type, Vector>& o) const {
    for (std::size_t i = 0; i < Dims; ++i) {
      if (!(min(i) >= o.min(i) && max(i) <= o.max(i))) {
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
  bool operator>=(const RangeXD<Dims, Type, Vector>& o) const {
    for (std::size_t i = 0; i < Dims; ++i) {
      if (!(min(i) <= o.min(i) && max(i) >= o.max(i))) {
        return false;
      }
    }

    return true;
  }

  /// @brief Compute the intersection of this range with another range
  ///
  /// The intersection of one orthogonal range with another orthogonal range
  /// is in itself an orthogonal range. This operation is commutative. This
  /// intersection between two n-dimensional ranges is defined simply as the
  /// intersection in each dimension of the two ranges.
  ///
  /// @param o The orthogonal range to compute the intersection with
  ///
  /// @return The intersection between the ranges
  RangeXD<Dims, Type, Vector> operator&(
      const RangeXD<Dims, Type, Vector>& o) const {
    RangeXD<Dims, Type> res;

    for (std::size_t i = 0; i < Dims; ++i) {
      res.min(i) = std::max(min(i), o.min(i));
      res.max(i) = std::min(max(i), o.max(i));
    }

    return res;
  }

  /// @brief Update the range to the intersection with another range
  ///
  /// This is the assignment version of the operator& method, meaning that it
  /// updates the object on which it is called rather than producing a new
  /// range.
  ///
  /// @param o The range to compute the intersection with
  ///
  /// @return This object
  RangeXD<Dims, Type, Vector>& operator&=(
      const RangeXD<Dims, Type, Vector>& o) {
    for (std::size_t i = 0; i < Dims; ++i) {
      min(i) = std::max(min(i), o.min(i));
      max(i) = std::min(max(i), o.max(i));
    }

    return *this;
  }

  /// @brief Determine whether this range intersects another
  ///
  /// Two n-dimensional ranges intersect if and only if they intersect in
  /// every one of their n dimensions. Otherwise, they are disjoint.
  ///
  /// @param r The other range to check
  ///
  /// @return true The ranges intersect
  /// @return false The ranges do not intersect
  bool operator&&(const RangeXD<Dims, Type, Vector>& r) const {
    for (std::size_t i = 0; i < Dims; ++i) {
      if (!(min(i) < r.max(i) && r.min(i) < max(i))) {
        return false;
      }
    }

    return true;
  }

  /// @brief Represent the range as a string
  ///
  /// This method produces a helpful string that can be used to debug the
  /// range if needed. Not really designed to be used in production code.
  ///
  /// @return A string representing the range
  std::string toString(void) const {
    std::stringstream s;

    for (std::size_t i = 0; i < Dims; ++i) {
      s << min(i) << " <= v[" << i << "] <= " << max(i);
      if (i != Dims - 1) {
        s << ", ";
      }
    }

    return s.str();
  }

  /// @brief Shrink a range by increasing the minimum value
  ///
  /// Shrink the range by increasing the minimum value. If the given value is
  /// smaller than the current minimum (in other words, if the proposed new
  /// range would be larger than the current range), this is a no-op.
  ///
  /// @param i The index of the dimension to shrink
  /// @param v The proposed new minimum for the range
  void shrinkMin(std::size_t i, const Type& v) { min(i) = std::max(min(i), v); }

  /// @brief Shrink a range by decreasing the maximum value
  ///
  /// Shrink the range by decreasing the maximum value. If the given value is
  /// larger than the current maximum (in other words, if the proposed new
  /// range would be larger than the current range), this is a no-op.
  ///
  /// @param i The index of the dimension to shrink
  /// @param v The proposed new maximum for the range
  void shrinkMax(std::size_t i, const Type& v) { max(i) = std::min(max(i), v); }

  /// @brief Shrink a range on both ends
  ///
  /// Shrink a range by increasing the minimum value as well as decreasing the
  /// maximum value. If either of the values are already smaller or larger
  /// (respectively) than the proposed values, then that particular boundary
  /// of the interval is not shrunk.
  ///
  /// @note After this operation, the range is always equal to or smaller than
  /// [min, max].
  ///
  /// @param i The index of the dimension to shrink
  /// @param min The proposed new minimum for the range
  /// @param max The proposed new maximum for the range
  void shrink(std::size_t i, const Type& min, const Type& max) {
    shrinkMin(i, min);
    shrinkMax(i, max);
  }

  /// @brief Expand a range by decreasing the minimum value
  ///
  /// Expand the range by decreasing the minimum value. If the given value is
  /// larger than the current minimum (in other words, if the proposed new
  /// range would be smaller than the current range), this is a no-op.
  ///
  /// @param i The index of the dimension to expand
  /// @param v The proposed new minimum for the range
  void expandMin(std::size_t i, const Type& v) { min(i) = std::min(min(i), v); }

  /// @brief Expand a range by increasing the maximum value
  ///
  /// Expand the range by increasing the maximum value. If the given value is
  /// smaller than the current maximum (in other words, if the proposed new
  /// range would be smaller than the current range), this is a no-op.
  ///
  /// @param i The index of the dimension to expand
  /// @param v The proposed new maximum for the range
  void expandMax(std::size_t i, const Type& v) { max(i) = std::max(max(i), v); }

  /// @brief Expand a range on both ends
  ///
  /// Expand a range by decreasing the minimum value as well as increasing the
  /// maximum value. If either of the values are already larger or smaller
  /// (respectively) than the proposed values, then that particular boundary
  /// of the interval is not expanded.
  ///
  /// @note After this operation, the range is always equal to or larger than
  /// [min, max].
  ///
  /// @param i The index of the dimension to expand
  /// @param min The proposed new minimum for the range
  /// @param max The proposed new maximum for the range
  void expand(std::size_t i, const Type& min, const Type& max) {
    expandMin(i, min);
    expandMax(i, max);
  }

  /// @brief Set the minimum value
  ///
  /// Override the minimum value of the range, regardless of what was already
  /// set.
  ///
  /// @note If you want to shrink or expand the range, use the shrink and
  /// expand methods.
  ///
  /// @param i The index of the dimension to set
  /// @param v The value to use as the new minimum
  void setMin(std::size_t i, const Type& v) { min(i) = v; }

  /// @brief Set the maximum value
  ///
  /// Override the maximum value of the range, regardless of what was already
  /// set.
  ///
  /// @note If you want to shrink or expand the range, use the shrink and
  /// expand methods.
  ///
  /// @param i The index of the dimension to set
  /// @param v The value to use as the new maximum
  void setMax(std::size_t i, const Type& v) { max(i) = v; }

  /// @brief Set the minimum and maximum value
  ///
  /// Override both the minimum and maximum value of the range, regardless of
  /// what they were set to.
  ///
  /// @note If you want to shrink or expand the range, use the shrink and
  /// expand methods.
  ///
  /// @note After this operation, the range should be exactly equal to [min,
  /// max]
  ///
  /// @param i The index of the dimension to set
  /// @param min The new minimum value of the range
  /// @param max The new maximum value of the range
  void set(std::size_t i, const Type& min, const Type& max) {
    setMin(i, min);
    setMax(i, max);
  }

  /// @brief Return the minimum value of the range @p i (inclusive)
  /// @param i The index of the dimension to access
  Type& min(std::size_t i) { return m_minima[i]; }

  /// @brief Return the maximum value of the range @p i (inclusive)
  /// @param i The index of the dimension to access
  Type& max(std::size_t i) { return m_maxima[i]; }

  /// @brief Return the minimum value of the range @p i (inclusive)
  /// @param i The index of the dimension to access
  Type min(std::size_t i) const { return m_minima[i]; }

  /// @brief Return the maximum value of the range @p i (inclusive)
  /// @param i The index of the dimension to access
  Type max(std::size_t i) const { return m_maxima[i]; }

  /// Methods for manipulating a range of dimension 1
  /// @{

  /// @brief Expand a range by decreasing the minimum value
  ///
  /// Expand the range by decreasing the minimum value. If the given value is
  /// larger than the current minimum (in other words, if the proposed new
  /// range would be smaller than the current range), this is a no-op.
  ///
  /// @param v The proposed new minimum for the range
  void expandMin(const Type& v)
    requires(Dims == 1)
  {
    min() = std::min(min(), v);
  }

  /// @brief Expand a range by increasing the maximum value
  ///
  /// Expand the range by increasing the maximum value. If the given value is
  /// smaller than the current maximum (in other words, if the proposed new
  /// range would be smaller than the current range), this is a no-op.
  ///
  /// @param v The proposed new maximum for the range
  void expandMax(const Type& v)
    requires(Dims == 1)
  {
    max() = std::max(max(), v);
  }

  /// @brief Expand a range on both ends
  ///
  /// Expand a range by decreasing the minimum value as well as increasing the
  /// maximum value. If either of the values are already larger or smaller
  /// (respectively) than the proposed values, then that particular boundary
  /// of the interval is not expanded.
  ///
  /// @note After this operation, the range is always equal to or larger than
  /// [min, max].
  ///
  /// @param min The proposed new minimum for the range
  /// @param max The proposed new maximum for the range
  void expand(const Type& min, const Type& max)
    requires(Dims == 1)
  {
    expandMin(min);
    expandMax(max);
  }

  /// @brief Shrink a range by increasing the minimum value
  ///
  /// Shrink the range by increasing the minimum value. If the given value is
  /// smaller than the current minimum (in other words, if the proposed new
  /// range would be larger than the current range), this is a no-op.
  ///
  /// @param v The proposed new minimum for the range
  void shrinkMin(const Type& v)
    requires(Dims == 1)
  {
    min() = std::max(min(), v);
  }

  /// @brief Shrink a range by decreasing the maximum value
  ///
  /// Shrink the range by decreasing the maximum value. If the given value is
  /// larger than the current maximum (in other words, if the proposed new
  /// range would be larger than the current range), this is a no-op.
  ///
  /// @param v The proposed new maximum for the range
  void shrinkMax(const Type& v)
    requires(Dims == 1)
  {
    max() = std::min(max(), v);
  }

  /// @brief Shrink a range on both ends
  ///
  /// Shrink a range by increasing the minimum value as well as decreasing the
  /// maximum value. If either of the values are already smaller or larger
  /// (respectively) than the proposed values, then that particular boundary
  /// of the interval is not shrunk.
  ///
  /// @note After this operation, the range is always equal to or smaller than
  /// [min, max].
  ///
  /// @param min The proposed new minimum for the range
  /// @param max The proposed new maximum for the range
  void shrink(const Type& min, const Type& max)
    requires(Dims == 1)
  {
    shrinkMin(min);
    shrinkMax(max);
  }

  /// @brief Set the minimum value
  ///
  /// Override the minimum value of the range, regardless of what was already
  /// set.
  ///
  /// @note If you want to shrink or expand the range, use the shrink and
  /// expand methods.
  ///
  /// @param v The value to use as the new minimum
  void setMin(const Type& v)
    requires(Dims == 1)
  {
    min() = v;
  }

  /// @brief Set the maximum value
  ///
  /// Override the maximum value of the range, regardless of what was already
  /// set.
  ///
  /// @note If you want to shrink or expand the range, use the shrink and
  /// expand methods.
  ///
  /// @param v The value to use as the new maximum
  void setMax(const Type& v)
    requires(Dims == 1)
  {
    max() = v;
  }

  /// @brief Set the minimum and maximum value
  ///
  /// Override both the minimum and maximum value of the range, regardless of
  /// what they were set to.
  ///
  /// @note If you want to shrink or expand the range, use the shrink and
  /// expand methods.
  ///
  /// @note After this operation, the range should be exactly equal to [min,
  /// max]
  ///
  /// @param min The new minimum value of the range
  /// @param max The new maximum value of the range
  void set(const Type& min, const Type& max)
    requires(Dims == 1)
  {
    setMin(min);
    setMax(max);
  }

  /// @brief Return the minimum value of the range (inclusive)
  Type min() const
    requires(Dims == 1)
  {
    return min(0);
  }

  /// @brief Return the minimum value of the range (inclusive)
  Type& min()
    requires(Dims == 1)
  {
    return min(0);
  }

  /// @brief Return the maximum value of the range (inclusive)
  Type max() const
    requires(Dims == 1)
  {
    return max(0);
  }

  /// @brief Return the maximum value of the range (inclusive)
  Type& max()
    requires(Dims == 1)
  {
    return max(0);
  }

  /// @brief Compute the size of the range
  ///
  /// The size of a range is defined as the difference between the minimum and
  /// the maximum. For degenerate ranges, this is zero.
  ///
  /// @warning Due to the nature of numbers, the result of this function can be
  /// somewhat ambiguous. For natural numbers, you could argue that the range
  /// [n, n] has size 0 or size 1. In this case we say it has size 0. The
  /// uncountable nature of the reals means this doesn't matter for them, but
  /// this can be awkward when working with integers.
  ///
  /// @return The size of the range
  Type size() const
    requires(Dims == 1)
  {
    return std::max(static_cast<Type>(0), max() - min());
  }

  /// @brief Determine if the range contains a given value
  ///
  /// A value is inside a range if and only if it is greater than the minimum
  /// and smaller than the maximum.
  ///
  /// @param v The value to check
  ///
  /// @return true The value is inside the range
  /// @return false The value is not inside the range
  bool contains(const Type& v) const
    requires(Dims == 1)
  {
    return min() <= v && v < max();
  }

  /// @}

 private:
  Vector<Type, Dims> m_minima{};
  Vector<Type, Dims> m_maxima{};
};

template <typename Type,
          template <typename, std::size_t> typename Vector = std::array>
using Range1D = RangeXD<1, Type, Vector>;

}  // namespace Acts
