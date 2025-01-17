// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <iterator>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Acts {

/// Store values mapped into the geometry hierarchy.
///
/// @tparam value_t stored value type
///
/// The core functionality is to find an equivalent element, i.e. an
/// identifier-value pair, for a given geometry identifier via
///
///     auto it = container.find(GeometryIdentifier(...));
///     if (it != container.end()) {
///         ...
///     }
///
/// Trailing zero levels of stored geometry identifiers are used as
/// broadcast values to refer to higher-level objects within the geometry, e.g.
/// a geometry identifier with vanishing approach and sensitive index identifies
/// a layer. An entry will all geometry identifier levels set to zero acts as
/// the global default value.
///
/// The container also supports range-based iteration over all stored elements
///
///     for (const auto& element : container) {
///         ...
///     }
///
/// and index-based access to stored elements and associated geometry
/// identifiers
///
///     GeometryIdentifier id3 = container.idAt(3);
///     const auto& element4 = container.valueAt(4);
///
/// @note No guarantees are given for the element order when using range-based
///   or index-based access. Any apparent ordering must be considered an
///   implementation detail and might change.
///
/// Adding elements is potentially expensive as the internal lookup structure
/// must be updated. In addition, modifying an element in-place could change its
/// identifier which would also break the lookup. Thus, the container can not be
/// modified after construction to prevent misuse.
template <typename value_t>
class GeometryHierarchyMap {
 public:
  /// Combined geometry identifier and value element. Only used for input.
  using InputElement = typename std::pair<GeometryIdentifier, value_t>;
  using Iterator = typename std::vector<value_t>::const_iterator;
  using Value = value_t;

  /// Construct the container from the given elements.
  ///
  /// @param elements input elements (must be unique with respect to identifier)
  GeometryHierarchyMap(std::vector<InputElement> elements);

  /// Construct the container from an initializer list.
  ///
  /// @param elements input initializer list
  GeometryHierarchyMap(std::initializer_list<InputElement> elements);

  // defaulted constructors and assignment operators
  GeometryHierarchyMap() = default;
  GeometryHierarchyMap(const GeometryHierarchyMap&) = default;
  GeometryHierarchyMap(GeometryHierarchyMap&&) noexcept = default;
  ~GeometryHierarchyMap() = default;
  GeometryHierarchyMap& operator=(const GeometryHierarchyMap&) = default;
  GeometryHierarchyMap& operator=(GeometryHierarchyMap&&) noexcept = default;

  /// Return an iterator pointing to the beginning of the stored values.
  Iterator begin() const { return m_values.begin(); }

  /// Return an iterator pointing to the end of the stored values.
  Iterator end() const { return m_values.end(); }

  /// Check if any elements are stored.
  bool empty() const { return m_values.empty(); }

  /// Return the number of stored elements.
  std::size_t size() const { return m_values.size(); }

  /// Access the geometry identifier for the i-th element with bounds check.
  ///
  /// @thhrows std::out_of_range for invalid indices
  GeometryIdentifier idAt(std::size_t index) const {
    return GeometryIdentifier(m_ids.at(index));
  }

  /// Access the value at the given index.
  ///
  /// @throws std::out_of_range for invalid indices
  const Value& valueAt(std::size_t index) const { return m_values.at(index); }

  /// Find the most specific value for a given geometry identifier.
  ///
  /// This can be either from the element matching exactly to the given geometry
  /// id, if it exists, or from the element for the next available higher level
  /// within the geometry hierarchy.
  ///
  /// @param id geometry identifier for which information is requested
  /// @retval iterator to an existing value
  /// @retval `.end()` iterator if no matching element exists
  Iterator find(const GeometryIdentifier& id) const;

  /// Check if the most specific value exists for a given geometry identifier.
  ///
  /// This function checks if there is an element matching exactly the given
  /// geometry id, or from the element for the next available higher level
  /// within the geometry hierarchy.
  ///
  /// @param id geometry identifier for which existence is being checked
  /// @retval `true` if a matching element exists
  /// @retval `false` if no matching element exists
  bool contains(const GeometryIdentifier& id) const;

 private:
  // NOTE this class assumes that it knows the ordering of the levels within
  //      the geometry id. if the geometry id changes, this code has to be
  //      adapted too. the asserts ensure that such a change is caught.
  static_assert(GeometryIdentifier().setVolume(1).value() <
                    GeometryIdentifier().setVolume(1).setBoundary(1).value(),
                "Incompatible GeometryIdentifier hierarchy");
  static_assert(GeometryIdentifier().setBoundary(1).value() <
                    GeometryIdentifier().setBoundary(1).setLayer(1).value(),
                "Incompatible GeometryIdentifier hierarchy");
  static_assert(GeometryIdentifier().setLayer(1).value() <
                    GeometryIdentifier().setLayer(1).setApproach(1).value(),
                "Incompatible GeometryIdentifier hierarchy");
  static_assert(GeometryIdentifier().setApproach(1).value() <
                    GeometryIdentifier().setApproach(1).setSensitive(1).value(),
                "Incompatible GeometryIdentifier hierarchy");

  using Identifier = GeometryIdentifier::Value;

  // encoded ids for all elements for faster lookup.
  std::vector<Identifier> m_ids;
  // validity bit masks for the ids: which parts to use for comparison
  std::vector<Identifier> m_masks;
  std::vector<Value> m_values;

  /// Construct a mask where all leading non-zero levels are set.
  static constexpr Identifier makeLeadingLevelsMask(GeometryIdentifier id) {
    // construct id from encoded value with all bits set
    auto allSet = GeometryIdentifier(~GeometryIdentifier::Value{0u});
    // manually iterate over identifier levels starting from the lowest
    if (id.sensitive() != 0u) {
      // all levels are valid; keep all bits set.
      return allSet.setExtra(0u).value();
    }
    if (id.approach() != 0u) {
      return allSet.setExtra(0u).setSensitive(0u).value();
    }
    if (id.layer() != 0u) {
      return allSet.setExtra(0u).setSensitive(0u).setApproach(0u).value();
    }
    if (id.boundary() != 0u) {
      return allSet.setExtra(0u)
          .setSensitive(0u)
          .setApproach(0u)
          .setLayer(0u)
          .value();
    }
    if (id.volume() != 0u) {
      return allSet.setExtra(0u)
          .setSensitive(0u)
          .setApproach(0u)
          .setLayer(0u)
          .setBoundary(0u)
          .value();
    }
    // no valid levels; all bits are zero.
    return Identifier{0u};
  }

  /// Construct a mask where only the highest level is set.
  static constexpr Identifier makeHighestLevelMask() {
    return makeLeadingLevelsMask(GeometryIdentifier(0u).setVolume(1u));
  }

  /// Compare the two identifiers only within the masked bits.
  static constexpr bool equalWithinMask(Identifier lhs, Identifier rhs,
                                        Identifier mask) {
    return (lhs & mask) == (rhs & mask);
  }

  /// Ensure identifier ordering and uniqueness.
  static void sortAndCheckDuplicates(std::vector<InputElement>& elements);

  /// Fill the container from the input elements.
  ///
  /// This assumes that the elements are ordered and unique with respect to
  /// their identifiers.
  void fill(const std::vector<InputElement>& elements);
};

// implementations

template <typename value_t>
inline GeometryHierarchyMap<value_t>::GeometryHierarchyMap(
    std::vector<InputElement> elements) {
  sortAndCheckDuplicates(elements);
  fill(elements);
}

template <typename value_t>
inline GeometryHierarchyMap<value_t>::GeometryHierarchyMap(
    std::initializer_list<InputElement> elements)
    : GeometryHierarchyMap(
          std::vector<InputElement>(elements.begin(), elements.end())) {}

template <typename value_t>
inline void GeometryHierarchyMap<value_t>::sortAndCheckDuplicates(
    std::vector<InputElement>& elements) {
  // ensure elements are sorted by identifier
  std::ranges::sort(elements, [=](const auto& lhs, const auto& rhs) {
    return lhs.first < rhs.first;
  });

  // Check that all elements have unique identifier
  auto dup = std::ranges::adjacent_find(
      elements,
      [](const auto& lhs, const auto& rhs) { return lhs.first == rhs.first; });

  if (dup != elements.end()) {
    throw std::invalid_argument("Input elements contain duplicates");
  }
}

template <typename value_t>
inline void GeometryHierarchyMap<value_t>::fill(
    const std::vector<InputElement>& elements) {
  m_ids.clear();
  m_masks.clear();
  m_values.clear();

  m_ids.reserve(elements.size());
  m_masks.reserve(elements.size());
  m_values.reserve(elements.size());

  for (const auto& element : elements) {
    m_ids.push_back(element.first.value());
    m_masks.push_back(
        makeLeadingLevelsMask(GeometryIdentifier(element.first.value())));
    m_values.push_back(std::move(element.second));
  }
}

template <typename value_t>
inline auto GeometryHierarchyMap<value_t>::find(
    const GeometryIdentifier& id) const -> Iterator {
  assert((m_ids.size() == m_values.size()) &&
         "Inconsistent container state: #ids != # values");
  assert((m_masks.size() == m_values.size()) &&
         "Inconsistent container state: #masks != #values");

  // we can not search for the element directly since the relevant one
  // might be stored at a higher level. ids for higher levels would always
  // be sorted before the requested id. searching for the first element
  // after the requested ensures that we include the full hierarchy.
  const auto it = std::upper_bound(m_ids.begin(), m_ids.end(), id.value());
  auto i = std::distance(m_ids.begin(), it);

  // now go up the hierarchy to find the first matching element.
  // example: the container stores four identifiers
  //
  //     2|x|x (volume-only)
  //     2|2|1 (volume, layer, and sensitive)
  //     2|3|x (volume and layer)
  //     2|3|4 (volume, layer, and sensitive)
  //
  // where | marks level boundaries. searching for either 2|3|4, 2|3|7, or
  // 2|4|x would first point to 2|3|4 and thus needs to go up the hierarchy.
  while (0 < i) {
    // index always starts after item of interest due to upper bound search
    --i;

    // if the input id does not even match at the highest hierarchy level
    // with the current comparison id, then have reached the end of this
    // hierarchy. having a special check for the highest level avoids an
    // unbounded search window all the way to the beginning of the container for
    // the global default entry.
    if (!equalWithinMask(id.value(), m_ids[i], makeHighestLevelMask())) {
      // check if a global default entry exists
      if (m_ids.front() == Identifier{0u}) {
        return begin();
      } else {
        return end();
      }
    }

    // since the search is going backwards in the sorted container, it
    // progresses from more specific to less specific elements. the first
    // match is automatically the appropriate one.
    if (equalWithinMask(id.value(), m_ids[i], m_masks[i])) {
      return std::next(begin(), i);
    }
  }

  // all options are exhausted and no matching element was found.
  return end();
}

template <typename value_t>
inline auto GeometryHierarchyMap<value_t>::contains(
    const GeometryIdentifier& id) const -> bool {
  return this->find(id) != this->end();
}

}  // namespace Acts
