// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <array>
#include <cstddef>
#include <iosfwd>
#include <vector>

#include <boost/container/small_vector.hpp>

namespace Acts {

/// @brief base class for convex polygon bounds
///
/// This class serves as a base class for the actual bounds class.
/// The only deriving type is the templated `ConvexPolygonBounds`.
///
class ConvexPolygonBoundsBase : public PlanarBounds {
 public:
  /// Output Method for std::ostream
  /// @param sl is the ostream to be written into
  std::ostream& toStream(std::ostream& sl) const final;

  /// Return the bounds type of this bounds object.
  /// @return The bounds type
  BoundsType type() const final { return eConvexPolygon; }

  /// Return the bound values as dynamically sized vector
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

 protected:
  /// Return a rectangle bounds instance that encloses a set of vertices.
  /// @param vertices A collection of vertices to enclose.
  /// @return Enclosing rectangle.
  template <typename coll_t>
  static RectangleBounds makeBoundingBox(const coll_t& vertices);

  /// Calculates whether a set of vertices forms a convex polygon. This is
  /// generic over the number of vertices, so it's factored out of the concrete
  /// classes and into this base class.
  /// @param vertices A collection of vertices.
  /// throws a logic error if this is not the case
  template <typename coll_t>
    requires std::same_as<typename coll_t::value_type, Acts::Vector2>
  static void convex_impl(const coll_t& vertices) noexcept(false);
};

/// This is the actual implementation of the bounds.
/// It is templated on the number of vertices, but there is a specialization for
/// *dynamic* number of vertices, where the underlying storage is then a vector.
///
/// @tparam N Number of vertices
template <int N>
class ConvexPolygonBounds : public ConvexPolygonBoundsBase {
 public:
  /// Expose number of vertices given as template parameter.
  ///
  static constexpr std::size_t num_vertices = N;
  /// Type that's used to store the vertices, in this case a fixed size array.
  ///
  using vertex_array = std::array<Vector2, num_vertices>;
  /// Expose number of parameters as a template parameter
  ///
  static constexpr std::size_t eSize = 2 * N;
  /// Type that's used to store the vertices, in this case a fixed size array.
  ///
  using value_array = std::array<double, eSize>;

  static_assert(N >= 3, "ConvexPolygonBounds needs at least 3 sides.");

  /// Constructor from a vector of vertices, to facilitate construction.
  /// This will throw if the vector size does not match `num_vertices`.
  /// This will throw if the vertices do not form a convex polygon.
  /// @param vertices The list of vertices.
  ConvexPolygonBounds(const std::vector<Vector2>& vertices) noexcept(false);

  /// Constructor from a fixed size array of vertices.
  /// This will throw if the vertices do not form a convex polygon.
  /// @param vertices The vertices
  ConvexPolygonBounds(const vertex_array& vertices) noexcept(false);

  /// Constructor from a fixed size array of parameters
  /// This will throw if the vertices do not form a convex polygon.
  /// @param values The values to build up the vertices
  ConvexPolygonBounds(const value_array& values) noexcept(false);

  /// @copydoc SurfaceBounds::inside
  bool inside(const Vector2& lposition) const final;

  /// @copydoc SurfaceBounds::closestPoint
  Vector2 closestPoint(const Vector2& lposition,
                       const std::optional<SquareMatrix2>& metric) const final;

  using SurfaceBounds::inside;

  /// Return the vertices
  ///
  /// @param ignoredSegments the number of segments used to approximate
  /// and eventually curved line
  ///
  /// @note the number of segments is ignored in this representation
  ///
  /// @return vector for vertices in 2D
  std::vector<Vector2> vertices(unsigned int ignoredSegments = 0u) const final;

  /// Return a rectangle bounds object that encloses this polygon.
  /// @return The rectangular bounds
  const RectangleBounds& boundingBox() const final;

 private:
  vertex_array m_vertices;
  RectangleBounds m_boundingBox;

  /// Return whether this bounds class is in fact convex
  /// throws a log error if not
  void checkConsistency() const noexcept(false);
};

/// Tag to trigger specialization of a dynamic polygon
constexpr int PolygonDynamic = -1;

/// This is the specialization handling a polygon with a dynamic number of
/// points. It can accept any number of points.
///
template <>
class ConvexPolygonBounds<PolygonDynamic> : public ConvexPolygonBoundsBase {
 public:
  constexpr static int eSize = -1;

  /// Constructor from a vector of vertices, to facilitate construction.
  /// This will throw if the vertices do not form a convex polygon.
  /// @param vertices The list of vertices.
  ConvexPolygonBounds(const std::vector<Vector2>& vertices);

  /// @copydoc SurfaceBounds::inside
  bool inside(const Vector2& lposition) const final;

  /// @copydoc SurfaceBounds::closestPoint
  Vector2 closestPoint(const Vector2& lposition,
                       const std::optional<SquareMatrix2>& metric) const final;

  using SurfaceBounds::inside;

  /// Return the vertices
  ///
  /// @param lseg the number of segments used to approximate
  /// and eventually curved line
  ///
  /// @note the number of segments is ignored in this representation
  ///
  /// @return vector for vertices in 2D
  std::vector<Vector2> vertices(unsigned int lseg = 1) const final;

  ///
  /// Return a rectangle bounds object that encloses this polygon.
  /// @return The rectangular bounds
  ///
  const RectangleBounds& boundingBox() const final;

 private:
  boost::container::small_vector<Vector2, 10> m_vertices;
  RectangleBounds m_boundingBox;

  /// Return whether this bounds class is in fact convex
  /// thorws a logic error if not
  void checkConsistency() const noexcept(false);
};

}  // namespace Acts

#include "Acts/Surfaces/ConvexPolygonBounds.ipp"
