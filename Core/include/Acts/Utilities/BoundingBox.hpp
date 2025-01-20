// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Frustum.hpp"
#include "Acts/Utilities/Ray.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"

#include <array>
#include <memory>
#include <tuple>
#include <vector>

namespace Acts {

/// Implementation of an Axis Aligned Bounding Box. This type is compatible
/// with 2D and 3D boxes
template <typename entity_t, typename value_t, std::size_t DIM>
class AxisAlignedBoundingBox {
 private:
  /// Private self type to capture template parameters
  using self_t = AxisAlignedBoundingBox<entity_t, value_t, DIM>;

  /// Strong type helper, not public
  /// This is only used to provide sensible tag-dispatch below.
  template <typename T, typename P>
  class NamedType {
   public:
    explicit NamedType(const T& value) : m_value(value) {}
    explicit NamedType(T&& value) : m_value(std::move(value)) {}
    T& get() { return m_value; }
    const T& get() const { return m_value; }

   private:
    T m_value;
  };

  /// SizeParameter Tag
  struct SizeParameter {};

 public:
  /// The value type used by this class
  using value_type = value_t;

  /// Re-export vertex type based on value type given
  using VertexType = Eigen::Matrix<value_t, DIM, 1>;

  /// Associated array value to `VertexType`
  using vertex_array_type = Eigen::Array<value_t, DIM, 1>;

  /// Type of stored entity
  using entity_type = entity_t;

  /// The transform type based on the `value_type`
  using transform_type = Eigen::Transform<value_type, DIM, Eigen::Affine>;

  /// Strong type to select the correct constructor
  using Size = NamedType<VertexType, struct SizeParameter>;

  /// Re-export dimension from template parameter
  static const std::size_t dim = DIM;

  /// Copy constructor from other bounding box.
  AxisAlignedBoundingBox(const self_t& other) = default;

  /// Copy assignment operator from other bounding box.
  /// @param other The other AABB
  AxisAlignedBoundingBox& operator=(const self_t& other) = default;

  /// Constructor from an entity pointer, and the min and max vertices.
  /// @param entity The entity to store
  /// @param vmin The minimum vertex.
  /// @param vmax The maximum vertex.
  AxisAlignedBoundingBox(const entity_t* entity, const VertexType& vmin,
                         const VertexType& vmax);

  /// Constructor from a center position, and a width and height.
  /// @param entity The entity to store
  /// @param center The center position
  /// @param size The size (width and height) of the box.
  /// @note The special type @c size is required to disambiguate this constructor
  /// from the other one above. It is a wrapper around a simple @c Vector3.
  AxisAlignedBoundingBox(const entity_t* entity, const VertexType& center,
                         const Size& size);

  /// Constructor from a list of child boxes. This box will wrap around all
  /// boxes
  /// contained in @p boxes, and additional envelope can be given.
  /// @param boxes Vector of child boxes to store in this bounding box.
  /// @param envelope Envelope that will be added/subtracted to the dimension.
  explicit AxisAlignedBoundingBox(
      const std::vector<self_t*>& boxes,
      vertex_array_type envelope = vertex_array_type::Zero());

  /// Helper function to calculate the size of a bounding box enclosing @p boxes.
  /// @param boxes The boxes to wrap (const pointers)
  /// @param envelope Optional envelop to add/subtract to dimension.
  /// @return Pair of vertices: min and max.
  static std::pair<VertexType, VertexType> wrap(
      const std::vector<const self_t*>& boxes,
      vertex_array_type envelope = vertex_array_type::Zero());

  /// Helper function to calculate the size of a bounding box enclosing @p boxes.
  /// Overload which accepts non-const boxes in @p boxes.
  /// @param boxes The boxes to wrap (non-const pointers)
  /// @param envelope Optional envelop to add/subtract to dimension.
  /// @return Pair of vertices: min and max.
  static std::pair<VertexType, VertexType> wrap(
      const std::vector<self_t*>& boxes,
      vertex_array_type envelope = vertex_array_type::Zero());

  /// Helper function to calculate the size of a bounding box enclosing @p boxes.
  /// Overload which accepts a vector in @p boxes which owns the instances
  /// @param boxes The boxes to wrap (by-value vector)
  /// @param envelope Optional envelop to add/subtract to dimension.
  /// @return Pair of vertices: min and max.
  static std::pair<VertexType, VertexType> wrap(
      const std::vector<self_t>& boxes,
      vertex_array_type envelope = vertex_array_type::Zero());

  /// Calculate whether a point is inside this box.
  /// @param point The point to test.
  /// @return Whether the point is inside or not.
  bool intersect(const VertexType& point) const;

  /// @brief Implements the slab method for Ray/AABB intersections.
  ///
  /// See https://tavianator.com/fast-branchless-raybounding-box-intersections/,
  /// https://tavianator.com/fast-branchless-raybounding-box-intersections-part-2-nans/,
  /// https://medium.com/@bromanz/another-view-on-the-classic-ray-aabb-intersection-algorithm-for-bvh-traversal-41125138b525
  /// The original algorithms is described in "Graphics Gems (1990)" [1]
  /// (https://doi.org/10.1016/B978-0-08-050753-8.50084-X)
  ///
  /// @note This implementation may treat parallel rays on any of the slabs
  ///       as **outside** due to how @c NaNs are handled by Eigen.
  ///       See https://eigen.tuxfamily.org/bz/show_bug.cgi?id=564
  /// @param ray The ray to intersect with
  /// @return Whether the ray intersects this AABB
  bool intersect(const Ray<value_type, DIM>& ray) const;

  /// Check if a frustum intersects with this bounding box.
  ///
  /// This method implements an algorithm similar to the one described in
  /// "Optimized View Frustum Culling Algorithms for Bounding Boxes (2012)" [2]
  /// (https://doi.org/10.1080/10867651.2000.10487517), but drops some of the
  /// more sophisticated optimization.
  ///
  /// @param fr The frustum
  /// @return Whether the frustum intersects this AABB
  template <std::size_t sides>
  bool intersect(const Frustum<value_type, DIM, sides>& fr) const;

  /// Set the skip node (bounding box)
  /// @param skip The target skip node pointer
  void setSkip(self_t* skip);

  /// Get the skip node for this box
  /// @return The skip node pointer
  const self_t* getSkip() const;

  /// Get the left child (i.e. the first of the children that are inside this
  /// bounding box).
  /// @return The lest most child.
  const self_t* getLeftChild() const;

  /// Check whether this node as an associated entity. If it does not have one,
  /// this is a purely abstract container box.
  /// @return Whether the box has an entity attached.
  bool hasEntity() const;

  /// Return the entity associated with this box. This might be nullptr if there
  /// is no entity attached.
  /// @return The entity pointer, might be nullptr
  const entity_t* entity() const;

  /// Set the entity associated with with this box.
  /// @param entity The entity
  void setEntity(const entity_t* entity);

  /// Get the center position of this bounding box.
  /// @return The center position
  const VertexType& center() const;

  /// Get the minimum vertex
  /// @return The minimum vertex
  const VertexType& min() const;

  /// Get the maximum vertex
  /// @return The maximum vertex
  const VertexType& max() const;

  /// Write information about this bounding box to a stream.
  /// @param os The output stream.
  /// @return The stream given as an argument.
  std::ostream& toStream(std::ostream& os) const;

  /// Transforms this bounding box using the given transform. This method
  /// modifies the box it is called on.
  /// @param trf The transform
  void transform(const transform_type& trf);

  /// Transforms this bounding box using the given transform. This method
  /// returns a copy of this box, with the transformation applied, and leaves
  /// this instance unchanged.
  /// @param trf The transform
  /// @return The transformed bounding box
  self_t transformed(const transform_type& trf) const;

  /// Draw this bounding box using the given visualization helper. This method
  /// is only available for the 3D case.
  /// @param helper The visualization helper to write to
  /// @param color The color to use for drawing
  /// @param trf An optional transform to apply first.
  void draw(IVisualization3D& helper, Color color = {120, 120, 120},
            const transform_type& trf = transform_type::Identity()) const
    requires(DIM == 3);

  /// Draw this bounding box as SVG. This method is only available for the 2D
  /// case.
  /// @param os The output stream to write to
  /// @param w The width of the output SVG.
  /// @param h The height of the output SVG.
  /// @param unit A scale factor to apply before drawing
  /// @param label A label to put next to the box.
  /// @param fillcolor Color to fill the box with.
  /// @return The outstream given in @p os.
  std::ostream& svg(std::ostream& os, value_type w, value_type h,
                    value_type unit = 10, const std::string& label = "",
                    const std::string& fillcolor = "grey") const
    requires(DIM == 2);

 private:
  std::pair<VertexType, VertexType> transformVertices(
      const transform_type& trf) const
    requires(DIM == 2);

  std::pair<VertexType, VertexType> transformVertices(
      const transform_type& trf) const
    requires(DIM == 3);

  const entity_t* m_entity;
  VertexType m_vmin;
  VertexType m_vmax;
  VertexType m_center;
  vertex_array_type m_width;
  vertex_array_type m_iwidth;

  self_t* m_left_child{nullptr};
  self_t* m_right_child{nullptr};
  self_t* m_skip{nullptr};
};

/// Build an octree from a list of bounding boxes.
/// @note @p store and @p prims do not need to contain the same objects. @p store
/// is only used to pass ownership back to the caller while preserving memory
/// location.
/// @tparam box_t Works will all box types.
/// @param store Owns the created boxes by means of `std::unique_ptr`.
/// @param prims Boxes to store. This is a read only vector.
/// @param max_depth No subdivisions beyond this level.
/// @param envelope1 Envelope to add/subtract to dimensions in all directions.
/// @return Pointer to the top most bounding box, containing the entire octree
template <typename box_t>
box_t* make_octree(std::vector<std::unique_ptr<box_t>>& store,
                   const std::vector<box_t*>& prims, std::size_t max_depth = 1,
                   typename box_t::value_type envelope1 = 0);

/// Overload of the << operator for bounding boxes.
/// @tparam T entity type
/// @tparam U value type
/// @tparam V dimension
/// @param os The output stream
/// @param box The bounding box
/// @return The given output stream.
template <typename T, typename U, std::size_t V>
std::ostream& operator<<(std::ostream& os,
                         const AxisAlignedBoundingBox<T, U, V>& box);

}  // namespace Acts

#include "Acts/Utilities/BoundingBox.ipp"
