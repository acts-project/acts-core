// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

#include <array>
#include <iosfwd>
#include <memory>
#include <vector>

namespace Acts {

class CylinderBounds;
class DiscBounds;

/// Class which implements a cutout cylinder. This shape is basically a
/// cylinder, with another, smaller cylinder subtracted from the center.
/// --------------------- rmax
/// |                   |
/// |    |---------|    | rmed
/// |    |         |    |
/// ------         ------ rmin
///       -- hlZc --
/// --------- hlZ -------
///
///
/// @todo add sectoral cutouts
class CutoutCylinderVolumeBounds : public VolumeBounds {
 public:
  /// @enum BoundValues for streaming and access
  enum BoundValues : int {
    eMinR = 0,
    eMedR = 1,
    eMaxR = 2,
    eHalfLengthZ = 3,
    eHalfLengthZcutout = 4,
    eSize
  };

  /// Constructor from defining parameters
  ///
  /// @param rmin Minimum radius at the "choke points"
  /// @param rmed The medium radius (outer radius of the cutout)
  /// @param rmax The outer radius of the overall shape
  /// @param hlZ The longer halflength of the shape
  /// @param hlZc The cutout halflength of the shape
  CutoutCylinderVolumeBounds(double rmin, double rmed, double rmax, double hlZ,
                             double hlZc) noexcept(false)
      : m_values({rmin, rmed, rmax, hlZ, hlZc}) {
    checkConsistency();
    buildSurfaceBounds();
  }

  /// Constructor - from a fixed size array
  ///
  /// @param values The bound values
  explicit CutoutCylinderVolumeBounds(
      const std::array<double, eSize>& values) noexcept(false)
      : m_values(values) {
    checkConsistency();
    buildSurfaceBounds();
  }

  VolumeBounds::BoundsType type() const final {
    return VolumeBounds::eCutoutCylinder;
  }

  /// Return the bound values as dynamically sized vector
  ///
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

  /// Inside method to test whether a point is inside the shape
  ///
  /// @param gpos The point to test
  /// @param tol The tolerance to test with
  /// @return Whether the point is inside or not.
  bool inside(const Vector3& gpos, double tol = 0) const override;

  /// Oriented surfaces, i.e. the decomposed boundary surfaces and the
  /// according navigation direction into the volume given the normal
  /// vector on the surface
  ///
  /// @param transform is the 3D transform to be applied to the boundary
  /// surfaces to position them in 3D space
  ///
  /// It will throw an exception if the orientation prescription is not adequate
  ///
  /// @return a vector of surfaces bounding this volume
  std::vector<OrientedSurface> orientedSurfaces(
      const Transform3& transform = Transform3::Identity()) const override;

  /// Construct bounding box for this shape
  ///
  /// @param trf Optional transform
  /// @param envelope Optional envelope to add / subtract from min/max
  /// @param entity Entity to associate this bounding box with
  /// @return Constructed bounding box
  Volume::BoundingBox boundingBox(const Transform3* trf = nullptr,
                                  const Vector3& envelope = {0, 0, 0},
                                  const Volume* entity = nullptr) const final;

  /// Get the canonical binning direction, i.e. the axis values
  /// that fully describe the shape's extent
  ///
  /// @return vector of canonical binning values
  std::vector<AxisDirection> canonicalAxes() const override {
    using enum AxisDirection;
    return {AxisR, AxisPhi, AxisZ};
  };

  /// Write information about this instance to an outstream
  ///
  /// @param sl The outstream
  /// @return The outstream
  std::ostream& toStream(std::ostream& sl) const override;

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  double get(BoundValues bValue) const { return m_values[bValue]; }

 private:
  std::array<double, eSize> m_values;

  // The surface bound objects
  std::shared_ptr<const CylinderBounds> m_innerCylinderBounds{nullptr};
  std::shared_ptr<const CylinderBounds> m_cutoutCylinderBounds{nullptr};
  std::shared_ptr<const CylinderBounds> m_outerCylinderBounds{nullptr};
  std::shared_ptr<const DiscBounds> m_outerDiscBounds{nullptr};
  std::shared_ptr<const DiscBounds> m_innerDiscBounds{nullptr};

  /// Create the surface bound objects
  void buildSurfaceBounds();

  /// Check the input values for consistency,
  /// will throw a logic_exception if consistency is not given
  void checkConsistency() noexcept(false);
};

}  // namespace Acts
