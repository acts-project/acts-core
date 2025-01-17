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
#include "Acts/Utilities/BoundingBox.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <iosfwd>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace Acts {

class RectangleBounds;

/// @class CuboidVolumeBounds
///
/// Bounds for a cubical Volume, the orientedSurfaces(...) method creates a
/// vector of 6 surfaces:
///
///  BoundarySurfaceFace [index]:
///
///    - negativeFaceXY [0] : Rectangular Acts::PlaneSurface, parallel to \f$ xy
/// \f$ plane at negative \f$ z \f$
///    - positiveFaceXY [1] : Rectangular Acts::PlaneSurface, parallel to \f$ xy
/// \f$ plane at positive \f$ z \f$
///    - negativeFaceXY [2] : Rectangular Acts::PlaneSurface, attached to \f$ yz
/// \f$ plane at negative \f$ x \f$
///    - positiveFaceXY [3] : Rectangular Acts::PlaneSurface, attached to \f$ yz
/// \f$ plane at negative \f$ x \f$
///    - negativeFaceXY [4] : Rectangular Acts::PlaneSurface, parallel to \f$ zx
/// \f$ plane at negative \f$ y \f$
///    - positiveFaceXY [5] : Rectangular Acts::PlaneSurface, parallel to \f$ zx
/// \f$ plane at positive \f$ y \f$
///
class CuboidVolumeBounds : public VolumeBounds {
 public:
  /// @enum BoundValues for streaming and access
  enum BoundValues : unsigned int {
    eHalfLengthX = 0,
    eHalfLengthY = 1,
    eHalfLengthZ = 2,
    eSize
  };

  CuboidVolumeBounds() = delete;

  /// Constructor - the box boundaries
  ///
  /// @param halex is the half length of the cube in x
  /// @param haley is the half length of the cube in y
  /// @param halez is the half length of the cube in z
  CuboidVolumeBounds(double halex, double haley, double halez) noexcept(false);

  /// Constructor - from a fixed size array
  ///
  /// @param values iw the bound values
  explicit CuboidVolumeBounds(const std::array<double, eSize>& values);

  /// Copy Constructor
  ///
  /// @param bobo is the source volume bounds to be copied
  CuboidVolumeBounds(const CuboidVolumeBounds& bobo) = default;

  /// Assignment operator
  ///
  /// @param bobo is the source volume bounds to be assigned
  CuboidVolumeBounds& operator=(const CuboidVolumeBounds& bobo) = default;

  ~CuboidVolumeBounds() override = default;

  VolumeBounds::BoundsType type() const final { return VolumeBounds::eCuboid; }

  /// Return the bound values as dynamically sized vector
  ///
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

  /// This method checks if position in the 3D volume
  /// frame is inside the cylinder
  ///
  /// @param pos is the position in volume frame to be checked
  /// @param tol is the absolute tolerance to be applied
  bool inside(const Vector3& pos, double tol = 0.) const override;

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
  /// @param trf Optional transform
  /// @param envelope Optional envelope to add / subtract from min/max
  /// @param entity Entity to associate this bounding box with
  /// @return Constructed bounding box
  Volume::BoundingBox boundingBox(const Transform3* trf = nullptr,
                                  const Vector3& envelope = {0, 0, 0},
                                  const Volume* entity = nullptr) const final;

  /// Get the canonical binning direction, i.e. the binning directions
  /// for that fully describe the shape's extent
  ///
  /// @return vector of canonical binning values
  std::vector<AxisDirection> canonicalAxes() const override {
    using enum AxisDirection;
    return {AxisX, AxisY, AxisZ};
  };

  /// Binning borders in double
  ///
  /// @param aDir is the axis direction for which the
  /// reference border is requested
  ///
  /// @return float offset to be used for the binning
  double referenceBorder(AxisDirection aDir) const final;

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  double get(BoundValues bValue) const { return m_values[bValue]; }

  /// Set a bound value
  /// @param bValue the bound value identifier
  /// @param value the value to be set
  void set(BoundValues bValue, double value);

  /// Set a range of bound values
  /// @param keyValues the initializer list of key value pairs
  void set(std::initializer_list<std::pair<BoundValues, double>> keyValues);

  /// Output Method for std::ostream
  ///
  /// @param os is ostream operator to be dumped into
  std::ostream& toStream(std::ostream& os) const override;

 private:
  /// The bound values ordered in a fixed size array
  std::array<double, eSize> m_values;

  std::shared_ptr<const RectangleBounds> m_xyBounds{nullptr};
  std::shared_ptr<const RectangleBounds> m_yzBounds{nullptr};
  std::shared_ptr<const RectangleBounds> m_zxBounds{nullptr};

  /// Create the surface bounds
  void buildSurfaceBounds();

  /// Check the input values for consistency,
  /// will throw a logic_exception if consistency is not given
  void checkConsistency() noexcept(false);
};
}  // namespace Acts
