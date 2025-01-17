// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/BoundarySurfaceFace.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

#include <array>
#include <initializer_list>
#include <iosfwd>
#include <memory>
#include <numbers>
#include <ostream>
#include <vector>

namespace Acts {

class CylinderBounds;
class RadialBounds;
class PlanarBounds;

/// @class CylinderVolumeBounds
///
/// Bounds for a cylindrical Volume, the orientedSurfaces(..) method creates a
/// vector of up to 6 surfaces:
///
/// case A) 3 Surfaces (full cylindrical tube):
///  BoundarySurfaceFace [index]:
///  - negativeFaceXY [0] : Acts::DiscSurface with \f$ r_{inner}=0 \f$,
///                         parallel to \f$ xy \f$ plane at negative \f$ z\f$
///  - positiveFaceXY [1] : Acts::DiscSurface with \f$ r_{inner}=0 \f$,
///                         parallel to \f$ xy \f$ plane at positive \f$ z\f$
///  - cylinderCover  [2] : Acts::CylinderSurface confining the Acts::Volume
///
/// case B) 4 Surfaces (tube with inner and outer radius):
///  BoundarySurfaceFace [index]:
///  - negativeFaceXY [0] : Acts::DiscSurface with \f$ r_{inner}>0 \f$,
///                         parallel to \f$ xy \f$ plane at negative \f$ z\f$
///  - positiveFaceXY [1] : Acts::DiscSurface with \f$ r_{inner}>0 \f$,
///                         parallel to \f$ xy \f$ plane at positive \f$ z\f$
///  - tubeOuterCover [2] : Acts::CylinderSurface with \f$ r = r_{outer} \f$
///  - tubeInnerCover [3] : Acts::CylinderSurface with \f$ r = r_{inner} \f$
///
/// case C) 6 Surfaces (sectoral tube with inner and outer radius):
///  BoundarySurfaceFace [index]:
///   - negativeFaceXY  [0] : Acts::DiscSurface with \f$ r_{inner}>0\f$
///                           and \f$ \phi < \pi \f$,
///                           parallel to \f$ xy \f$ plane at negative \f$z\f$
///   - positiveFaceXY  [1] : Acts::DiscSurface with \f$ r_{inner}>0 \f$
///                           and \f$ \phi < \pi \f$,
///                           parallel to \f$ xy \f$ plane at positive \f$z\f$
///   - tubeSectorOuterCover  [2] : Acts::CylinderSurface with
///                                 \f$ r = r_{outer}\f$
///   - tubeSectorInnerCover  [3] : Acts::CylinderSurface with
///                                 \f$ r = r_{inner} \f$
///   - tubeSectorNegativePhi [4] : Rectangular Acts::PlaneSurface attached to
///                 [0] and [1] at negative \f$ \phi \f$
///                      - tubeSectorNegativePhi [5] :
//                          Rectangular Acts::PlaneSurface attached to
///                 [0] and [1] at positive \f$ \phi \f$
///
class CylinderVolumeBounds : public VolumeBounds {
 public:
  /// @enum BoundValues for streaming and access
  enum BoundValues : unsigned int {
    eMinR = 0,
    eMaxR = 1,
    eHalfLengthZ = 2,
    eHalfPhiSector = 3,
    eAveragePhi = 4,
    eBevelMinZ = 5,
    eBevelMaxZ = 6,
    eSize
  };

  /// Enum describing the possible faces of a cylinder volume
  /// @note These values are synchronized with the BoundarySurfaceFace enum.
  ///       Once Gen1 is removed, this can be changed.
  enum class Face : unsigned int {
    PositiveDisc = BoundarySurfaceFace::positiveFaceXY,
    NegativeDisc = BoundarySurfaceFace::negativeFaceXY,
    OuterCylinder = BoundarySurfaceFace::tubeOuterCover,
    InnerCylinder = BoundarySurfaceFace::tubeInnerCover,
    NegativePhiPlane = BoundarySurfaceFace::tubeSectorNegativePhi,
    PositivePhiPlane = BoundarySurfaceFace::tubeSectorPositivePhi
  };

  CylinderVolumeBounds() = delete;

  /// Constructor
  ///
  /// @param rmin The inner radius of the cylinder
  /// @param rmax The outer radius of the cylinder
  /// @param halfz The half length in z
  /// @param halfphi The half lopening angle
  /// @param avgphi The average phi value
  /// @param bevelMinZ The bevel angle, in radians, for the negative side
  /// @param bevelMaxZ The bevel angle, in radians, for the positive side
  CylinderVolumeBounds(double rmin, double rmax, double halfz,
                       double halfphi = std::numbers::pi, double avgphi = 0.,
                       double bevelMinZ = 0., double bevelMaxZ = 0.);

  /// Constructor - from a fixed size array
  ///
  /// @param values The bound values
  explicit CylinderVolumeBounds(const std::array<double, eSize>& values);

  /// Constructor - extruded from cylinder bounds and thickness
  ///
  /// @param cBounds the cylinder bounds
  /// @param thickness of the extrusion
  CylinderVolumeBounds(const CylinderBounds& cBounds, double thickness);

  /// Constructor - extruded from radial bounds and thickness
  ///
  /// @param rBounds the Radial bounds
  /// @param thickness
  CylinderVolumeBounds(const RadialBounds& rBounds, double thickness);

  /// Copy Constructor
  ///
  /// @param cylbo is the source cylinder volume bounds for the copy
  CylinderVolumeBounds(const CylinderVolumeBounds& cylbo);

  ~CylinderVolumeBounds() override = default;
  CylinderVolumeBounds& operator=(const CylinderVolumeBounds& cylbo) = default;

  VolumeBounds::BoundsType type() const final {
    return VolumeBounds::eCylinder;
  }

  /// Return the bound values as dynamically sized vector
  ///
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

  /// This method checks if position in the 3D volume
  /// frame is inside the cylinder
  ///
  /// @param pos is a global position to be checked
  /// @param tol is the tolerance for the check
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

  /// Get the canonical binning directions, i.e. the axis directions
  /// that fully describe the shape's extent
  ///
  /// @return vector of canonical binning values
  std::vector<AxisDirection> canonicalAxes() const override {
    using enum AxisDirection;
    return {AxisR, AxisPhi, AxisZ};
  };

  /// Binning offset - overloaded for some R-binning types
  ///
  /// @param aDir is the axis direction used for the binning
  Vector3 referenceOffset(AxisDirection aDir) const override;

  /// Binning borders in double
  ///
  /// @param aDir is the axis direction used for the binning
  double referenceBorder(AxisDirection aDir) const override;

  /// Output Method for std::ostream
  /// @param os is the output stream
  std::ostream& toStream(std::ostream& os) const override;

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

 private:
  /// The internal version of the bounds can be float/double
  std::array<double, eSize> m_values{};

  /// Bounds of the inner CylinderBounds
  std::shared_ptr<const CylinderBounds> m_innerCylinderBounds{nullptr};
  /// Bounds of the inner CylinderBounds
  std::shared_ptr<const CylinderBounds> m_outerCylinderBounds{nullptr};
  /// Bounds of the bottom/top Radial
  std::shared_ptr<const RadialBounds> m_discBounds{nullptr};
  /// Bounds of the sector planes
  std::shared_ptr<const PlanarBounds> m_sectorPlaneBounds{nullptr};

  /// Check the input values for consistency,
  /// will throw a logic_exception if consistency is not given
  void checkConsistency();

  /// Helper method to create the surface bounds
  void buildSurfaceBounds();
};

}  // namespace Acts
