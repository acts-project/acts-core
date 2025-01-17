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

#include <array>
#include <iosfwd>
#include <memory>
#include <ostream>
#include <vector>

namespace Acts {

class RectangleBounds;
class TrapezoidBounds;

/// @class TrapezoidVolumeBounds
///
/// Bounds for a trapezoidal shaped Volume, the orientedSurface(...) method
/// creates a vector of 6 surfaces:
///
///  BoundarySurfaceFace [index]:
///
///  - negativeFaceXY     [0] : Trazpezoidal Acts::PlaneSurface,
///                             parallel to \f$ xy \f$ plane at negative \f$z\f$
///  - positiveFaceXY     [1] : Trazpezoidal Acts::PlaneSurface,
///                             parallel to \f$ xy \f$ plane at positive \f$z\f$
///  - trapezoidFaceAlpha [2] : Rectangular  Acts::PlaneSurface,
///                             attached to [0] and [1] at negative \f$x\f$
///  (associated to alpha)
///  - trapezoidFaceBeta  [3] : Rectangular  Acts::PlaneSurface,
///                             attached to [0] and [1] at positive \f$ x\f$
///  (associated to beta)
///  - negativeFaceZX     [4] : Rectangular  Acts::PlaneSurface,
///                             parallel to \f$ zx \f$ plane at negative \f$y\f$
///  - positiveFaceZX     [5] : Rectangular  Acts::PlaneSurface,
///                             parallel to \f$ zx \f$ plane at positive \f$y\f$
///
class TrapezoidVolumeBounds : public VolumeBounds {
 public:
  /// @enum BoundValues for access / streaming
  enum BoundValues : unsigned int {
    eHalfLengthXnegY = 0,  //!< halflength in x at negative y
    eHalfLengthXposY = 1,  //!< halflength in x at positive y
    eHalfLengthY = 2,      //!< halflength in y
    eHalfLengthZ = 3,      //!< halflength in z
    eAlpha = 4,            //!< opening angle alpha (in point A)
    eBeta = 5,             //!< opening angle beta  (in point B)
    eSize                  //!< length of the bounds vector
  };

  TrapezoidVolumeBounds() = delete;

  /// Constructor - the trapezoid boundaries (symmetric trapezoid)
  ///
  /// @param minhalex is the half length in x at minimal y
  /// @param maxhalex is the half length in x at maximal y
  /// @param haley is the half length in y
  /// @param halez is the half length in z
  TrapezoidVolumeBounds(double minhalex, double maxhalex, double haley,
                        double halez) noexcept(false);

  /// Constructor - the trapezoid boundaries (arbitrary trapezoid)
  ///
  /// @param minhalex is the half length in x at minimal y
  /// @param haley is the half length in y
  /// @param halez is the half length in z
  /// @param alpha is the opening angle at -x,-y
  /// @param beta is the opening angle at +x,-y
  TrapezoidVolumeBounds(double minhalex, double haley, double halez,
                        double alpha, double beta) noexcept(false);

  /// Constructor - from a fixed size array
  ///
  /// @param values The bound values
  TrapezoidVolumeBounds(const std::array<double, eSize>& values) noexcept(false)
      : m_values(values) {
    checkConsistency();
    buildSurfaceBounds();
  }

  TrapezoidVolumeBounds(const TrapezoidVolumeBounds& trabo) = default;
  ~TrapezoidVolumeBounds() override = default;
  TrapezoidVolumeBounds& operator=(const TrapezoidVolumeBounds& trabo) =
      default;

  VolumeBounds::BoundsType type() const final {
    return VolumeBounds::eTrapezoid;
  }

  /// Return the bound values as dynamically sized vector
  ///
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

  /// This method checks if position in the 3D volume frame
  /// is inside the cylinder
  ///
  /// @param pos is the global position to be checked
  /// @param tol is the tolerance applied
  ///
  /// @return boolean indicator if position is inside
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

  /// Output Method for std::ostream
  /// @param os is the output stream
  std::ostream& toStream(std::ostream& os) const override;

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  double get(BoundValues bValue) const { return m_values[bValue]; }

 private:
  /// The internal version of the bounds can be float/double
  std::array<double, eSize> m_values{};
  /// The face PlaneSurface parallel to local xy plane
  std::shared_ptr<const TrapezoidBounds> m_faceXYTrapezoidBounds{nullptr};
  /// Thhe face PlaneSurface attached to alpha (negative local x)
  std::shared_ptr<const RectangleBounds> m_faceAlphaRectangleBounds{nullptr};
  /// The face PlaneSurface attached to beta (positive local x)
  std::shared_ptr<const RectangleBounds> m_faceBetaRectangleBounds{nullptr};
  /// The face PlaneSurface parallel to local zx plane, negative local y
  std::shared_ptr<const RectangleBounds> m_faceZXRectangleBoundsBottom{nullptr};
  /// The face PlaneSurface parallel to local zx plane, positive local y
  std::shared_ptr<const RectangleBounds> m_faceZXRectangleBoundsTop{nullptr};

  /// Check the input values for consistency,
  /// will throw a logic_exception if consistency is not given
  void checkConsistency() noexcept(false);

  /// Helper method to create the surface bounds
  void buildSurfaceBounds();
};

}  // namespace Acts
