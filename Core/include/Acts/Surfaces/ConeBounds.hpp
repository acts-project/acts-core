// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iosfwd>
#include <numbers>
#include <vector>

namespace Acts {

///  @class ConeBounds
///
///  Bounds for a conical surface,
///  the opening angle is stored in \f$ \tan(\alpha) \f$ and always positively
/// defined.
///  The cone can open to both sides, steered by \f$ z_min \f$ and \f$ z_max
///  \f$.
///
///  @image html ConeBounds.gif
///
class ConeBounds : public SurfaceBounds {
 public:
  enum BoundValues : int {
    eAlpha = 0,
    eMinZ = 1,
    eMaxZ = 2,
    eHalfPhiSector = 3,
    eAveragePhi = 4,
    eSize = 5
  };

  /// Constructor - open cone with alpha, by default a full cone
  /// but optionally can make a conical section
  ///
  /// @param alpha is the opening angle of the cone
  /// @param symm is the boolean indicating if the cone is symmetric in +/- z
  /// @param halfphi is the half opening angle (default is pi)
  /// @param avphi is the phi value around which the bounds are opened
  /// (default=0)
  ConeBounds(double alpha, bool symm, double halfphi = std::numbers::pi,
             double avphi = 0.) noexcept(false);

  /// Constructor - open cone with alpha, minz and maxz, by
  /// default a full cone but can optionally make it a conical section
  ///
  /// @param alpha is the opening angle of the cone
  /// @param minz cone expanding from minimal z
  /// @param maxz cone expanding to maximal z
  /// @param halfphi is the half opening angle (default is pi)
  /// @param avphi is the phi value around which the bounds are opened
  /// (default=0)
  ConeBounds(double alpha, double minz, double maxz,
             double halfphi = std::numbers::pi,
             double avphi = 0.) noexcept(false);

  /// Constructor - from parameters array
  ///
  /// @param values The parameter array
  ConeBounds(const std::array<double, eSize>& values) noexcept(false);

  /// @copydoc SurfaceBounds::type
  BoundsType type() const final { return eCone; }

  /// @copydoc SurfaceBounds::isCartesian
  bool isCartesian() const final { return true; }

  /// @copydoc SurfaceBounds::boundToCartesianJacobian
  SquareMatrix2 boundToCartesianJacobian(const Vector2& lposition) const final {
    (void)lposition;
    return SquareMatrix2::Identity();
  }

  /// @copydoc SurfaceBounds::boundToCartesianMetric
  SquareMatrix2 boundToCartesianMetric(const Vector2& lposition) const final {
    (void)lposition;
    return SquareMatrix2::Identity();
  }

  /// @copydoc SurfaceBounds::values
  std::vector<double> values() const final;

  /// @copydoc SurfaceBounds::inside
  bool inside(const Vector2& lposition) const final;

  /// @copydoc SurfaceBounds::closestPoint
  Vector2 closestPoint(const Vector2& lposition,
                       const std::optional<SquareMatrix2>& metric) const final;

  using SurfaceBounds::inside;

  /// Output Method for std::ostream
  /// @param sl is the ostrea into which the dump is done
  /// @return is the input object
  std::ostream& toStream(std::ostream& sl) const final;

  /// Return the radius at a specific z values
  ///
  /// @param z is the z value for which r is requested
  /// @return is the r value associated with z
  double r(double z) const { return std::abs(z * m_tanAlpha); }

  /// Return tangent of alpha (pre-computed)
  double tanAlpha() const { return m_tanAlpha; }

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  double get(BoundValues bValue) const { return m_values[bValue]; }

 private:
  std::array<double, eSize> m_values;
  double m_tanAlpha;

  /// Check the input values for consistency, will throw a logic_exception
  /// if consistency is not given
  void checkConsistency() noexcept(false);

  /// Private helper function to shift a local 2D position
  ///
  /// Shift r-phi coordinate to be centered around the average phi.
  ///
  /// @param lposition The original local position
  Vector2 shifted(const Vector2& lposition) const;
};

}  // namespace Acts
