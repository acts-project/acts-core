// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iosfwd>
#include <numbers>
#include <vector>

namespace Acts {

/// @class DiscTrapezoidBounds
///
/// Class to describe the bounds for a planar DiscSurface.
/// By providing an argument for hphisec, the bounds can
/// be restricted to a phi-range around the center position.
///
class DiscTrapezoidBounds : public DiscBounds {
 public:
  enum BoundValues : int {
    eHalfLengthXminR = 0,
    eHalfLengthXmaxR = 1,
    eMinR = 2,
    eMaxR = 3,
    eAveragePhi = 4,
    eStereo = 5,
    eSize = 6
  };

  /// Constructor for a symmetric Trapezoid giving min X length, max X length,
  /// Rmin and R max
  /// @param halfXminR half length in X at min radius
  /// @param halfXmaxR half length in X at maximum radius
  /// @param minR inner radius
  /// @param maxR outer radius
  /// @param avgPhi average phi value
  /// @param stereo optional stero angle applied
  DiscTrapezoidBounds(double halfXminR, double halfXmaxR, double minR,
                      double maxR, double avgPhi = std::numbers::pi / 2.,
                      double stereo = 0.) noexcept(false);

  /// Constructor - from fixed size array
  ///
  /// @param values The parameter values
  DiscTrapezoidBounds(const std::array<double, eSize>& values) noexcept(false)
      : m_values(values) {
    checkConsistency();
  }

  /// @copydoc SurfaceBounds::type
  BoundsType type() const final { return eDiscTrapezoid; }

  /// @copydoc SurfaceBounds::isCartesian
  bool isCartesian() const final { return false; }

  /// @copydoc SurfaceBounds::boundToCartesianJacobian
  SquareMatrix2 boundToCartesianJacobian(const Vector2& lposition) const final;

  /// @copydoc SurfaceBounds::boundToCartesianMetric
  SquareMatrix2 boundToCartesianMetric(const Vector2& lposition) const final;

  /// Return the bound values as dynamically sized vector
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

  /// @copydoc SurfaceBounds::inside
  bool inside(const Vector2& lposition) const final;

  /// @copydoc SurfaceBounds::closestPoint
  Vector2 closestPoint(const Vector2& lposition,
                       const std::optional<SquareMatrix2>& metric) const final;

  using SurfaceBounds::inside;

  /// Output Method for std::ostream
  std::ostream& toStream(std::ostream& sl) const final;

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  double get(BoundValues bValue) const { return m_values[bValue]; }

  /// This method returns inner radius
  double rMin() const final { return get(eMinR); }

  /// This method returns outer radius
  double rMax() const final { return get(eMaxR); }

  /// This method returns the center radius
  double rCenter() const {
    double rmin = get(eMinR);
    double rmax = get(eMaxR);
    double hxmin = get(eHalfLengthXminR);
    double hxmax = get(eHalfLengthXmaxR);
    auto hmin = std::sqrt(rmin * rmin - hxmin * hxmin);
    auto hmax = std::sqrt(rmax * rmax - hxmax * hxmax);
    return 0.5 * (hmin + hmax);
  }

  /// This method returns the stereo angle
  double stereo() const { return get(eStereo); }

  /// This method returns the halfPhiSector which is covered by the disc
  double halfPhiSector() const {
    auto minHalfPhi = std::asin(get(eHalfLengthXminR) / get(eMinR));
    auto maxHalfPhi = std::asin(get(eHalfLengthXmaxR) / get(eMaxR));
    return std::max(minHalfPhi, maxHalfPhi);
  }

  /// This method returns the half length in Y (this is Rmax -Rmin)
  double halfLengthY() const {
    double rmin = get(eMinR);
    double rmax = get(eMaxR);
    double hxmin = get(eHalfLengthXminR);
    double hxmax = get(eHalfLengthXmaxR);
    auto hmin = std::sqrt(rmin * rmin - hxmin * hxmin);
    auto hmax = std::sqrt(rmax * rmax - hxmax * hxmax);
    return 0.5 * (hmax - hmin);
  }

  /// Returns true for full phi coverage - obviously false here
  bool coversFullAzimuth() const final { return false; }

  /// Checks if this is inside the radial coverage
  /// given the a tolerance
  bool insideRadialBounds(double R, double tolerance = 0.) const final {
    return (R + tolerance > get(eMinR) && R - tolerance < get(eMaxR));
  }

  /// Return a reference radius for binning
  double binningValueR() const final { return 0.5 * (get(eMinR) + get(eMaxR)); }

  /// Return a reference phi for binning
  double binningValuePhi() const final { return get(eAveragePhi); }

  /// This method returns the xy coordinates of the four corners of the
  /// bounds in module coorindates (in xy)
  ///
  /// @param ignoredSegments is an ignored parameter only used for
  /// curved bound segments
  ///
  /// @return vector for vertices in 2D
  std::vector<Vector2> vertices(unsigned int ignoredSegments = 0u) const final;

 private:
  std::array<double, eSize> m_values;

  /// Dreived maximum y value
  double m_ymax = 0;

  /// Check the input values for consistency, will throw a logic_exception
  /// if consistency is not given
  void checkConsistency() noexcept(false);

  /// Private helper method to convert a local position
  /// into its Cartesian representation
  ///
  /// @param lposition The local position in polar coordinates
  Vector2 toLocalCartesian(const Vector2& lposition) const;
};

}  // namespace Acts
