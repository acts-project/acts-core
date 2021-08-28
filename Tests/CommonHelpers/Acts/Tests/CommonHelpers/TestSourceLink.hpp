// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <array>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <stdexcept>

namespace Acts {
namespace Test {

/// A minimal source link implementation for testing.
///
/// Instead of storing a reference to a measurement or raw data, the measurement
/// data is stored inline directly in the source link. Only 1d or 2d
/// measurements are supported to limit the overhead. Additionaly, a source
/// identifier is stored that can be used to store additional information. How
/// this is interpreted depends on the specific tests.
struct TestSourceLink {
  GeometryIdentifier geoId;
  size_t sourceId = 0u;
  // use eBoundSize to indicate unused indices
  std::array<BoundIndices, 2> indices = {eBoundSize, eBoundSize};
  Acts::ActsVector<2> parameters;
  Acts::ActsSymMatrix<2> covariance;

  /// Construct a source link for a 1d measurement.
  TestSourceLink(BoundIndices idx, ActsScalar val, ActsScalar var,
                 GeometryIdentifier gid = GeometryIdentifier(), size_t sid = 0u)
      : geoId(gid),
        sourceId(sid),
        indices{idx, eBoundSize},
        parameters(val, 0),
        covariance(Acts::ActsVector<2>(var, 0).asDiagonal()) {}
  /// Construct a source link for a 2d measurement.
  TestSourceLink(BoundIndices idx0, BoundIndices idx1,
                 const Acts::ActsVector<2>& params,
                 const Acts::ActsSymMatrix<2>& cov,
                 GeometryIdentifier gid = GeometryIdentifier(), size_t sid = 0u)
      : geoId(gid),
        sourceId(sid),
        indices{idx0, idx1},
        parameters(params),
        covariance(cov) {}
  /// Default-construct an invalid source link to satisfy SourceLinkConcept.
  TestSourceLink() = default;
  TestSourceLink(const TestSourceLink&) = default;
  TestSourceLink(TestSourceLink&&) = default;
  TestSourceLink& operator=(const TestSourceLink&) = default;
  TestSourceLink& operator=(TestSourceLink&&) = default;

  constexpr GeometryIdentifier geometryId() const { return geoId; }
  std::size_t index() const { return sourceId; }
};

bool operator==(const TestSourceLink& lhs, const TestSourceLink& rhs);
bool operator!=(const TestSourceLink& lhs, const TestSourceLink& rhs);
std::ostream& operator<<(std::ostream& os, const TestSourceLink& sourceLink);

/// Extract measurements from TestSourceLinks.
struct TestSourceLinkCalibrator {
  /// Extract the measurement from a TestSourceLink.
  ///
  /// @tparam parameters_t Track parameters type
  /// @param sourceLink Input source link
  /// @param parameters Input track parameters (unused)
  ///
  /// Since the TestSourceLink stores the necessary data inline, this just
  /// constructs the correct type from the stored data. Consequently, it does
  /// not depend on the track parameters, but they still must be part of the
  /// interface.
  template <typename parameters_t>
  BoundVariantMeasurement<TestSourceLink> operator()(
      const TestSourceLink& sl, const parameters_t& /* parameters */) const {
    if ((sl.indices[0] != eBoundSize) and (sl.indices[1] != eBoundSize)) {
      return makeMeasurement(sl, sl.parameters, sl.covariance, sl.indices[0],
                             sl.indices[1]);
    } else if (sl.indices[0] != eBoundSize) {
      return makeMeasurement(sl, sl.parameters.head<1>(),
                             sl.covariance.topLeftCorner<1, 1>(),
                             sl.indices[0]);
    } else {
      throw std::runtime_error(
          "Tried to extract measurement from invalid TestSourceLink");
    }
  }
};

}  // namespace Test
}  // namespace Acts
