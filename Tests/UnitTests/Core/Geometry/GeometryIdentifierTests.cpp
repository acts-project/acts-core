// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryIdentifier.hpp"

namespace Acts::Test {

BOOST_AUTO_TEST_CASE(GeometryIdentifier_construct_default) {
  GeometryIdentifier id;
  BOOST_CHECK_EQUAL(id.volume(), 0u);
  BOOST_CHECK_EQUAL(id.boundary(), 0u);
  BOOST_CHECK_EQUAL(id.layer(), 0u);
  BOOST_CHECK_EQUAL(id.approach(), 0u);
  BOOST_CHECK_EQUAL(id.sensitive(), 0u);
}

BOOST_AUTO_TEST_CASE(GeometryIdentifier_construct_encoded) {
  // not sure if it is a good idea to test for the encoding since it should be
  // an implementation detail. only the resulting functionality is relevant.
  GeometryIdentifier id{0xa0b00c00d00affe0u};
  BOOST_CHECK_EQUAL(id.volume(), 0xa0u);
  BOOST_CHECK_EQUAL(id.boundary(), 0xb0u);
  BOOST_CHECK_EQUAL(id.layer(), 0x0c0u);
  BOOST_CHECK_EQUAL(id.approach(), 0x0du);
  BOOST_CHECK_EQUAL(id.sensitive(), 0x00affu);
  BOOST_CHECK_EQUAL(id.extra(), 0xe0u);
}

BOOST_AUTO_TEST_CASE(GeometryIdentifier_max_values) {
  // compute maximum value for each component
  constexpr GeometryIdentifier::Value volumeMax = (1u << 8) - 1;
  constexpr GeometryIdentifier::Value boundaryMax = (1u << 8) - 1;
  constexpr GeometryIdentifier::Value layerMax = (1u << 12) - 1;
  constexpr GeometryIdentifier::Value approachMax = (1u << 8) - 1;
  constexpr GeometryIdentifier::Value sensitiveMax = (1u << 20) - 1;
  constexpr GeometryIdentifier::Value extraMax = (1u << 8) - 1;
  // reference values non-zero values everywhere
  constexpr GeometryIdentifier ref{0xdeadaffe01234567};
  // values above the maximum are truncated
  // max+1 has all available bits zeroed
  BOOST_CHECK_EQUAL(GeometryIdentifier(ref).setVolume(volumeMax + 1),
                    GeometryIdentifier(ref).setVolume(0u));
  BOOST_CHECK_EQUAL(GeometryIdentifier(ref).setBoundary(boundaryMax + 1),
                    GeometryIdentifier(ref).setBoundary(0u));
  BOOST_CHECK_EQUAL(GeometryIdentifier(ref).setLayer(layerMax + 1),
                    GeometryIdentifier(ref).setLayer(0u));
  BOOST_CHECK_EQUAL(GeometryIdentifier(ref).setApproach(approachMax + 1),
                    GeometryIdentifier(ref).setApproach(0u));
  BOOST_CHECK_EQUAL(GeometryIdentifier(ref).setSensitive(sensitiveMax + 1),
                    GeometryIdentifier(ref).setSensitive(0u));
  BOOST_CHECK_EQUAL(GeometryIdentifier(ref).setExtra(extraMax + 1),
                    GeometryIdentifier(ref).setExtra(0u));
}

BOOST_AUTO_TEST_CASE(GeometryIdentifier_order) {
  auto vol1 = GeometryIdentifier()
                  .setVolume(1u)
                  .setLayer(14u)
                  .setSensitive(5u)
                  .setExtra(42u);
  auto vol2 = GeometryIdentifier()
                  .setVolume(2u)
                  .setLayer(13u)
                  .setSensitive(3u)
                  .setExtra(43u);
  // order uses volume first even if other components are larger
  BOOST_CHECK_LT(vol1, vol2);
  BOOST_CHECK_LT(GeometryIdentifier(vol1).setBoundary(64u), vol2);
  BOOST_CHECK_LT(GeometryIdentifier(vol1).setLayer(64u), vol2);
  BOOST_CHECK_LT(GeometryIdentifier(vol1).setApproach(64u), vol2);
  BOOST_CHECK_LT(GeometryIdentifier(vol1).setSensitive(64u), vol2);
  BOOST_CHECK_LT(GeometryIdentifier(vol1).setSensitive(64u), vol2);
  BOOST_CHECK_LT(vol2, GeometryIdentifier(vol1).setVolume(3u));
  // other components are hierarchical
  BOOST_CHECK_LT(GeometryIdentifier(vol1).setVolume(1u).setBoundary(2u),
                 GeometryIdentifier(vol1).setVolume(2u).setBoundary(1u));
  BOOST_CHECK_LT(GeometryIdentifier(vol1).setBoundary(1u).setLayer(2u),
                 GeometryIdentifier(vol1).setBoundary(2u).setLayer(1u));
  BOOST_CHECK_LT(GeometryIdentifier(vol1).setLayer(1u).setApproach(2u),
                 GeometryIdentifier(vol1).setLayer(2u).setApproach(1u));
  BOOST_CHECK_LT(GeometryIdentifier(vol1).setApproach(1u).setSensitive(2u),
                 GeometryIdentifier(vol1).setApproach(2u).setSensitive(1u));
  BOOST_CHECK_LT(GeometryIdentifier(vol1).setSensitive(1u).setExtra(2u),
                 GeometryIdentifier(vol1).setSensitive(2u).setExtra(1u));
}

}  // namespace Acts::Test
