// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Utilities/DelegateChainBuilder.hpp"

using namespace Acts;

struct AddTo {
  int value = 0;

  void add(int &x) const { x += value; }
};

void addFive(int &x) {
  x += 5;
}

BOOST_AUTO_TEST_SUITE(DelegateChainBuilderTests)

BOOST_AUTO_TEST_CASE(DelegateChainBuilderAdd) {
  AddTo a1{1}, a2{2}, a3{3};
  int x = 0;

  // Basic building
  OwningDelegate<void(int &)> chain{DelegateChainBuilder<void(int &)>{}
                                        .add<&AddTo::add>(&a1)
                                        .add<&addFive>()
                                        .add<&AddTo::add>(&a2)
                                        .add<&AddTo::add>(&a3)
                                        .build()};
  chain(x);
  BOOST_CHECK_EQUAL(x, 11);

  chain.disconnect();

  // In case of no return types, we can rebind the owning delegate with a chain
  // of different size
  chain = DelegateChainBuilder<void(int &)>{}
              .add<&AddTo::add>(&a1)
              .add<&addFive>()
              .add<&AddTo::add>(&a3)
              .build();

  x = 0;

  chain(x);
  BOOST_CHECK_EQUAL(x, 9);

  // CTAD helper from delegate type
  chain = DelegateChainBuilder{chain}
              .add<&AddTo::add>(&a1)
              .add<&addFive>()
              .add<&AddTo::add>(&a3)
              .build();

  x = 0;

  chain(x);
  BOOST_CHECK_EQUAL(x, 9);

  Delegate<void(int &)> nonOwning;

  // In case of a single callable, we can store it in a non-owning delegate
  DelegateChainBuilder<void(int &)>{}.add<&AddTo::add>(&a1).store(nonOwning);

  x = 0;
  nonOwning(x);
  BOOST_CHECK_EQUAL(x, 1);
}

struct GetInt {
  int value;

  int get() const { return value; }
};

int getSix() {
  return 6;
}

BOOST_AUTO_TEST_CASE(DelegateChainBuilderReturn) {
  GetInt g1{1}, g2{2}, g3{3};

  Delegate<std::array<int, 4>(), void, DelegateType::Owning> chain =
      DelegateChainBuilder<int()>{}
          .add<&GetInt::get>(&g1)
          .add<&getSix>()
          .add<&GetInt::get>(&g2)
          .add<&GetInt::get>(&g3)
          .build();

  auto results = chain();
  std::vector<int> expected = {1, 6, 2, 3};
  BOOST_CHECK_EQUAL_COLLECTIONS(results.begin(), results.end(),
                                expected.begin(), expected.end());

  Delegate<std::array<int, 3>(), void, DelegateType::Owning> delegate;
  DelegateChainBuilder<int()>{}
      .add<&GetInt::get>(&g1)
      .add<&getSix>()
      .add<&GetInt::get>(&g3)
      .store(delegate);

  auto results2 = delegate();
  expected = {1, 6, 3};
  BOOST_CHECK_EQUAL_COLLECTIONS(results2.begin(), results2.end(),
                                expected.begin(), expected.end());
}

BOOST_AUTO_TEST_SUITE_END()
