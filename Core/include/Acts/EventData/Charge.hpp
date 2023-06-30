// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"

#include <cassert>
#include <cmath>

namespace Acts {

/// @defgroup eventdata-charge Charge interpretation for track parameters
///
/// Track parameters store a single coefficient that describes charge and
/// momentum. This is either charge/momentum or 1/momentum, but the
/// interpretation depends on what type of particle is described. In this code
/// base this coefficient is always referred to as `qOverP` (or
/// charge-over-momentum) even for uncharged particles. The following types are
/// used to restrict the particle charge magnitude (at compile time) and support
/// the umambigous extraction of charge and absolute momentum from said track
/// parameter coefficient.
///
/// All types are designed to be interchangeable. Each one can be
/// constructed with the input charge magnitude
///
/// ```cpp
/// Charge c(1_e);
/// ```
///
/// and can then be used to extract the charge value
///
/// ```cpp
/// auto q = c.extractCharge(qOverP);
/// ```
///
/// or the absolute momentum
///
/// ```cpp
/// auto p = c.extractMomentum(qOverP);
/// ```
///
/// from the charge-over-momentum track parameter.
///
/// @{

/// Charge and momentum interpretation for neutral particles.
struct Neutral {
  static constexpr bool canHoldNeutral = true;
  static constexpr bool canHoldSinglyCharged = false;
  static constexpr bool canHoldMultiCharged = false;

  Neutral() = default;
  /// Construct and verify the input charge magnitude (in debug builds).
  ///
  /// This constructor is only provided to allow consistent construction.
  constexpr Neutral(float absQ) noexcept {
    assert((absQ == 0) and "Input charge must be zero");
    (void)absQ;
  }

  constexpr float absQ() const noexcept { return 0; }

  template <typename T>
  constexpr auto qFromQOP(T /*qOverP*/) const noexcept {
    return 0.0f;
  }
  template <typename T>
  constexpr auto pFromQOP(T qOverP) const noexcept {
    assert(qOverP >= 0 && "qOverP cannot be negative");
    return 1.0f / qOverP;
  }

  template <typename P, typename Q>
  constexpr auto qopFromPQ(P p, Q q) const noexcept {
    assert((q != 0) and "charge must be 0");
    (void)q;
    return 1.0f / p;
  }

  /// Compare for equality.
  ///
  /// This is always `true` as `Neutral` has no internal state.
  /// Must be available to provide a consistent interface.
  friend constexpr bool operator==(Neutral /*lhs*/, Neutral /*rhs*/) noexcept {
    return true;
  }
};

/// Charge and momentum interpretation for particles with +-e charge.
struct SinglyCharged {
  static constexpr bool canHoldNeutral = false;
  static constexpr bool canHoldSinglyCharged = true;
  static constexpr bool canHoldMultiCharged = false;

  SinglyCharged() = default;
  /// Construct and verify the input charge magnitude (in debug builds).
  ///
  /// This constructor is only provided to allow consistent construction.
  constexpr SinglyCharged(float absQ) noexcept {
    assert((absQ == UnitConstants::e) and "Input charge magnitude must be e");
    (void)absQ;
  }

  constexpr float absQ() const noexcept { return UnitConstants::e; }

  template <typename T>
  constexpr auto qFromQOP(T qOverP) const noexcept {
    // using because of autodiff
    using std::copysign;
    return copysign(UnitConstants::e, qOverP);
  }
  template <typename T>
  constexpr auto pFromQOP(T qOverP) const noexcept {
    // using because of autodiff
    using std::abs;
    return qFromQOP(qOverP) / qOverP;
  }

  template <typename P, typename Q>
  constexpr auto qopFromPQ(P p, Q q) const noexcept {
    using std::abs;
    assert((abs(q) == UnitConstants::e) && "absolute charge must be e");
    return q / p;
  }

  /// Compare for equality.
  ///
  /// This is always `true` as `SinglyCharged` has no internal state.
  /// Must be available to provide a consistent interface.
  friend constexpr bool operator==(SinglyCharged /*lhs*/,
                                   SinglyCharged /*rhs*/) noexcept {
    return true;
  }
};

/// Charge and momentum interpretation for arbitrarily charged but not neutral
/// particles.
class NonNeutralCharge {
 public:
  static constexpr bool canHoldNeutral = false;
  static constexpr bool canHoldSinglyCharged = true;
  static constexpr bool canHoldMultiCharged = true;

  /// Construct with the magnitude of the input charge.
  constexpr NonNeutralCharge(float absQ) noexcept : m_absQ{absQ} {
    assert((0 < absQ) and "Input charge magnitude must be positive");
  }
  constexpr NonNeutralCharge(SinglyCharged) noexcept
      : m_absQ{UnitConstants::e} {}

  constexpr float absQ() const noexcept { return m_absQ; }

  template <typename T>
  constexpr auto qFromQOP(T qOverP) const noexcept {
    // using because of autodiff
    using std::copysign;
    return copysign(m_absQ, qOverP);
  }
  template <typename T>
  constexpr auto pFromQOP(T qOverP) const noexcept {
    // using because of autodiff
    using std::abs;
    return qFromQOP(qOverP) / qOverP;
  }

  template <typename P, typename Q>
  constexpr auto qopFromPQ(P p, Q q) const noexcept {
    // using because of autodiff
    using std::abs;
    assert(abs(q) == m_absQ && "inconsistent charge");
    return q / p;
  }

  /// Compare for equality.
  friend constexpr bool operator==(NonNeutralCharge lhs,
                                   NonNeutralCharge rhs) noexcept {
    return lhs.m_absQ == rhs.m_absQ;
  }

 private:
  float m_absQ{};
};

/// Charge and momentum interpretation for arbitrarily charged particles.
///
/// Only a charge magnitude identical to zero is interpreted as representing a
/// neutral particle. This avoids ambiguities that might arise from using an
/// approximate comparison with an arbitrary epsilon.
class AnyCharge {
 public:
  static constexpr bool canHoldNeutral = true;
  static constexpr bool canHoldSinglyCharged = true;
  static constexpr bool canHoldMultiCharged = true;

  /// Construct with the magnitude of the input charge.
  constexpr AnyCharge(float absQ) noexcept : m_absQ{absQ} {
    assert((0 <= absQ) and "Input charge magnitude must be zero or positive");
  }
  constexpr AnyCharge(SinglyCharged) noexcept : m_absQ{UnitConstants::e} {}
  constexpr AnyCharge(Neutral) noexcept : m_absQ{0} {}

  constexpr float absQ() const noexcept { return m_absQ; }

  template <typename T>
  constexpr auto qFromQOP(T qOverP) const noexcept {
    // using because of autodiff
    using std::copysign;
    return copysign(m_absQ, qOverP);
  }
  template <typename T>
  constexpr auto pFromQOP(T qOverP) const noexcept {
    // using because of autodiff
    using std::abs;
    return (m_absQ != 0.0f) ? qFromQOP(qOverP) / qOverP : 1.0f / qOverP;
  }

  template <typename P, typename Q>
  constexpr auto qopFromPQ(P p, Q q) const noexcept {
    // using because of autodiff
    using std::abs;
    assert(abs(q) == m_absQ && "inconsistent charge");
    return (m_absQ != 0.0f) ? q / p : 1.0f / p;
  }

  /// Compare for equality.
  friend constexpr bool operator==(AnyCharge lhs, AnyCharge rhs) noexcept {
    return lhs.m_absQ == rhs.m_absQ;
  }

 private:
  float m_absQ{};
};

/// @}

}  // namespace Acts
