// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <type_traits>

namespace Acts {

/// The following operators have to be inplemented in order to satisfy
/// as an actor in the propagation
///
/// clang-format off
///
/// @code
///
/// template <typename propagator_state_t, typename result_t>
/// void
/// operator()(propagator_state_t& state, result_t& result)
/// const
/// {
///   return false;
/// }
///
/// template <typename propagator_state_t>
/// void
/// operator()(propagator_state_t& state)
/// const
/// {
///   return false;
/// }
///
/// @endcode
///
/// clang-format off
namespace detail {

namespace {

template <typename T, typename propagator_state_t, typename stepper_t,
          typename result_t,
          typename = decltype(std::declval<T>().operator()(
              std::declval<propagator_state_t&>(), std::declval<stepper_t&>(),
              std::declval<result_t&>()))>
std::true_type test_action_with_result(int);

template <typename, typename, typename, typename>
std::false_type test_action_with_result(...);

template <typename T, typename propagator_state_t, typename stepper_t,
          typename = decltype(std::declval<T>().operator()(
              std::declval<propagator_state_t&>(), std::declval<stepper_t&>()))>
std::true_type test_action_without_result(int);

template <typename>
std::false_type test_action_without_result(...);

template <typename T, typename propagator_state_t, typename stepper_t,
          bool has_result = false>
struct action_signature_check_impl
    : decltype(test_action_without_result<T, propagator_state_t, stepper_t>(
          0)) {};

template <typename T, typename propagator_state_t, typename stepper_t>
struct action_signature_check_impl<T, propagator_state_t, stepper_t, true>
    : decltype(test_action_with_result<T, propagator_state_t, stepper_t,
                                       detail::result_type_t<T>>(0)) {};

template <typename T, typename propagator_state_t, typename stepper_t>
struct action_signature_check
    : action_signature_check_impl<T, propagator_state_t, stepper_t,
                                  detail::has_result_type_v<T>> {};
}  // end of anonymous namespace

template <typename T, typename propagator_state_t, typename stepper_t>
constexpr bool action_signature_check_v =
    action_signature_check<T, propagator_state_t, stepper_t>::value;
}  // namespace detail

}  // namespace Acts