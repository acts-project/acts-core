// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <concepts>
#include <optional>
#include <sstream>
#include <system_error>
#include <type_traits>
#include <utility>
#include <variant>

namespace Acts {

/// Class which encapsulates either a valid result, or an error
/// @tparam T The valid result value
/// @tparam E The error, defaults to `std::error_code`
///
template <typename T, typename E = std::error_code>
class Result {
  /// Private constructor which accepts an external variant.
  /// This is used by the factory static methods to set up
  /// the variant unambiguously in all cases.
  explicit Result(std::variant<T, E>&& var) : m_var(std::move(var)) {}

 public:
  using ValueType = T;
  using ErrorType = E;

  /// Default construction is disallowed.
  Result() = delete;

  /// Copy construction is disallowed
  Result(const Result<T, E>& other) = delete;

  /// Assignment is disallowed
  Result<T, E>& operator=(const Result<T, E>& other) = delete;

  /// Move construction is allowed
  Result(Result<T, E>&& other) noexcept : m_var(std::move(other.m_var)) {}

  /// Move assignment is allowed
  /// @param other The other result instance, rvalue reference
  /// @return The assigned instance
  Result<T, E>& operator=(Result<T, E>&& other) noexcept {
    m_var = std::move(other.m_var);
    return *this;
  }

  /// @brief Constructor from arbitrary value
  /// This constructor allows construction from any value. This constructor is
  /// only enabled if T and E are unambiguous, meaning they cannot be implicitly
  /// converted and there is T cannot be constructed from E and vice-versa.
  /// This means that when this is invoked, the value can be propagated to the
  /// underlying variant, and the assignment will be correct, and error will be
  /// an error, and a value will be a value.
  /// @note If T and E are ambiguous, use the `success` and `failure` static
  /// factory methods.
  /// @tparam T2 Type of the potential assignment
  /// @param value The potential value, could be an actual valid value or an
  /// error.
  template <typename T2>
  Result(T2 value) noexcept  // NOLINT <- suppress clang-tidy warning.
                             // Conversion here is crucial for ergonomics
    requires(!std::same_as<T, E> && !std::constructible_from<T, E> &&
             !std::convertible_to<T, E> && !std::constructible_from<E, T> &&
             !std::convertible_to<E, T> &&
             !(std::convertible_to<T2, T> && std::convertible_to<T2, E>))
      : m_var(std::conditional_t<std::is_convertible_v<T2, T>, T, E>{
            std::move(value)}) {}

  /// @brief Assignment operator from arbitrary value
  /// This operator allows construction from any value. The same rules as for
  /// the `Result(T2 value)` constructor apply.
  /// * @tparam T2 Type of the potential assignment
  /// @param value The potential value, could be an actual valid value or an
  /// error.
  /// @return The assigned instance
  template <typename T2>
  Result<T, E>& operator=(T2 value) noexcept
    requires(!std::same_as<T, E> && !std::constructible_from<T, E> &&
             !std::convertible_to<T, E> && !std::constructible_from<E, T> &&
             !std::convertible_to<E, T> &&
             !(std::convertible_to<T2, T> && std::convertible_to<T2, E>))
  {
    m_var = std::move(std::conditional_t<std::is_convertible_v<T2, T>, T, E>{
        std::move(value)});
    return *this;
  }

  /// Static helper factory which forces assignment as valid value.
  /// @param value The valid value to assign. Will not be converted to E.
  /// @return Initialized result object
  static Result<T, E> success(T value) {
    return Result<T, E>(
        std::variant<T, E>{std::in_place_index<0>, std::move(value)});
  }

  /// Static helper factory which forces assignment as an error.
  /// @param error The error to assign. Will not be converted to T.
  /// @return Initialized result object
  static Result<T, E> failure(E error) {
    return Result<T, E>(
        std::variant<T, E>{std::in_place_index<1>, std::move(error)});
  }

  /// Checks whether this result contains a valid value, and no error.
  /// @return bool Whether result contains an error or not.
  bool ok() const noexcept { return m_var.index() == 0; }

  /// Returns a reference into the variant to the valid value.
  /// @note If `!res.ok()`, this method will abort (noexcept)
  /// @return Reference to value stored in the variant.
  T& operator*() noexcept { return std::get<T>(m_var); }

  /// Returns a reference into the variant to the valid value.
  /// @note If `!res.ok()`, this method will abort (noexcept)
  /// @return Reference to value stored in the variant.
  const T& operator*() const noexcept { return std::get<T>(m_var); }

  /// Allows to access members of the stored object with `res->foo`
  /// similar to `std::optional`.
  /// @note If `!res.ok()`, this method will abort (noexcept)
  /// @return Pointer to value stored in the variant.
  T* operator->() noexcept { return &std::get<T>(m_var); }

  /// Allows to access members of the stored object with `res->foo`
  /// similar to `std::optional`.
  /// @note If `!res.ok()`, this method will abort (noexcept)
  /// @return Pointer to value stored in the variant.
  const T* operator->() const noexcept { return &std::get<T>(m_var); }

  /// Returns a reference to the error stored in the result.
  /// @note If `res.ok()` this method will abort (noexcept)
  /// @return Reference to the error
  E& error() & noexcept { return std::get<E>(m_var); }

  /// Returns a reference to the error stored in the result.
  /// @note If `res.ok()` this method will abort (noexcept)
  /// @return Reference to the error
  const E& error() const& noexcept { return std::get<E>(m_var); }

  /// Returns the error by-value.
  /// @note If `res.ok()` this method will abort (noexcept)
  /// @return The error
  E error() && noexcept { return std::move(std::get<E>(m_var)); }

  /// Retrieves the valid value from the result object.
  /// @note This is the lvalue version, returns a reference to the value
  /// @return The valid value as a reference
  T& value() & {
    checkValueAccess();
    return std::get<T>(m_var);
  }

  /// Retrieves the valid value from the result object.
  /// @note This is the lvalue version, returns a reference to the value
  /// @return The valid value as a reference
  const T& value() const& {
    checkValueAccess();
    return std::get<T>(m_var);
  }

  /// Retrieves the valid value from the result object.
  /// @note This is the rvalue version, returns the value
  /// by-value and moves out of the variant.
  /// @return The valid value by value, moved out of the variant.
  T value() && {
    checkValueAccess();
    return std::move(std::get<T>(m_var));
  }

  /// Retrieves the valid value from the result object, or returns a default
  /// value if no valid value exists.
  ///
  /// @param[in] v The default value to use if no valid value exists.
  /// @note This is the lvalue version.
  /// @note This function always returns by value.
  /// @return Either the valid value, or the given substitute.
  template <typename U>
  std::conditional_t<std::is_reference_v<U>, const T&, T> value_or(U&& v) const&
    requires(std::same_as<std::decay_t<U>, T>)
  {
    if (ok()) {
      return value();
    } else {
      return std::forward<U>(v);
    }
  }

  /// Retrieves the valid value from the result object, or returns a default
  /// value if no valid value exists.
  ///
  /// @param[in] v The default value to use if no valid value exists.
  /// @note This is the rvalue version which moves the value out.
  /// @note This function always returns by value.
  /// @return Either the valid value, or the given substitute.
  template <typename U>
  T value_or(U&& v) &&
    requires(std::same_as<std::decay_t<U>, T>)
  {
    if (ok()) {
      return std::move(*this).value();
    } else {
      return std::forward<U>(v);
    }
  }

  /// Transforms the value contained in this result.
  ///
  /// Applying a function `f` to a valid value `x` returns `f(x)`, while
  /// applying `f` to an invalid value returns another invalid value.
  ///
  /// @param[in] callable The transformation function to apply.
  /// @note This is the lvalue version.
  /// @note This functions is `fmap` on the functor in `A` of `Result<A, E>`.
  /// @return The modified valid value if exists, or an error otherwise.
  template <typename C>
  auto transform(C&& callable) const&
    requires std::invocable<C, const T&>
  {
    using CallableReturnType = decltype(std::declval<C>()(std::declval<T>()));
    using R = Result<std::decay_t<CallableReturnType>, E>;
    if (ok()) {
      return R::success(callable(value()));
    } else {
      return R::failure(error());
    }
  }

  /// Transforms the value contained in this result.
  ///
  /// Applying a function `f` to a valid value `x` returns `f(x)`, while
  /// applying `f` to an invalid value returns another invalid value.
  ///
  /// @param[in] callable The transformation function to apply.
  /// @note This is the rvalue version.
  /// @note This functions is `fmap` on the functor in `A` of `Result<A, E>`.
  /// @return The modified valid value if exists, or an error otherwise.
  template <typename C>
  auto transform(C&& callable) &&
    requires std::invocable<C, T&&>
  {
    using CallableReturnType = decltype(std::declval<C>()(std::declval<T>()));
    using R = Result<std::decay_t<CallableReturnType>, E>;
    if (ok()) {
      return R::success(callable(std::move(*this).value()));
    } else {
      return R::failure(std::move(*this).error());
    }
  }

  /// Bind a function to this result monadically.
  ///
  /// This function takes a function `f` and, if this result contains a valid
  /// value `x`, returns `f(x)`. If the type of `x` is `T`, then `f` is
  /// expected to accept type `T` and return `Result<U>`. In this case,
  /// `transform` would return the unhelpful type `Result<Result<U>>`, so
  /// `and_then` strips away the outer layer to return `Result<U>`. If the
  /// value is invalid, this returns an invalid value in `Result<U>`.
  ///
  /// @param[in] callable The transformation function to apply.
  /// @note This is the lvalue version.
  /// @note This functions is `>>=` on the functor in `A` of `Result<A, E>`.
  /// @return The modified valid value if exists, or an error otherwise.
  template <typename C>
  auto and_then(C&& callable) const&
    requires std::invocable<C, const T&>
  {
    using R = decltype(std::declval<C>()(std::declval<T>()));

    static_assert(std::same_as<typename R::ErrorType, ErrorType>,
                  "bind must take a callable with the same error type");

    if (ok()) {
      return callable(value());
    } else {
      return R::failure(error());
    }
  }

  /// Bind a function to this result monadically.
  ///
  /// This function takes a function `f` and, if this result contains a valid
  /// value `x`, returns `f(x)`. If the type of `x` is `T`, then `f` is
  /// expected to accept type `T` and return `Result<U>`. In this case,
  /// `transform` would return the unhelpful type `Result<Result<U>>`, so
  /// `and_then` strips away the outer layer to return `Result<U>`. If the
  /// value is invalid, this returns an invalid value in `Result<U>`.
  ///
  /// @param[in] callable The transformation function to apply.
  /// @note This is the rvalue version.
  /// @note This functions is `>>=` on the functor in `A` of `Result<A, E>`.
  /// @return The modified valid value if exists, or an error otherwise.
  template <typename C>
  auto and_then(C&& callable) &&
    requires std::invocable<C, T&&>
  {
    using R = decltype(std::declval<C>()(std::declval<T>()));

    static_assert(std::same_as<typename R::ErrorType, ErrorType>,
                  "bind must take a callable with the same error type");

    if (ok()) {
      return callable(std::move(*this).value());
    } else {
      return R::failure(std::move(*this).error());
    }
  }

 private:
  std::variant<T, E> m_var;

  void checkValueAccess() const {
    if (m_var.index() != 0) {
      if constexpr (std::is_same_v<E, std::error_code>) {
        std::stringstream ss;
        const auto& e = std::get<E>(m_var);
        ss << "Value called on error value: " << e.category().name() << ": "
           << e.message() << " [" << e.value() << "]";
        throw std::runtime_error(ss.str());
      } else {
        throw std::runtime_error("Value called on error value");
      }
    }
  }
};

/// Template specialization for the void case.
/// This specialization handles the case where there is no actual return value,
/// but
/// an error might be returned. Returning the error directly would make handling
/// different from other functions using the `Result<T, E>` mechanism.
/// `Result<void, E>` does not have the dereference operator, and value methods.
/// The static `success` factory does not accept a value.
/// @note To ease usage, this `Result<void, E>` is default constructible in the
/// *ok*
/// state, whereas `Result<T, E>` is not.
/// @tparam E The type of the error
///
template <typename E>
class Result<void, E> {
 public:
  /// Default constructor which initializes the result in the ok state.
  Result() = default;

  /// The copy constructor is deleted.
  Result(const Result<void, E>& other) = default;

  /// The (self) assignment operator is deleted.
  Result<void, E>& operator=(const Result<void, E>& other) = default;

  /// Move constructor
  /// @param other The other result object, rvalue ref
  Result(Result<void, E>&& other) noexcept : m_opt(std::move(other.m_opt)) {}

  /// Move assignment operator
  /// @param other The other result object, rvalue ref
  Result<void, E>& operator=(Result<void, E>&& other) noexcept {
    m_opt = std::move(other.m_opt);
    return *this;
  }

  /// Constructor from error. This implicitly requires E2 to be convertible to
  /// E.
  /// @tparam E2 The type of the actual error
  /// @param error The instance of the actual error
  template <typename E2>
  Result(E2 error) noexcept  // NOLINT <- suppress clang-tidy warning.
                             // Conversion here is crucial for ergonomics
      : m_opt(std::move(error)) {}

  /// Assignment operator from an error.
  /// @tparam E2 The type of the actual error
  /// @param error The instance of the actual error
  /// @return The assigned instance
  template <typename E2>
  Result<void, E>& operator=(E2 error) {
    m_opt = std::move(error);
    return *this;
  }

  /// Static factory function to initialize the result in the ok state.
  /// @return Result object, in ok state
  static Result<void, E> success() { return Result<void, E>(); }

  /// Static factory function to initialize the result in the error state.
  /// @param error The error to initialize with.
  /// @return Result object, in error state.
  static Result<void, E> failure(E error) {
    return Result<void, E>(std::move(error));
  }

  /// Checks whether this result is in the ok state, and no error.
  /// @return bool Whether result contains an error or not.
  bool ok() const noexcept { return !m_opt; }

  /// Returns a reference to the error stored in the result.
  /// @note If `res.ok()` this method will abort (noexcept)
  /// @return Reference to the error
  E& error() & noexcept { return m_opt.value(); }

  /// Returns a reference to the error stored in the result.
  /// @note If `res.ok()` this method will abort (noexcept)
  /// @return Reference to the error
  const E& error() const& noexcept { return m_opt.value(); }

  /// Returns the error by-value.
  /// @note If `res.ok()` this method will abort (noexcept)
  /// @return Reference to the error
  E error() && noexcept { return std::move(m_opt.value()); }

 private:
  std::optional<E> m_opt;
};

}  // namespace Acts
