// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Concepts.hpp"

#include <cassert>
#include <functional>
#include <memory>
#include <type_traits>

namespace Acts {
/// Ownership enum for @c Delegate
enum class DelegateType { Owning, NonOwning };

template <auto C>
struct DelegateFuncTag {
  explicit constexpr DelegateFuncTag() = default;
};

// Specialization needed for defaulting ownership and for R(Args...) syntax
template <typename, typename H = void, DelegateType = DelegateType::NonOwning>
class Delegate;

/// Delegate type that allows type erasure of a callable without allocation
/// and with a single level of indirection. This type can support:
/// - a free function pointer
/// - a pointer to a member function alongside an instance pointer
/// @note @c Delegate by default does not assume ownership of the instance.
///          In that case You need to ensure that the lifetime of the callable
///          instance is longer than that of the @c Delegate. If you set @c O
///          to @c DelegateType::Owning, it will assume ownership.
/// @note Currently @c Delegate only supports callables that are ``const``
/// @tparam R Return type of the function signature
/// @tparam H Holder type that is used to store an instance
/// @tparam O Ownership type of the delegate: Owning or NonOwning
/// @tparam Args Types of the arguments of the function signatures
///
template <typename R, typename H, DelegateType O, typename... Args>
class Delegate<R(Args...), H, O> {
 public:
  static constexpr DelegateType kOwnership = O;

  /// Alias of the return type
  using return_type = R;
  using holder_type = H;
  using function_type = return_type (*)(const holder_type *, Args...);
  using function_ptr_type = return_type (*)(Args...);
  using signature_type = R(Args...);

  using deleter_type = void (*)(const holder_type *);

 private:
  template <typename T, typename C>
  using isSignatureCompatible =
      decltype(std::declval<T &>() = std::declval<C>());

  using OwningDelegate =
      Delegate<R(Args...), holder_type, DelegateType::Owning>;
  using NonOwningDelegate =
      Delegate<R(Args...), holder_type, DelegateType::NonOwning>;

  template <typename T>
  using isNoFunPtr = std::conjunction<
      std::negation<std::is_convertible<std::decay_t<T>, function_type>>,
      std::negation<std::is_same<std::decay_t<T>, OwningDelegate>>,
      std::negation<std::is_same<std::decay_t<T>, NonOwningDelegate>>>;

 public:
  Delegate() = default;

  Delegate(Delegate &&) noexcept = default;
  Delegate &operator=(Delegate &&) noexcept = default;
  Delegate(const Delegate &) noexcept = default;
  Delegate &operator=(const Delegate &) noexcept = default;

  /// @cond
  /// Constructor with an explicit runtime callable
  /// @param callable The runtime value of the callable
  /// @note The function signature requires the first argument of the callable is `const void*`.
  ///       i.e. if the signature of the delegate is `void(int)`, the
  ///       callable's signature has to be `void(const void*, int)`.
  explicit Delegate(function_type callable) { connect(callable); }
  /// @endcond

  /// Constructor with a possibly stateful function object.
  /// @tparam Callable Type of the callable
  /// @param callable The callable (function object or lambda)
  /// @note @c Delegate does not assume owner ship over @p callable. You need to ensure
  ///       it's lifetime is longer than that of @c Delegate.
  template <typename Callable>
  explicit Delegate(Callable &callable)
    requires(isNoFunPtr<Callable>::value)
  {
    connect(callable);
  }

  /// Constructor with a compile-time free function pointer
  /// @tparam Callable The compile-time free function pointer
  /// @note @c DelegateFuncTag is used to communicate the callable type
  template <auto Callable>
  explicit Delegate(DelegateFuncTag<Callable> /*tag*/) {
    connect<Callable>();
  }

  /// Constructor with a compile-time member function pointer and instance
  /// @tparam Callable The compile-time member function pointer
  /// @tparam Type The type of the instance the member function should be called on
  /// @param instance The instance on which the member function pointer should be called on
  /// @note @c Delegate does not assume owner ship over @p instance.
  /// @note @c DelegateFuncTag is used to communicate the callable type
  template <auto Callable, typename Type>

  Delegate(DelegateFuncTag<Callable> /*tag*/, const Type *instance)
    requires(kOwnership == DelegateType::NonOwning)
  {
    connect<Callable>(instance);
  }

  /// Constructor from rvalue reference is deleted, should catch construction
  /// with temporary objects and thus invalid pointers
  template <typename Callable>
  Delegate(Callable &&)
    requires(isNoFunPtr<Callable>::value)
  = delete;

  /// @cond
  /// Assignment operator with an explicit runtime callable
  /// @param callable The runtime value of the callable
  /// @note The function signature requires the first argument of the callable is `const void*`.
  ///       i.e. if the signature of the delegate is `void(int)`, the
  ///       callable's signature has to be `void(const void*, int)`.
  void operator=(function_type callable) { connect(callable); }

  /// Assignment operator with possibly stateful function object.
  /// @tparam Callable Type of the callable
  /// @param callable The callable (function object or lambda)
  /// @note @c Delegate does not assume owner ship over @p callable. You need to ensure
  ///       it's lifetime is longer than that of @c Delegate.
  template <typename Callable>
  void operator=(Callable &callable)
    requires(isNoFunPtr<Callable>::value)
  {
    connect(callable);
  }
  /// @endcond

  /// Assignment operator from rvalue reference is deleted, should catch
  /// assignment from temporary objects and thus invalid pointers
  template <typename Callable>
  void operator=(Callable &&)
    requires(isNoFunPtr<Callable>::value)
  = delete;

  /// @cond
  /// Connect a free function pointer.
  /// @note The function pointer must be ``constexpr`` for @c Delegate to accept it
  /// @tparam Callable The compile-time free function pointer
  template <auto Callable>
  void connect()
    requires(
        Concepts::invocable_and_returns<Callable, return_type, Args && ...>)
  {
    m_payload.payload = nullptr;

    static_assert(
        Concepts::invocable_and_returns<Callable, return_type, Args &&...>,
        "Callable given does not correspond exactly to required call "
        "signature");

    m_function = [](const holder_type * /*payload*/,
                    Args... args) -> return_type {
      return std::invoke(Callable, std::forward<Args>(args)...);
    };
  }
  /// @endcond

  /// Assignment operator with possibly stateful function object.
  /// @tparam Callable Type of the callable
  /// @param callable The callable (function object or lambda)
  /// @note @c Delegate does not assume owner ship over @p callable. You need to ensure
  ///       it's lifetime is longer than that of @c Delegate.
  template <typename Callable>
  void connect(Callable &callable)
    requires(isNoFunPtr<Callable>::value)
  {
    connect<&Callable::operator(), Callable>(&callable);
  }

  /// Connection with rvalue reference is deleted, should catch assignment
  /// from temporary objects and thus invalid pointers
  template <typename Callable>
  void connect(Callable &&)
    requires(isNoFunPtr<Callable>::value)
  = delete;

  /// @cond
  /// Connect anything that is assignable to the function pointer
  /// @param callable The runtime value of the callable
  /// @note The function signature requires the first argument of the callable is `const void*`.
  ///       i.e. if the signature of the delegate is `void(int)`, the
  ///       callable's signature has to be `void(const void*, int)`.
  void connect(function_type callable) {
    if constexpr (kOwnership == DelegateType::NonOwning) {
      m_payload.payload = nullptr;
    }
    m_function = callable;
  }
  /// @endcond

  template <typename Type>
  void connect(function_type callable, const Type *instance)
    requires(kOwnership == DelegateType::NonOwning)
  {
    m_payload.payload = instance;
    m_function = callable;
  }

  /// Connect a member function to be called on an instance
  /// @tparam Callable The compile-time member function pointer
  /// @tparam Type The type of the instance the member function should be called on
  /// @param instance The instance on which the member function pointer should be called on
  /// @note @c Delegate does not assume owner ship over @p instance. You need to ensure
  ///       it's lifetime is longer than that of @c Delegate.
  template <auto Callable, typename Type>
  void connect(const Type *instance)
    requires(kOwnership == DelegateType::NonOwning &&
             Concepts::invocable_and_returns<Callable, return_type, Type,
                                             Args && ...>)

  {
    m_payload.payload = instance;

    m_function = [](const holder_type *payload, Args... args) -> return_type {
      assert(payload != nullptr && "Payload is required, but not set");
      const auto *concretePayload = static_cast<const Type *>(payload);
      return std::invoke(Callable, concretePayload,
                         std::forward<Args>(args)...);
    };
  }

  /// Connect a member function to be called on an instance
  /// @tparam Callable The compile-time member function pointer
  /// @tparam Type The type of the instance the member function should be called on
  /// @param instance The instance on which the member function pointer should be called on
  /// @note @c Delegate assumes owner ship over @p instance.
  template <auto Callable, typename Type>
  void connect(std::unique_ptr<const Type> instance)
    requires(kOwnership == DelegateType::Owning &&
             Concepts::invocable_and_returns<Callable, return_type, Type,
                                             Args && ...>)
  {
    m_payload.payload = std::unique_ptr<const holder_type, deleter_type>(
        instance.release(), [](const holder_type *payload) {
          const auto *concretePayload = static_cast<const Type *>(payload);
          delete concretePayload;
        });

    m_function = [](const holder_type *payload, Args &&...args) -> return_type {
      assert(payload != nullptr && "Payload is required, but not set");
      const auto *concretePayload = static_cast<const Type *>(payload);
      return std::invoke(Callable, concretePayload,
                         std::forward<Args>(args)...);
    };
  }

  /// The call operator that exposes the functionality of the @c Delegate type.
  /// @param args The arguments to call the contained function with
  /// @return Return value of the contained function
  template <typename... Ts>
  return_type operator()(Ts &&...args) const
    requires(std::is_invocable_v<function_type, const holder_type *, Ts...>)
  {
    assert(connected() && "Delegate is not connected");
    return std::invoke(m_function, m_payload.ptr(), std::forward<Ts>(args)...);
  }

  /// Return whether this delegate is currently connected
  /// @return True if this delegate is connected
  bool connected() const { return m_function != nullptr; }

  /// Return whether this delegate is currently connected
  /// @return True if this delegate is connected
  explicit operator bool() const { return connected(); }

  /// Disconnect this delegate, meaning it cannot be called anymore
  void disconnect() {
    m_payload.clear();
    m_function = nullptr;
  }

  const holder_type *instance() const
    requires(!std::same_as<holder_type, void>)
  {
    return m_payload.ptr();
  }

 private:
  // Deleter that does not do anything
  static void noopDeleter(const holder_type * /*unused*/) {
    // we do not own the payload
  }

  /// @cond

  // Payload object without a deleter
  struct NonOwningPayload {
    void clear() { payload = nullptr; }

    const holder_type *ptr() const { return payload; }

    const holder_type *payload{nullptr};
  };

  // Payload object with a deleter
  struct OwningPayload {
    void clear() { payload.reset(); }

    const holder_type *ptr() const { return payload.get(); }

    std::unique_ptr<const holder_type, deleter_type> payload{nullptr,
                                                             &noopDeleter};
  };

  /// Stores the instance pointer and maybe a deleter
  std::conditional_t<kOwnership == DelegateType::NonOwning, NonOwningPayload,
                     OwningPayload>
      m_payload;

  /// @endcond

  /// Stores the function pointer wrapping the compile time function pointer given in @c connect().
  function_type m_function{nullptr};
};

template <typename, typename H = void>
class OwningDelegate;

/// Alias for an owning delegate
template <typename R, typename H, typename... Args>
class OwningDelegate<R(Args...), H>
    : public Delegate<R(Args...), H, DelegateType::Owning> {
 public:
  OwningDelegate() = default;
  explicit OwningDelegate(
      Delegate<R(Args...), H, DelegateType::Owning> &&delegate)
      : Delegate<R(Args...), H, DelegateType::Owning>(std::move(delegate)) {}

  OwningDelegate &operator=(
      Delegate<R(Args...), H, DelegateType::Owning> &&delegate) {
    *this = OwningDelegate{std::move(delegate)};
    return *this;
  }
};

}  // namespace Acts
