// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Utilities/Logger.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <unordered_map>
#include <vector>

namespace ActsExamples {

/// A container to store arbitrary objects with ownership transfer.
///
/// This is an append-only container that takes ownership of the objects
/// added to it. Once an object has been added, it can only be read but not
/// be modified. Trying to replace an existing object is considered an error.
/// Its lifetime is bound to the liftime of the white board.
class WhiteBoard {
 public:
  WhiteBoard(std::unique_ptr<const Acts::Logger> logger =
                 Acts::getDefaultLogger("WhiteBoard", Acts::Logging::INFO));

  // A WhiteBoard holds unique elements and can not be copied
  WhiteBoard(const WhiteBoard& other) = delete;
  WhiteBoard& operator=(const WhiteBoard&) = delete;

  /// Store an object on the white board and transfer ownership.
  ///
  /// @param name Non-empty identifier to store it under
  /// @param object Movable reference to the transferable object
  /// @throws std::invalid_argument on empty or duplicate name
  template <typename T>
  void add(const std::string& name, T&& object);

  /// Get access to a stored object.
  ///
  /// @param[in] name Identifier for the object
  /// @return reference to the stored object
  /// @throws std::out_of_range if no object is stored under the requested name
  template <typename T>
  const T& get(const std::string& name) const;

  bool exists(const std::string& name) const;

 private:
  // type-erased value holder for move-constructible types
  struct IHolder {
    virtual ~IHolder() = default;
    virtual const std::type_info& type() const = 0;
  };
  template <typename T,
            typename =
                std::enable_if_t<std::is_nothrow_move_constructible<T>::value>>
  struct HolderT : public IHolder {
    T value;

    HolderT(T&& v) : value(std::move(v)) {}
    const std::type_info& type() const { return typeid(T); }
  };

  std::unique_ptr<const Acts::Logger> m_logger;
  std::unordered_map<std::string, std::unique_ptr<IHolder>> m_store;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples

inline ActsExamples::WhiteBoard::WhiteBoard(
    std::unique_ptr<const Acts::Logger> logger)
    : m_logger(std::move(logger)) {}

template <typename T>
inline void ActsExamples::WhiteBoard::add(const std::string& name, T&& object) {
  if (name.empty()) {
    throw std::invalid_argument("Object can not have an empty name");
  }
  if (0 < m_store.count(name)) {
    throw std::invalid_argument("Object '" + name + "' already exists");
  }
  m_store.emplace(name, std::make_unique<HolderT<T>>(std::forward<T>(object)));
  ACTS_VERBOSE("Added object '" << name << "'");
}

template <typename T>
inline const T& ActsExamples::WhiteBoard::get(const std::string& name) const {
  auto it = m_store.find(name);
  if (it == m_store.end()) {
    throw std::out_of_range("Object '" + name + "' does not exists");
  }
  const IHolder* holder = it->second.get();
  if (typeid(T) != holder->type()) {
    throw std::out_of_range("Type mismatch for object '" + name + "'");
  }
  ACTS_VERBOSE("Retrieved object '" << name << "'");
  return reinterpret_cast<const HolderT<T>*>(holder)->value;
}

inline bool ActsExamples::WhiteBoard::exists(const std::string& name) const {
  return m_store.find(name) != m_store.end();
}