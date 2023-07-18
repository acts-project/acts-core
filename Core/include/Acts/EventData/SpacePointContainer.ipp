// -*- C++ -*-
// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <math.h>

namespace Acts {

template <typename container_t, template <typename> class holder_t>
template <template <typename> class, typename>
SpacePointContainer<container_t, holder_t>::SpacePointContainer(
    const Acts::SpacePointContainerConfig& config,
    const Acts::SpacePointContainerOptions& options, container_t& container)
    : m_config(config.toInternalUnits()),
      m_options(options.toInternalUnits()),
      m_container(container) {
  this->initialize();
}

template <typename container_t, template <typename> class holder_t>
template <template <typename> class, typename>
SpacePointContainer<container_t, holder_t>::SpacePointContainer(
    const Acts::SpacePointContainerConfig& config,
    const Acts::SpacePointContainerOptions& options, container_t&& container)
    : m_config(config.toInternalUnits()),
      m_options(options.toInternalUnits()),
      m_container(std::move(container)) {
  this->initialize();
}

template <typename container_t, template <typename> class holder_t>
void SpacePointContainer<container_t, holder_t>::initialize() {
  m_data.resize(this->size(), m_config.useDetailedDoubleMeasurementInfo);
  m_proxies.reserve(this->size());
  const auto& external_container = this->container();
  for (std::size_t i(0); i < this->size(); ++i) {
    m_data.setX(i, external_container.x_impl(i) - m_options.beamPos[0]);
    m_data.setY(i, external_container.y_impl(i) - m_options.beamPos[1]);
    m_data.setZ(i, external_container.z_impl(i));
    m_data.setRadius(
        i, std::sqrt(m_data.x(i) * m_data.x(i) + m_data.y(i) * m_data.y(i)));
    m_data.setPhi(i, atan2f(m_data.y(i), m_data.x(i)));
    m_data.setVarianceR(i, external_container.varianceR_impl(i));
    m_data.setVarianceZ(i, external_container.varianceZ_impl(i));

    m_proxies.emplace_back(*this, i);
  }

  // Dynamic variables
  if (m_config.useDetailedDoubleMeasurementInfo) {
    using namespace Acts::HashedStringLiteral;
    for (std::size_t i(0); i < this->size(); ++i) {
      m_data.setTopStripVector(
          i, std::any_cast<Acts::Vector3>(
                 external_container.component_impl("TopStripVector"_hash, i)));
      m_data.setBottomStripVector(
          i, std::any_cast<Acts::Vector3>(external_container.component_impl(
                 "BottomStripVector"_hash, i)));
      m_data.setStripCenterDistance(
          i, std::any_cast<Acts::Vector3>(external_container.component_impl(
                 "StripCenterDistance"_hash, i)));
      m_data.setTopStripCenterPosition(
          i, std::any_cast<Acts::Vector3>(external_container.component_impl(
                 "TopStripCenterPosition"_hash, i)));
    }
  }
}

template <typename container_t, template <typename> class holder_t>
template <template <typename> class, typename>
SpacePointContainer<container_t, holder_t>::SpacePointContainer(
    SpacePointContainer<container_t, holder_t>& other)
    : m_config(other.m_config),
      m_options(other.m_options),
      m_container(*other.m_container.ptr) {}

template <typename container_t, template <typename> class holder_t>
template <template <typename> class, typename>
SpacePointContainer<container_t, holder_t>&
SpacePointContainer<container_t, holder_t>::operator=(
    SpacePointContainer<container_t, holder_t>& other) {
  m_config = other.m_config;
  m_options = other.m_options;
  m_container.ptr = other.m_container.ptr;
  return *this;
}

template <typename container_t, template <typename> class holder_t>
SpacePointContainer<container_t, holder_t>::SpacePointContainer(
    SpacePointContainer<container_t, holder_t>&& other) noexcept
    : m_config(other.m_config),
      m_options(other.m_options),
      m_container(std::exchange(other.m_container.ptr, nullptr)) {}

template <typename container_t, template <typename> class holder_t>
inline const Acts::Vector3&
SpacePointContainer<container_t, holder_t>::topStripVector(
    const std::size_t n) const {
  return m_data.topStripVector(n);
}

template <typename container_t, template <typename> class holder_t>
inline const Acts::Vector3&
SpacePointContainer<container_t, holder_t>::bottomStripVector(
    const std::size_t n) const {
  return m_data.bottomStripVector(n);
}

template <typename container_t, template <typename> class holder_t>
inline const Acts::Vector3&
SpacePointContainer<container_t, holder_t>::stripCenterDistance(
    const std::size_t n) const {
  return m_data.stripCenterDistance(n);
}

template <typename container_t, template <typename> class holder_t>
inline const Acts::Vector3&
SpacePointContainer<container_t, holder_t>::topStripCenterPosition(
    const std::size_t n) const {
  return m_data.topStripCenterPosition(n);
}

template <typename container_t, template <typename> class holder_t>
SpacePointContainer<container_t, holder_t>&
SpacePointContainer<container_t, holder_t>::operator=(
    SpacePointContainer<container_t, holder_t>&& other) noexcept {
  m_container = std::exchange(other.m_container.ptr, nullptr);
  return *this;
}

template <typename container_t, template <typename> class holder_t>
inline std::size_t SpacePointContainer<container_t, holder_t>::size() const {
  return container().size_impl();
}

template <typename container_t, template <typename> class holder_t>
template <bool, typename>
inline typename SpacePointContainer<container_t, holder_t>::iterator
SpacePointContainer<container_t, holder_t>::begin() {
  return {*this, 0};
}

template <typename container_t, template <typename> class holder_t>
template <bool, typename>
inline typename SpacePointContainer<container_t, holder_t>::iterator
SpacePointContainer<container_t, holder_t>::end() {
  return {*this, size()};
}

template <typename container_t, template <typename> class holder_t>
inline typename SpacePointContainer<container_t, holder_t>::const_iterator
SpacePointContainer<container_t, holder_t>::begin() const {
  return {*this, 0};
}

template <typename container_t, template <typename> class holder_t>
inline typename SpacePointContainer<container_t, holder_t>::const_iterator
SpacePointContainer<container_t, holder_t>::end() const {
  return {*this, size()};
}

template <typename container_t, template <typename> class holder_t>
template <bool, typename>
inline container_t& SpacePointContainer<container_t, holder_t>::container() {
  return *m_container;
}

template <typename container_t, template <typename> class holder_t>
inline const container_t&
SpacePointContainer<container_t, holder_t>::container() const {
  return *m_container;
}

template <typename container_t, template <typename> class holder_t>
template <bool, typename>
inline typename SpacePointContainer<container_t, holder_t>::ValueType&
SpacePointContainer<container_t, holder_t>::sp(const std::size_t n) {
  return container().get_impl(n);
}

template <typename container_t, template <typename> class holder_t>
inline typename SpacePointContainer<container_t, holder_t>::ValueType&
SpacePointContainer<container_t, holder_t>::sp(const std::size_t n) const {
  return container().get_impl(n);
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::x(
    const std::size_t n) const {
  return m_data.x(n);
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::y(
    const std::size_t n) const {
  return m_data.y(n);
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::z(
    const std::size_t n) const {
  return m_data.z(n);
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::phi(
    const std::size_t n) const {
  return m_data.phi(n);
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::radius(
    const std::size_t n) const {
  return m_data.radius(n);
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::varianceR(
    const std::size_t n) const {
  return m_data.varianceR(n);
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::varianceZ(
    const std::size_t n) const {
  return m_data.varianceZ(n);
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::quality(
    const std::size_t n) const {
  return m_data.quality(n);
}

template <typename container_t, template <typename> class holder_t>
inline float SpacePointContainer<container_t, holder_t>::deltaR(
    const std::size_t n) const {
  return m_data.deltaR(n);
}

template <typename container_t, template <typename> class holder_t>
inline void SpacePointContainer<container_t, holder_t>::setQuality(
    const std::size_t n, const float value) const {
  m_data.setQuality(n, value);
}

template <typename container_t, template <typename> class holder_t>
inline void SpacePointContainer<container_t, holder_t>::setDeltaR(
    const std::size_t n, const float value) const {
  m_data.setDeltaR(n, value);
}

template <typename container_t, template <typename> class holder_t>
template <bool, typename>
inline typename SpacePointContainer<container_t, holder_t>::ProxyType&
SpacePointContainer<container_t, holder_t>::proxy(const std::size_t n) {
  return proxies()[n];
}

template <typename container_t, template <typename> class holder_t>
inline const typename SpacePointContainer<container_t, holder_t>::ProxyType&
SpacePointContainer<container_t, holder_t>::proxy(const std::size_t n) const {
  return proxies()[n];
}

template <typename container_t, template <typename> class holder_t>
const std::vector<
    typename SpacePointContainer<container_t, holder_t>::ProxyType>&
SpacePointContainer<container_t, holder_t>::proxies() const {
  return m_proxies;
}

template <typename container_t, template <typename> class holder_t>
template <bool, typename>
std::vector<typename SpacePointContainer<container_t, holder_t>::ProxyType>&
SpacePointContainer<container_t, holder_t>::proxies() {
  return m_proxies;
}

}  // namespace Acts
