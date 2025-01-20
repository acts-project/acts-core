// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryBackendConcept.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/EventData/detail/DynamicColumn.hpp"
#include "Acts/EventData/detail/DynamicKeyIterator.hpp"
#include "Acts/Utilities/EigenConcepts.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

#include <any>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/histogram.hpp>

namespace Acts {
class Surface;
template <typename T>
struct IsReadOnlyMultiTrajectory;

namespace detail_vmt {

using MultiTrajectoryTraits::IndexType;
constexpr auto kInvalid = MultiTrajectoryTraits::kInvalid;
constexpr auto MeasurementSizeMax = MultiTrajectoryTraits::MeasurementSizeMax;

template <typename T>
struct NonInitializingAllocator {
  using value_type = T;

  NonInitializingAllocator() noexcept = default;

  template <class U>
  explicit NonInitializingAllocator(
      const NonInitializingAllocator<U>& /*other*/) noexcept {}

  template <class U>
  bool operator==(const NonInitializingAllocator<U>& /*other*/) const noexcept {
    return true;
  }

  T* allocate(std::size_t n) const { return std::allocator<T>{}.allocate(n); }

  void deallocate(T* const p, std::size_t n) const noexcept {
    std::allocator<T>{}.deallocate(p, n);
  }

  void construct(T* /*p*/) const {
    // This construct function intentionally does not initialize the object!
    // Be very careful when using this allocator.
  }
};

class VectorMultiTrajectoryBase {
 public:
  struct Statistics {
    using axis_t = boost::histogram::axis::variant<
        boost::histogram::axis::category<std::string>,
        boost::histogram::axis::category<>>;

    using axes_t = std::vector<axis_t>;
    using hist_t = boost::histogram::histogram<axes_t>;

    hist_t hist;

    void toStream(std::ostream& os, std::size_t n = 1);
  };

  template <typename T>
  Statistics statistics(T& instance) const {
    using namespace boost::histogram;
    using cat = axis::category<std::string>;

    Statistics::axes_t axes;
    axes.emplace_back(cat({
        "count",
        "index",
        "parPred",
        "covPred",
        "parFilt",
        "covFilt",
        "parSmth",
        "covSmth",
        "meas",
        "measOffset",
        "measCov",
        "measCovOffset",
        "jac",
        "sourceLinks",
        "projectors",
    }));

    axes.emplace_back(axis::category<>({0, 1}));

    auto h = make_histogram(axes);

    for (IndexType i = 0; i < instance.size(); i++) {
      auto ts = instance.getTrackState(i);

      bool isMeas = ts.typeFlags().test(TrackStateFlag::MeasurementFlag);

      h("count", isMeas);

      h("index", isMeas, weight(sizeof(IndexData)));

      using scalar = typename decltype(ts.predicted())::Scalar;
      std::size_t par_size = eBoundSize * sizeof(scalar);
      std::size_t cov_size = eBoundSize * eBoundSize * sizeof(scalar);

      const IndexData& index = m_index[i];
      if (ts.hasPredicted() &&
          ACTS_CHECK_BIT(index.allocMask, TrackStatePropMask::Predicted)) {
        h("parPred", isMeas, weight(par_size));
        h("covPred", isMeas, weight(cov_size));
      }
      if (ts.hasFiltered() &&
          ACTS_CHECK_BIT(index.allocMask, TrackStatePropMask::Filtered)) {
        h("parFilt", isMeas, weight(par_size));
        h("covFilt", isMeas, weight(cov_size));
      }
      if (ts.hasSmoothed() &&
          ACTS_CHECK_BIT(index.allocMask, TrackStatePropMask::Smoothed)) {
        h("parSmth", isMeas, weight(par_size));
        h("covSmth", isMeas, weight(cov_size));
      }
      h("sourceLinks", isMeas, weight(sizeof(SourceLink)));
      h("measOffset", isMeas,
        weight(sizeof(decltype(m_measOffset)::value_type)));
      h("measCovOffset", isMeas,
        weight(sizeof(decltype(m_measCovOffset)::value_type)));
      if (ts.hasCalibrated() &&
          ACTS_CHECK_BIT(index.allocMask, TrackStatePropMask::Calibrated)) {
        std::size_t meas_size = ts.calibratedSize() * sizeof(scalar);
        std::size_t meas_cov_size =
            ts.calibratedSize() * ts.calibratedSize() * sizeof(scalar);

        h("meas", isMeas, weight(meas_size));
        h("measCov", isMeas, weight(meas_cov_size));
        h("sourceLinks", isMeas, weight(sizeof(const SourceLink)));
        h("projectors", isMeas, weight(sizeof(SerializedSubspaceIndices)));
      }

      if (ts.hasJacobian() &&
          ACTS_CHECK_BIT(index.allocMask, TrackStatePropMask::Jacobian)) {
        h("jac", isMeas, weight(cov_size));
      }
    }

    return Statistics{h};
  }

 protected:
  struct IndexData {
    IndexType ipredicted = kInvalid;
    IndexType ifiltered = kInvalid;
    IndexType ismoothed = kInvalid;
    IndexType ijacobian = kInvalid;
    IndexType iprojector = kInvalid;

    float chi2 = 0;
    double pathLength = 0;
    TrackStateType::raw_type typeFlags{};

    IndexType iUncalibrated = kInvalid;
    IndexType iCalibratedSourceLink = kInvalid;
    IndexType measdim = kInvalid;

    TrackStatePropMask allocMask = TrackStatePropMask::None;
  };

  VectorMultiTrajectoryBase() noexcept = default;

  VectorMultiTrajectoryBase(const VectorMultiTrajectoryBase& other)
      : m_index{other.m_index},
        m_previous{other.m_previous},
        m_next{other.m_next},
        m_params{other.m_params},
        m_cov{other.m_cov},
        m_meas{other.m_meas},
        m_measOffset{other.m_measOffset},
        m_measCov{other.m_measCov},
        m_measCovOffset{other.m_measCovOffset},
        m_jac{other.m_jac},
        m_sourceLinks{other.m_sourceLinks},
        m_projectors{other.m_projectors},
        m_referenceSurfaces{other.m_referenceSurfaces} {
    for (const auto& [key, value] : other.m_dynamic) {
      m_dynamic.insert({key, value->clone()});
    }
    m_dynamicKeys = other.m_dynamicKeys;
  };

  VectorMultiTrajectoryBase(VectorMultiTrajectoryBase&& other) = default;

  // BEGIN INTERFACE HELPER
  template <typename T>
  static constexpr bool has_impl(T& instance, HashedString key,
                                 IndexType istate) {
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      case "predicted"_hash:
        return instance.m_index[istate].ipredicted != kInvalid;
      case "filtered"_hash:
        return instance.m_index[istate].ifiltered != kInvalid;
      case "smoothed"_hash:
        return instance.m_index[istate].ismoothed != kInvalid;
      case "calibrated"_hash:
        return instance.m_measOffset[istate] != kInvalid;
      case "calibratedCov"_hash:
        return instance.m_measCovOffset[istate] != kInvalid;
      case "jacobian"_hash:
        return instance.m_index[istate].ijacobian != kInvalid;
      case "projector"_hash:
        return instance.m_index[istate].iprojector != kInvalid;
      case "uncalibratedSourceLink"_hash:
        return instance.m_sourceLinks[instance.m_index[istate].iUncalibrated]
            .has_value();
      case "previous"_hash:
      case "next"_hash:
      case "referenceSurface"_hash:
      case "measdim"_hash:
      case "chi2"_hash:
      case "pathLength"_hash:
      case "typeFlags"_hash:
        return true;
      default:
        return instance.m_dynamic.contains(key);
    }
  }

  template <bool EnsureConst, typename T>
  static std::any component_impl(T& instance, HashedString key,
                                 IndexType istate) {
    if constexpr (EnsureConst) {
      static_assert(std::is_const_v<std::remove_reference_t<T>>,
                    "Is not const");
    }
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      case "previous"_hash:
        return &instance.m_previous[istate];
      case "next"_hash:
        return &instance.m_next[istate];
      case "predicted"_hash:
        return &instance.m_index[istate].ipredicted;
      case "filtered"_hash:
        return &instance.m_index[istate].ifiltered;
      case "smoothed"_hash:
        return &instance.m_index[istate].ismoothed;
      case "projector"_hash:
        return &instance.m_projectors[instance.m_index[istate].iprojector];
      case "measdim"_hash:
        return &instance.m_index[istate].measdim;
      case "chi2"_hash:
        return &instance.m_index[istate].chi2;
      case "pathLength"_hash:
        return &instance.m_index[istate].pathLength;
      case "typeFlags"_hash:
        return &instance.m_index[istate].typeFlags;
      default:
        auto it = instance.m_dynamic.find(key);
        if (it == instance.m_dynamic.end()) {
          throw std::runtime_error("Unable to handle this component");
        }
        std::conditional_t<EnsureConst, const detail::DynamicColumnBase*,
                           detail::DynamicColumnBase*>
            col = it->second.get();
        assert(col && "Dynamic column is null");
        return col->get(istate);
    }
  }

  template <typename T>
  static bool hasColumn_impl(T& instance, HashedString key) {
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      case "predicted"_hash:
      case "filtered"_hash:
      case "smoothed"_hash:
      case "calibrated"_hash:
      case "calibratedCov"_hash:
      case "jacobian"_hash:
      case "projector"_hash:
      case "previous"_hash:
      case "next"_hash:
      case "uncalibratedSourceLink"_hash:
      case "referenceSurface"_hash:
      case "measdim"_hash:
      case "chi2"_hash:
      case "pathLength"_hash:
      case "typeFlags"_hash:
        return true;
      default:
        return instance.m_dynamic.contains(key);
    }
  }

 public:
  detail::DynamicKeyRange<detail::DynamicColumnBase> dynamicKeys_impl() const {
    return {m_dynamic.begin(), m_dynamic.end()};
  }

  // END INTERFACE HELPER

 public:
  IndexType calibratedSize_impl(IndexType istate) const {
    return m_index[istate].measdim;
  }

  SourceLink getUncalibratedSourceLink_impl(IndexType istate) const {
    return m_sourceLinks[m_index[istate].iUncalibrated].value();
  }

  const Surface* referenceSurface_impl(IndexType istate) const {
    return m_referenceSurfaces[istate].get();
  }

 protected:
  /// index to map track states to the corresponding
  std::vector<IndexData> m_index;
  std::vector<IndexType> m_previous;
  std::vector<IndexType> m_next;
  std::vector<typename detail_lt::FixedSizeTypes<eBoundSize>::Coefficients>
      m_params;
  std::vector<typename detail_lt::FixedSizeTypes<eBoundSize>::Covariance> m_cov;

  std::vector<double, NonInitializingAllocator<double>> m_meas;
  std::vector<MultiTrajectoryTraits::IndexType> m_measOffset;
  std::vector<double, NonInitializingAllocator<double>> m_measCov;
  std::vector<MultiTrajectoryTraits::IndexType> m_measCovOffset;

  std::vector<typename detail_lt::FixedSizeTypes<eBoundSize>::Covariance> m_jac;
  std::vector<std::optional<SourceLink>> m_sourceLinks;
  std::vector<SerializedSubspaceIndices> m_projectors;

  // owning vector of shared pointers to surfaces
  //
  // This might be problematic when appending a large number of surfaces
  // trackstates, because vector has to reallocated and thus copy. This might
  // be handled in a smart way by moving but not sure.
  std::vector<std::shared_ptr<const Surface>> m_referenceSurfaces;

  std::vector<HashedString> m_dynamicKeys;
  std::unordered_map<HashedString, std::unique_ptr<detail::DynamicColumnBase>>
      m_dynamic;
};

}  // namespace detail_vmt

class VectorMultiTrajectory;

template <>
struct IsReadOnlyMultiTrajectory<VectorMultiTrajectory> : std::false_type {};

class VectorMultiTrajectory final
    : public detail_vmt::VectorMultiTrajectoryBase,
      public MultiTrajectory<VectorMultiTrajectory> {
#ifndef DOXYGEN
  friend MultiTrajectory<VectorMultiTrajectory>;
#endif

 public:
  VectorMultiTrajectory() = default;
  VectorMultiTrajectory(const VectorMultiTrajectory& other)
      : VectorMultiTrajectoryBase{other} {}

  VectorMultiTrajectory(VectorMultiTrajectory&& other) noexcept
      : VectorMultiTrajectoryBase{std::move(other)} {}

  Statistics statistics() const {
    return detail_vmt::VectorMultiTrajectoryBase::statistics(*this);
  }

  // BEGIN INTERFACE
  TrackStateProxy::Parameters parameters_impl(IndexType parIdx) {
    return TrackStateProxy::Parameters{m_params[parIdx].data()};
  }

  ConstTrackStateProxy::Parameters parameters_impl(IndexType parIdx) const {
    return ConstTrackStateProxy::Parameters{m_params[parIdx].data()};
  }

  TrackStateProxy::Covariance covariance_impl(IndexType parIdx) {
    return TrackStateProxy::Covariance{m_cov[parIdx].data()};
  }

  ConstTrackStateProxy::Covariance covariance_impl(IndexType parIdx) const {
    return ConstTrackStateProxy::Covariance{m_cov[parIdx].data()};
  }

  TrackStateProxy::Covariance jacobian_impl(IndexType istate) {
    IndexType jacIdx = m_index[istate].ijacobian;
    return TrackStateProxy::Covariance{m_jac[jacIdx].data()};
  }

  ConstTrackStateProxy::Covariance jacobian_impl(IndexType istate) const {
    IndexType jacIdx = m_index[istate].ijacobian;
    return ConstTrackStateProxy::Covariance{m_jac[jacIdx].data()};
  }

  template <std::size_t measdim>
  TrackStateProxy::Calibrated<measdim> calibrated_impl(IndexType istate) {
    IndexType offset = m_measOffset[istate];
    return TrackStateProxy::Calibrated<measdim>{&m_meas[offset]};
  }

  template <std::size_t measdim>
  ConstTrackStateProxy::Calibrated<measdim> calibrated_impl(
      IndexType istate) const {
    IndexType offset = m_measOffset[istate];
    return ConstTrackStateProxy::Calibrated<measdim>{&m_meas[offset]};
  }

  template <std::size_t measdim>
  TrackStateProxy::CalibratedCovariance<measdim> calibratedCovariance_impl(
      IndexType istate) {
    IndexType offset = m_measCovOffset[istate];
    return TrackStateProxy::CalibratedCovariance<measdim>{&m_measCov[offset]};
  }

  template <std::size_t measdim>
  ConstTrackStateProxy::CalibratedCovariance<measdim> calibratedCovariance_impl(
      IndexType istate) const {
    IndexType offset = m_measCovOffset[istate];
    return ConstTrackStateProxy::CalibratedCovariance<measdim>{
        &m_measCov[offset]};
  }

  IndexType addTrackState_impl(
      TrackStatePropMask mask = TrackStatePropMask::All,
      IndexType iprevious = kInvalid);

  void addTrackStateComponents_impl(IndexType istate, TrackStatePropMask mask);

  void reserve(std::size_t n);

  void shareFrom_impl(IndexType iself, IndexType iother,
                      TrackStatePropMask shareSource,
                      TrackStatePropMask shareTarget);

  void unset_impl(TrackStatePropMask target, IndexType istate);

  bool has_impl(HashedString key, IndexType istate) const {
    return detail_vmt::VectorMultiTrajectoryBase::has_impl(*this, key, istate);
  }

  IndexType size_impl() const { return m_index.size(); }

  void clear_impl();

  std::any component_impl(HashedString key, IndexType istate) {
    return detail_vmt::VectorMultiTrajectoryBase::component_impl<false>(
        *this, key, istate);
  }

  std::any component_impl(HashedString key, IndexType istate) const {
    return detail_vmt::VectorMultiTrajectoryBase::component_impl<true>(
        *this, key, istate);
  }

  template <typename T>
  void addColumn_impl(std::string_view key) {
    HashedString hashedKey = hashStringDynamic(key);
    m_dynamic.insert({hashedKey, std::make_unique<detail::DynamicColumn<T>>()});
  }

  bool hasColumn_impl(HashedString key) const {
    return detail_vmt::VectorMultiTrajectoryBase::hasColumn_impl(*this, key);
  }

  template <typename val_t, typename cov_t>
  void allocateCalibrated_impl(IndexType istate,
                               const Eigen::DenseBase<val_t>& val,
                               const Eigen::DenseBase<cov_t>& cov)
    requires(Concepts::eigen_base_is_fixed_size<val_t> &&
             Concepts::eigen_bases_have_same_num_rows<val_t, cov_t> &&
             Concepts::eigen_base_is_square<cov_t> &&
             Eigen::PlainObjectBase<val_t>::RowsAtCompileTime <=
                 toUnderlying(eBoundSize))
  {
    constexpr std::size_t measdim = val_t::RowsAtCompileTime;

    if (m_index[istate].measdim != kInvalid &&
        m_index[istate].measdim != measdim) {
      throw std::invalid_argument{
          "Measurement dimension does not match the allocated dimension"};
    }

    if (m_measOffset[istate] == kInvalid ||
        m_measCovOffset[istate] == kInvalid) {
      m_measOffset[istate] = static_cast<IndexType>(m_meas.size());
      m_meas.resize(m_meas.size() + measdim);

      m_measCovOffset[istate] = static_cast<IndexType>(m_measCov.size());
      m_measCov.resize(m_measCov.size() + measdim * measdim);
    }

    m_index[istate].measdim = measdim;

    double* measPtr = &m_meas[m_measOffset[istate]];
    Eigen::Map<ActsVector<measdim>> valMap(measPtr);
    valMap = val;

    double* covPtr = &m_measCov[m_measCovOffset[istate]];
    Eigen::Map<ActsSquareMatrix<measdim>> covMap(covPtr);
    covMap = cov;
  }

  void setUncalibratedSourceLink_impl(IndexType istate,
                                      SourceLink&& sourceLink) {
    m_sourceLinks[m_index[istate].iUncalibrated] = std::move(sourceLink);
  }

  void setReferenceSurface_impl(IndexType istate,
                                std::shared_ptr<const Surface> surface) {
    m_referenceSurfaces[istate] = std::move(surface);
  }

  void copyDynamicFrom_impl(IndexType dstIdx, HashedString key,
                            const std::any& srcPtr);

  // END INTERFACE
};

static_assert(
    MutableMultiTrajectoryBackend<VectorMultiTrajectory>,
    "VectorMultiTrajectory does not fulfill MutableMultiTrajectoryBackend");

class ConstVectorMultiTrajectory;

template <>
struct IsReadOnlyMultiTrajectory<ConstVectorMultiTrajectory> : std::true_type {
};

class ConstVectorMultiTrajectory final
    : public detail_vmt::VectorMultiTrajectoryBase,
      public MultiTrajectory<ConstVectorMultiTrajectory> {
#ifndef DOXYGEN
  friend MultiTrajectory<ConstVectorMultiTrajectory>;
#endif

 public:
  ConstVectorMultiTrajectory() = default;

  ConstVectorMultiTrajectory(const ConstVectorMultiTrajectory& other)
      : VectorMultiTrajectoryBase{other} {}

  ConstVectorMultiTrajectory(const VectorMultiTrajectory& other)
      : VectorMultiTrajectoryBase{other} {}

  ConstVectorMultiTrajectory(VectorMultiTrajectory&& other)
      : VectorMultiTrajectoryBase{std::move(other)} {}

  ConstVectorMultiTrajectory(ConstVectorMultiTrajectory&&) = default;

  Statistics statistics() const {
    return detail_vmt::VectorMultiTrajectoryBase::statistics(*this);
  }

  // BEGIN INTERFACE

  ConstTrackStateProxy::Parameters parameters_impl(IndexType parIdx) const {
    return ConstTrackStateProxy::Parameters{m_params[parIdx].data()};
  }

  ConstTrackStateProxy::Covariance covariance_impl(IndexType parIdx) const {
    return ConstTrackStateProxy::Covariance{m_cov[parIdx].data()};
  }

  ConstTrackStateProxy::Covariance jacobian_impl(IndexType istate) const {
    IndexType jacIdx = m_index[istate].ijacobian;
    return ConstTrackStateProxy::Covariance{m_jac[jacIdx].data()};
  }

  template <std::size_t measdim>
  ConstTrackStateProxy::Calibrated<measdim> calibrated_impl(
      IndexType istate) const {
    IndexType offset = m_measOffset[istate];
    return ConstTrackStateProxy::Calibrated<measdim>{&m_meas[offset]};
  }

  template <std::size_t measdim>
  ConstTrackStateProxy::CalibratedCovariance<measdim> calibratedCovariance_impl(
      IndexType istate) const {
    IndexType offset = m_measCovOffset[istate];
    return ConstTrackStateProxy::CalibratedCovariance<measdim>{
        &m_measCov[offset]};
  }

  bool has_impl(HashedString key, IndexType istate) const {
    return detail_vmt::VectorMultiTrajectoryBase::has_impl(*this, key, istate);
  }

  IndexType size_impl() const { return m_index.size(); }

  std::any component_impl(HashedString key, IndexType istate) const {
    return detail_vmt::VectorMultiTrajectoryBase::component_impl<true>(
        *this, key, istate);
  }

  bool hasColumn_impl(HashedString key) const {
    return detail_vmt::VectorMultiTrajectoryBase::hasColumn_impl(*this, key);
  }

  // END INTERFACE
};

static_assert(
    ConstMultiTrajectoryBackend<ConstVectorMultiTrajectory>,
    "ConctVectorMultiTrajectory does not fulfill ConstMultiTrajectoryBackend");

}  // namespace Acts
