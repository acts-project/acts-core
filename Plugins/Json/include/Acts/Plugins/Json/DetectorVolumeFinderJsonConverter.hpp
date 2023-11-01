// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/DetectorVolumeUpdators.hpp"
#include "Acts/Plugins/Json/DetrayJsonHelper.hpp"
#include "Acts/Plugins/Json/IndexedGridJsonHelper.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <tuple>
#include <vector>

namespace Acts {

namespace DetectorVolumeFinderJsonConverter {

/// @brief Convert the single delegate if it is of the type of the reference
///
/// @note It will do nothing if the type does not match
///
/// @param jIndexedVolumes the json object to be filled
/// @param delegate the delegate to be translated
/// @param detray indicates if detray json format is to be written out
/// @param refInstance is a reference instance of potential type casting
template <typename instance_type>
void convert(nlohmann::json& jIndexedVolumes,
             const Experimental::DetectorVolumeUpdator& delegate, bool detray,
             [[maybe_unused]] const instance_type& refInstance) {
  using GridType = typename instance_type::template grid_type<std::size_t>;
  // Defining a Delegate type
  using DelegateType = Experimental::IndexedUpdatorImpl<
      GridType, Acts::Experimental::IndexedDetectorVolumeExtractor,
      Acts::Experimental::DetectorVolumeFiller>;
  // Get the instance
  const auto* instance = delegate.instance();
  auto castedDelegate = dynamic_cast<const DelegateType*>(instance);
  if (castedDelegate != nullptr) {
    jIndexedVolumes = IndexedGridJsonHelper::convertImpl<DelegateType>(
        *castedDelegate, detray);
    if (detray) {
      jIndexedVolumes["volume_link"] = 1;
      nlohmann::json jAccLink;
      jAccLink["type"] =
          DetrayJsonHelper::accelerationLink(castedDelegate->casts);
      jAccLink["index"] = 0;
      jIndexedVolumes["acc_link"] = jAccLink;
    } else {
      jIndexedVolumes["type"] = "IndexedVolumes";
    }
  }
}

/// @brief Unrolling function for catching the right instance
///
/// @param jIndexedVolumes the json object to be filled
/// @param delegate the delegate to be translated
/// @param detray indicate if this is a detray json to be written out
/// @param axesTuple the tuple of axes to be unrolled
///
/// @note parameters are as of the `convertImpl` method
template <typename tuple_type, std::size_t... I>
void unrollConvert(nlohmann::json& jIndexedVolumes,
                   const Experimental::DetectorVolumeUpdator& delegate,
                   bool detray, const tuple_type& axesTuple,
                   std::index_sequence<I...> /*unused*/) {
  (convert(jIndexedVolumes, delegate, detray, std::get<I>(axesTuple)), ...);
}

/// Convert a volume finder
///
/// @param delegate the delegate to be translated
/// @param detray indicate if this is a detray json to be written out
///
/// @note this is the entry point of the conversion, i.e. top of the
/// unrolling loop
///
/// @note detray json format can not be read back in with Acts
///
/// @return a json object
static inline nlohmann::json toJson(
    const Experimental::DetectorVolumeUpdator& delegate, bool detray = false) {
  // Convert if dynamic cast happens to work
  nlohmann::json jIndexedVolumes;
  unrollConvert(jIndexedVolumes, delegate, detray,
                IndexedGridJsonHelper::s_possibleAxes,
                std::make_index_sequence<std::tuple_size<
                    decltype(IndexedGridJsonHelper::s_possibleAxes)>::value>());
  // Return the newly filled ones
  return jIndexedVolumes;
}

/// @brief Convert to a delegate from json
///
/// @param jVolumeFinder the json file to read from
///
/// @return the connected navigation delegate
Experimental::DetectorVolumeUpdator fromJson(
    const nlohmann::json& jVolumeFinder);

}  // namespace DetectorVolumeFinderJsonConverter
}  // namespace Acts
