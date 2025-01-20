// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialSlab.hpp"

#include <utility>

namespace Acts {

/// Accumulate material properties from multiple hits/track and multiple tracks.
///
/// This is a helper class for the `SurfaceMaterialMapper` to handle material
/// accumulation and averaging for one surface bin. The accumulation procedure
/// is done in two steps:
///
/// 1.  The per-track store accumulates material steps from one track/particle.
///     Multiple material steps can be assigned to the same bin time by
///     one particle, e.g. if the simulation has created more than one step in
///     the material or if several components are compressed into one
///     description. Multiple steps are treated as if they are passed one after
///     the other.
/// 2.  The total store averages the accumulated material properties over all
///     tracks. Each track contributes equally.
class AccumulatedMaterialSlab {
 public:
  // this class does not have a custom default constructor and thus should not
  // provide any custom default cstors, dstor, or assignment. see ISOCPP C.20.

  /// Add the material to the current per-track store.
  ///
  /// @param slabAlongTrack Recorded equivalent material slab for this step
  /// @param pathCorrection Correction factor due to non-perpendicular incident
  ///
  /// The recoded material slab is assumed to be defined along the track
  /// direction. The track can have non-perpendicular incidence on the surface
  /// and the recorded slab has to be projected along the surface normal. The
  /// path correction gives the scaling factor from normal incidence to the
  /// recorded incidence as provided by the `Surface` interface.
  ///
  ///  Vacuum steps with a non-zero thickness can be added to account for holes
  ///  in material structures.
  void accumulate(MaterialSlab slabAlongTrack, double pathCorrection = 1.);

  /// Use the accumulated material to update the material variance
  ///
  /// @param slabReference reference slab (from the map) used to compute the variance
  /// @param useEmptyTrack indicate whether to consider an empty track store
  ///
  /// The material variance can be used to optimised the mapping process as it
  /// should be inversely proportional to the map quality
  void trackVariance(MaterialSlab slabReference, bool useEmptyTrack = false);

  /// Add the accumulated material for the current track to the total average.
  ///
  /// @param useEmptyTrack indicate whether to consider an empty track store
  ///
  /// This finishes the material accumulation for the current track and resets
  /// the per-track store. Subsequent calls to `.accumulate(...)` will start
  /// accumulating material for a new track.
  ///
  /// Each track contributes equally to the total average regardless of its
  /// measured path within the material. An empty per-track store, i.e.
  /// vanishing per-track material thickness, does not contribute to the total
  /// unless explicitly requested.
  void trackAverage(bool useEmptyTrack = false);

  /// Return the average material properties from all accumulated tracks.
  ///
  /// @returns Average material properties and the number of contributing tracks
  ///
  /// Only contains the information up to the last `.trackAverage(...)` call. If
  /// there have been additional calls to `.accumulate(...)` afterwards, the
  /// information is not part of the total average. The thickness corresponds to
  /// the average thickness seen by the tracks.
  std::pair<MaterialSlab, unsigned int> totalAverage() const;

  /// Return the material variance from all accumulated tracks.
  ///
  /// @returns Average material properties and the number of contributing tracks
  ///
  /// Only contains the information up to the last `.trackVariance(...)` call.
  /// If there have been additional calls to `.accumulate(...)` afterwards, the
  /// information is not part of the total average. The number of tracks is only
  /// updated on the call of `.trackAverage(...)`
  std::pair<float, unsigned int> totalVariance() const;

 private:
  /// Averaged properties for a single track.
  MaterialSlab m_trackAverage;
  /// Averaged properties over multiple tracks.
  MaterialSlab m_totalAverage;
  /// Averaged variance over multiple tracks.
  float m_totalVariance = 0.0;
  // Number of tracks contributing to the total average.
  unsigned int m_totalCount = 0u;
};

}  // namespace Acts
