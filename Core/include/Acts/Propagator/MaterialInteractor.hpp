// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/Propagator/detail/VolumeMaterialInteraction.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <sstream>

namespace Acts {

/// Material interactor propagator action.
///
/// Apply material interactions at a surface and update the track state.
struct MaterialInteractor {
  /// Whether to consider multiple scattering.
  bool multipleScattering = true;
  /// Whether to consider energy loss.
  bool energyLoss = true;
  /// Whether to record all material interactions.
  bool recordInteractions = false;

  using result_type = RecordedMaterial;

  /// @brief Interaction with detector material for the ActionList
  /// of the Propagator
  ///
  /// It checks if the state has a current surface, in which case
  /// the action is performed: the covariance is transported to the position,
  /// multiple scattering and energy loss is applied  according to the
  /// configuration.
  ///
  /// @tparam propagator_state_t is the type of Propagagor state
  /// @tparam stepper_t Type of the stepper of the propagation
  ///
  /// @param state is the mutable propagator state object
  /// @param stepper The stepper in use
  /// @param result is the mutable result state object
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    const auto& logger = state.options.logger;
    // In case of Volume material update the result of the previous step
    if (recordInteractions && !result.materialInteractions.empty() &&
        result.materialInteractions.back().volume != nullptr &&
        result.materialInteractions.back().updatedVolumeStep == false) {
      UpdateResult(state, stepper, result);
    }

    // If we are on target, everything should have been done
    if (state.navigation.targetReached) {
      return;
    }
    // Do nothing if nothing is what is requested.
    if (not(multipleScattering or energyLoss or recordInteractions)) {
      return;
    }
    // We only have material interactions if there is potential material
    const Surface* surface = state.navigation.currentSurface;
    const TrackingVolume* volume = state.navigation.currentVolume;

    if (not(surface and surface->surfaceMaterial()) and
        not(volume and volume->volumeMaterial())) {
      return;
    }

    if (surface and surface->surfaceMaterial()) {
      // Prepare relevant input particle properties
      detail::PointwiseMaterialInteraction d(surface, state, stepper);

      // Determine the effective traversed material and its properties
      // Material exists but it's not real, i.e. vacuum; there is nothing to do
      if (not d.evaluateMaterialSlab(state)) {
        return;
      }

      // Evaluate the material effects
      d.evaluatePointwiseMaterialInteraction(multipleScattering, energyLoss);

      if (energyLoss) {
        using namespace UnitLiterals;
        ACTS_VERBOSE(d.slab << " pdg=" << d.pdg << " mass=" << d.mass / 1_MeV
                            << "MeV"
                            << " momentum=" << d.momentum / 1_GeV << "GeV"
                            << " energyloss=" << d.Eloss / 1_MeV << "MeV");
      }

      // To integrate process noise, we need to transport
      // the covariance to the current position in space
      if (d.performCovarianceTransport) {
        stepper.transportCovarianceToCurvilinear(state.stepping);
      }
      // Apply the material interactions
      d.updateState(state, stepper);
      // Record the result
      recordResult(d, result);
    } else if (recordInteractions && volume and volume->volumeMaterial()) {
      // Prepare relevant input particle properties
      detail::VolumeMaterialInteraction d(volume, state, stepper);
      // Determine the effective traversed material and its properties
      // Material exists but it's not real, i.e. vacuum; there is nothing to do
      if (not d.evaluateMaterialSlab(state)) {
        return;
      }
      // Record the result
      recordResult(d, result);
    }
  }

  /// Material interaction has no pure observer.
  template <typename propagator_state_t>
  void operator()(propagator_state_t& /* unused */) const {}

 private:
  /// @brief This function records the material effect
  ///
  /// @param [in] d Data cache container
  /// @param [in, out] result Result storage
  void recordResult(const detail::PointwiseMaterialInteraction& d,
                    result_type& result) const {
    result.materialInX0 += d.slab.thicknessInX0();
    result.materialInL0 += d.slab.thicknessInL0();
    // Record the interaction if requested
    if (recordInteractions) {
      MaterialInteraction mi;
      mi.position = d.pos;
      mi.time = d.time;
      mi.direction = d.dir;
      mi.deltaP = d.nextP - d.momentum;
      mi.sigmaPhi2 = d.variancePhi;
      mi.sigmaTheta2 = d.varianceTheta;
      mi.sigmaQoP2 = d.varianceQoverP;
      mi.surface = d.surface;
      mi.volume = nullptr;
      mi.pathCorrection = d.pathCorrection;
      mi.materialSlab = d.slab;
      result.materialInteractions.push_back(std::move(mi));
    }
  }

  /// @brief This function records the material effect
  ///
  /// @param [in] d Data cache container
  /// @param [in, out] result Result storage
  void recordResult(const detail::VolumeMaterialInteraction& d,
                    result_type& result) const {
    // Record the interaction
    MaterialInteraction mi;
    mi.position = d.pos;
    mi.time = d.time;
    mi.direction = d.dir;
    mi.surface = nullptr;
    mi.volume = d.volume;
    mi.pathCorrection = d.pathCorrection;
    mi.materialSlab = d.slab;
    result.materialInteractions.push_back(std::move(mi));
  }

  /// @brief This function update the previous material step
  ///
  /// @param [in,out] state The state object
  /// @param [in] stepper The stepper instance
  /// @param [in, out] result Result storage
  template <typename propagator_state_t, typename stepper_t>
  void UpdateResult(propagator_state_t& state, const stepper_t& stepper,
                    result_type& result) const {
    // Update the previous interaction
    Vector3 shift = stepper.position(state.stepping) -
                    result.materialInteractions.back().position;
    double momentum = stepper.direction(state.stepping).norm();
    result.materialInteractions.back().deltaP =
        momentum - result.materialInteractions.back().direction.norm();
    result.materialInteractions.back().materialSlab.scaleThickness(
        shift.norm());
    result.materialInteractions.back().updatedVolumeStep = true;
    result.materialInX0 +=
        result.materialInteractions.back().materialSlab.thicknessInX0();
    result.materialInL0 +=
        result.materialInteractions.back().materialSlab.thicknessInL0();
  }
};  // namespace Acts

}  // namespace Acts
