// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiComponentTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/StepperOptions.hpp"
#include "Acts/Propagator/StepperStatistics.hpp"
#include "Acts/Propagator/detail/LoopStepperUtils.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <vector>

#include <boost/container/small_vector.hpp>

namespace Acts {

using namespace Acts::UnitLiterals;

namespace detail {

struct MaxMomentumComponent {
  template <typename component_range_t>
  auto operator()(const component_range_t& cmps) const {
    return std::max_element(cmps.begin(), cmps.end(),
                            [&](const auto& a, const auto& b) {
                              return std::abs(a.state.pars[eFreeQOverP]) >
                                     std::abs(b.state.pars[eFreeQOverP]);
                            });
  }
};

struct MaxWeightComponent {
  template <typename component_range_t>
  auto operator()(const component_range_t& cmps) {
    return std::max_element(
        cmps.begin(), cmps.end(),
        [&](const auto& a, const auto& b) { return a.weight < b.weight; });
  }
};

template <typename component_chooser_t>
struct SingleComponentReducer {
  template <typename stepper_state_t>
  static Vector3 position(const stepper_state_t& s) {
    return component_chooser_t{}(s.components)
        ->state.pars.template segment<3>(eFreePos0);
  }

  template <typename stepper_state_t>
  static Vector3 direction(const stepper_state_t& s) {
    return component_chooser_t{}(s.components)
        ->state.pars.template segment<3>(eFreeDir0);
  }

  template <typename stepper_state_t>
  static double qOverP(const stepper_state_t& s) {
    const auto cmp = component_chooser_t{}(s.components);
    return cmp->state.pars[eFreeQOverP];
  }

  template <typename stepper_state_t>
  static double absoluteMomentum(const stepper_state_t& s) {
    const auto cmp = component_chooser_t{}(s.components);
    return s.particleHypothesis.extractMomentum(cmp->state.pars[eFreeQOverP]);
  }

  template <typename stepper_state_t>
  static Vector3 momentum(const stepper_state_t& s) {
    const auto cmp = component_chooser_t{}(s.components);
    return s.particleHypothesis.extractMomentum(cmp->state.pars[eFreeQOverP]) *
           cmp->state.pars.template segment<3>(eFreeDir0);
  }

  template <typename stepper_state_t>
  static double charge(const stepper_state_t& s) {
    const auto cmp = component_chooser_t{}(s.components);
    return s.particleHypothesis.extractCharge(cmp->state.pars[eFreeQOverP]);
  }

  template <typename stepper_state_t>
  static double time(const stepper_state_t& s) {
    return component_chooser_t{}(s.components)->state.pars[eFreeTime];
  }

  template <typename stepper_state_t>
  static FreeVector pars(const stepper_state_t& s) {
    return component_chooser_t{}(s.components)->state.pars;
  }

  template <typename stepper_state_t>
  static FreeVector cov(const stepper_state_t& s) {
    return component_chooser_t{}(s.components)->state.cov;
  }
};

}  // namespace detail

using MaxMomentumReducerLoop =
    detail::SingleComponentReducer<detail::MaxMomentumComponent>;
using MaxWeightReducerLoop =
    detail::SingleComponentReducer<detail::MaxWeightComponent>;

/// @brief Stepper based on the EigenStepper, but handles Multi-Component Tracks
/// (e.g., for the GSF). Internally, this only manages a vector of
/// EigenStepper::States. This simplifies implementation, but has several
/// drawbacks:
/// * There are certain redundancies between the global State and the
/// component states
/// * The components do not share a single magnetic-field-cache
/// @tparam extension_t See EigenStepper for details
/// @tparam component_reducer_t How to map the multi-component state to a single
/// component
/// @tparam small_vector_size A size-hint how much memory should be allocated
/// by the small vector
template <typename extension_t = EigenStepperDefaultExtension,
          typename component_reducer_t = MaxWeightReducerLoop>
class MultiEigenStepperLoop : public EigenStepper<extension_t> {
  /// Limits the number of steps after at least one component reached the
  /// surface
  std::size_t m_stepLimitAfterFirstComponentOnSurface = 50;

  /// The logger (used if no logger is provided by caller of methods)
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Small vector type for speeding up some computations where we need to
  /// accumulate stuff of components. We think 16 is a reasonable amount here.
  template <typename T>
  using SmallVector = boost::container::small_vector<T, 16>;

 public:
  /// @brief Typedef to the Single-Component Eigen Stepper
  using SingleStepper = EigenStepper<extension_t>;

  /// @brief Typedef to the Single-Component Stepper Options
  using SingleOptions = typename SingleStepper::Options;

  /// @brief Typedef to the State of the single component Stepper
  using SingleState = typename SingleStepper::State;

  /// @brief Use the definitions from the Single-stepper
  using typename SingleStepper::Covariance;
  using typename SingleStepper::Jacobian;

  /// @brief Define an own bound state
  using BoundState =
      std::tuple<MultiComponentBoundTrackParameters, Jacobian, double>;

  /// @brief Define an own curvilinear state
  using CurvilinearState =
      std::tuple<MultiComponentCurvilinearTrackParameters, Jacobian, double>;

  /// @brief The reducer type
  using Reducer = component_reducer_t;

  /// @brief How many components can this stepper manage?
  static constexpr int maxComponents = std::numeric_limits<int>::max();

  struct Config {
    std::shared_ptr<const MagneticFieldProvider> bField;
  };

  struct Options : public SingleOptions {
    Options(const GeometryContext& gctx, const MagneticFieldContext& mctx)
        : SingleOptions(gctx, mctx) {}

    void setPlainOptions(const StepperPlainOptions& options) {
      static_cast<StepperPlainOptions&>(*this) = options;
    }
  };

  struct State {
    /// The struct that stores the individual components
    struct Component {
      SingleState state;
      double weight;
      IntersectionStatus status;
    };

    Options options;

    /// Particle hypothesis
    ParticleHypothesis particleHypothesis = ParticleHypothesis::pion();

    /// The components of which the state consists
    SmallVector<Component> components;

    bool covTransport = false;
    double pathAccumulated = 0.;
    std::size_t steps = 0;

    /// Step-limit counter which limits the number of steps when one component
    /// reached a surface
    std::optional<std::size_t> stepCounterAfterFirstComponentOnSurface;

    /// The stepper statistics
    StepperStatistics statistics;

    /// Constructor from the initial bound track parameters
    ///
    /// @param [in] optionsIn is the options object for the stepper
    ///
    /// @note the covariance matrix is copied when needed
    explicit State(const Options& optionsIn) : options(optionsIn) {}
  };

  /// Constructor from a magnetic field and a optionally provided Logger
  MultiEigenStepperLoop(std::shared_ptr<const MagneticFieldProvider> bField,
                        std::unique_ptr<const Logger> logger =
                            getDefaultLogger("GSF", Logging::INFO))
      : EigenStepper<extension_t>(std::move(bField)),
        m_logger(std::move(logger)) {}

  /// Constructor from a configuration and optionally provided Logger
  MultiEigenStepperLoop(const Config& config,
                        std::unique_ptr<const Logger> logger =
                            getDefaultLogger("GSF", Logging::INFO))
      : EigenStepper<extension_t>(config), m_logger(std::move(logger)) {}

  /// Construct and initialize a state
  State makeState(const Options& options,
                  const MultiComponentBoundTrackParameters& par) const {
    if (par.components().empty()) {
      throw std::invalid_argument(
          "Cannot construct MultiEigenStepperLoop::State with empty "
          "multi-component parameters");
    }

    State state(options);

    state.particleHypothesis = par.particleHypothesis();

    const auto surface = par.referenceSurface().getSharedPtr();

    for (auto i = 0ul; i < par.components().size(); ++i) {
      const auto& [weight, singlePars] = par[i];
      state.components.push_back({SingleStepper::makeState(options, singlePars),
                                  weight, IntersectionStatus::onSurface});
    }

    if (std::get<2>(par.components().front())) {
      state.covTransport = true;
    }

    return state;
  }

  /// @brief Resets the state
  ///
  /// @param [in, out] state State of the stepper
  /// @param [in] boundParams Parameters in bound parametrisation
  /// @param [in] cov Covariance matrix
  /// @param [in] surface The reference surface of the bound parameters
  /// @param [in] stepSize Step size
  void resetState(
      State& state, const BoundVector& boundParams,
      const BoundSquareMatrix& cov, const Surface& surface,
      const double stepSize = std::numeric_limits<double>::max()) const {
    for (auto& component : state.components) {
      SingleStepper::resetState(component.state, boundParams, cov, surface,
                                stepSize);
    }
  }

  /// A proxy struct which allows access to a single component of the
  /// multi-component state. It has the semantics of a const reference, i.e.
  /// it requires a const reference of the single-component state it
  /// represents
  using ConstComponentProxy =
      detail::LoopComponentProxyBase<const typename State::Component,
                                     MultiEigenStepperLoop>;

  /// A proxy struct which allows access to a single component of the
  /// multi-component state. It has the semantics of a mutable reference, i.e.
  /// it requires a mutable reference of the single-component state it
  /// represents
  using ComponentProxy = detail::LoopComponentProxy<typename State::Component,
                                                    MultiEigenStepperLoop>;

  /// Creates an iterable which can be plugged into a range-based for-loop to
  /// iterate over components
  /// @note Use a for-loop with by-value semantics, since the Iterable returns a
  /// proxy internally holding a reference
  auto componentIterable(State& state) const {
    struct Iterator {
      using difference_type [[maybe_unused]] = std::ptrdiff_t;
      using value_type [[maybe_unused]] = ComponentProxy;
      using reference [[maybe_unused]] = ComponentProxy;
      using pointer [[maybe_unused]] = void;
      using iterator_category [[maybe_unused]] = std::forward_iterator_tag;

      typename decltype(state.components)::iterator it;
      const State& s;

      // clang-format off
      auto& operator++() { ++it; return *this; }
      auto operator==(const Iterator& other) const { return it == other.it; }
      auto operator*() const { return ComponentProxy(*it, s); }
      // clang-format on
    };

    struct Iterable {
      State& s;

      // clang-format off
      auto begin() { return Iterator{s.components.begin(), s}; }
      auto end() { return Iterator{s.components.end(), s}; }
      // clang-format on
    };

    return Iterable{state};
  }

  /// Creates an constant iterable which can be plugged into a range-based
  /// for-loop to iterate over components
  /// @note Use a for-loop with by-value semantics, since the Iterable returns a
  /// proxy internally holding a reference
  auto constComponentIterable(const State& state) const {
    struct ConstIterator {
      using difference_type [[maybe_unused]] = std::ptrdiff_t;
      using value_type [[maybe_unused]] = ConstComponentProxy;
      using reference [[maybe_unused]] = ConstComponentProxy;
      using pointer [[maybe_unused]] = void;
      using iterator_category [[maybe_unused]] = std::forward_iterator_tag;

      typename decltype(state.components)::const_iterator it;
      const State& s;

      // clang-format off
      auto& operator++() { ++it; return *this; }
      auto operator==(const ConstIterator& other) const { return it == other.it; }
      auto operator*() const { return ConstComponentProxy{*it}; }
      // clang-format on
    };

    struct Iterable {
      const State& s;

      // clang-format off
      auto begin() const { return ConstIterator{s.components.cbegin(), s}; }
      auto end() const { return ConstIterator{s.components.cend(), s}; }
      // clang-format on
    };

    return Iterable{state};
  }

  /// Get the number of components
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  std::size_t numberComponents(const State& state) const {
    return state.components.size();
  }

  /// Remove missed components from the component state
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  void removeMissedComponents(State& state) const {
    auto new_end = std::remove_if(
        state.components.begin(), state.components.end(), [](const auto& cmp) {
          return cmp.status == IntersectionStatus::unreachable;
        });

    state.components.erase(new_end, state.components.end());
  }

  /// Reweight the components
  ///
  /// @param [in,out] state The stepping state (thread-local cache)
  void reweightComponents(State& state) const {
    double sumOfWeights = 0.0;
    for (const auto& cmp : state.components) {
      sumOfWeights += cmp.weight;
    }
    for (auto& cmp : state.components) {
      cmp.weight /= sumOfWeights;
    }
  }

  /// Reset the number of components
  ///
  /// @param [in,out] state  The stepping state (thread-local cache)
  void clearComponents(State& state) const { state.components.clear(); }

  /// Add a component to the Multistepper
  ///
  /// @param [in,out] state  The stepping state (thread-local cache)
  /// @param [in] pars Parameters of the component to add
  /// @param [in] weight Weight of the component to add
  ///
  /// @note: It is not ensured that the weights are normalized afterwards
  /// @note This function makes no garantuees about how new components are
  /// initialized, it is up to the caller to ensure that all components are
  /// valid in the end.
  /// @note The returned component-proxy is only garantueed to be valid until
  /// the component number is again modified
  Result<ComponentProxy> addComponent(State& state,
                                      const BoundTrackParameters& pars,
                                      double weight) const {
    state.components.push_back({SingleStepper::makeState(state.options, pars),
                                weight, IntersectionStatus::onSurface});

    return ComponentProxy{state.components.back(), state};
  }

  /// Get the field for the stepping, it checks first if the access is still
  /// within the Cell, and updates the cell if necessary.
  ///
  /// @param [in,out] state is the propagation state associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  ///
  /// @note This uses the cache of the first component stored in the state
  Result<Vector3> getField(State& state, const Vector3& pos) const {
    return SingleStepper::getField(state.components.front().state, pos);
  }

  /// Global particle position accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  Vector3 position(const State& state) const {
    return Reducer::position(state);
  }

  /// Momentum direction accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  Vector3 direction(const State& state) const {
    return Reducer::direction(state);
  }

  /// QoP access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double qOverP(const State& state) const { return Reducer::qOverP(state); }

  /// Absolute momentum accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double absoluteMomentum(const State& state) const {
    return Reducer::absoluteMomentum(state);
  }

  /// Momentum accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  Vector3 momentum(const State& state) const {
    return Reducer::momentum(state);
  }

  /// Charge access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double charge(const State& state) const { return Reducer::charge(state); }

  /// Particle hypothesis
  ///
  /// @param state [in] The stepping state (thread-local cache)
  ParticleHypothesis particleHypothesis(const State& state) const {
    return state.particleHypothesis;
  }

  /// Time access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double time(const State& state) const { return Reducer::time(state); }

  /// Update surface status
  ///
  /// It checks the status to the reference surface & updates
  /// the step size accordingly
  ///
  /// @param [in,out] state The stepping state (thread-local cache)
  /// @param [in] surface The surface provided
  /// @param [in] index The surface intersection index
  /// @param [in] navDir The navigation direction
  /// @param [in] boundaryTolerance The boundary check for this status update
  /// @param [in] surfaceTolerance Surface tolerance used for intersection
  /// @param [in] stype The step size type to be set
  /// @param [in] logger A @c Logger instance
  IntersectionStatus updateSurfaceStatus(
      State& state, const Surface& surface, std::uint8_t index,
      Direction navDir, const BoundaryTolerance& boundaryTolerance,
      double surfaceTolerance, ConstrainedStep::Type stype,
      const Logger& logger = getDummyLogger()) const {
    using Status = IntersectionStatus;

    std::array<int, 3> counts = {0, 0, 0};

    for (auto& component : state.components) {
      component.status = detail::updateSingleSurfaceStatus<SingleStepper>(
          *this, component.state, surface, index, navDir, boundaryTolerance,
          surfaceTolerance, stype, logger);
      ++counts[static_cast<std::size_t>(component.status)];
    }

    // If at least one component is on a surface, we can remove all missed
    // components before the step. If not, we must keep them for the case that
    // all components miss and we need to retarget
    if (counts[static_cast<std::size_t>(Status::onSurface)] > 0) {
      removeMissedComponents(state);
      reweightComponents(state);
    }

    ACTS_VERBOSE("Component status wrt "
                 << surface.geometryId() << " at {"
                 << surface.center(state.options.geoContext).transpose()
                 << "}:\t" << [&]() {
                      std::stringstream ss;
                      for (auto& component : state.components) {
                        ss << component.status << "\t";
                      }
                      return ss.str();
                    }());

    // Switch on stepCounter if one or more components reached a surface, but
    // some are still in progress of reaching the surface
    if (!state.stepCounterAfterFirstComponentOnSurface &&
        counts[static_cast<std::size_t>(Status::onSurface)] > 0 &&
        counts[static_cast<std::size_t>(Status::reachable)] > 0) {
      state.stepCounterAfterFirstComponentOnSurface = 0;
      ACTS_VERBOSE("started stepCounterAfterFirstComponentOnSurface");
    }

    // If there are no components onSurface, but the counter is switched on
    // (e.g., if the navigator changes the target surface), we need to switch it
    // off again
    if (state.stepCounterAfterFirstComponentOnSurface &&
        counts[static_cast<std::size_t>(Status::onSurface)] == 0) {
      state.stepCounterAfterFirstComponentOnSurface.reset();
      ACTS_VERBOSE("switch off stepCounterAfterFirstComponentOnSurface");
    }

    // This is a 'any_of' criterium. As long as any of the components has a
    // certain state, this determines the total state (in the order of a
    // somewhat importance)
    if (counts[static_cast<std::size_t>(Status::reachable)] > 0) {
      return Status::reachable;
    } else if (counts[static_cast<std::size_t>(Status::onSurface)] > 0) {
      state.stepCounterAfterFirstComponentOnSurface.reset();
      return Status::onSurface;
    } else {
      return Status::unreachable;
    }
  }

  /// Update step size
  ///
  /// This method intersects the provided surface and update the navigation
  /// step estimation accordingly (hence it changes the state). It also
  /// returns the status of the intersection to trigger onSurface in case
  /// the surface is reached.
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param oIntersection [in] The ObjectIntersection to layer, boundary, etc
  /// @param direction [in] The propagation direction
  /// @param stype [in] The step size type to be set
  template <typename object_intersection_t>
  void updateStepSize(State& state, const object_intersection_t& oIntersection,
                      Direction direction, ConstrainedStep::Type stype) const {
    const Surface& surface = *oIntersection.object();

    for (auto& component : state.components) {
      auto intersection = surface.intersect(
          component.state.options.geoContext,
          SingleStepper::position(component.state),
          direction * SingleStepper::direction(component.state),
          BoundaryTolerance::None())[oIntersection.index()];

      SingleStepper::updateStepSize(component.state, intersection, direction,
                                    stype);
    }
  }

  /// Update step size - explicitly with a double
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param stepSize [in] The step size value
  /// @param stype [in] The step size type to be set
  void updateStepSize(State& state, double stepSize,
                      ConstrainedStep::Type stype) const {
    for (auto& component : state.components) {
      SingleStepper::updateStepSize(component.state, stepSize, stype);
    }
  }

  /// Get the step size
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @param stype [in] The step size type to be returned
  /// @note This returns the smallest step size of all components. It uses
  /// std::abs for comparison to handle backward propagation and negative
  /// step sizes correctly.
  double getStepSize(const State& state, ConstrainedStep::Type stype) const {
    return std::min_element(state.components.begin(), state.components.end(),
                            [=](const auto& a, const auto& b) {
                              return std::abs(a.state.stepSize.value(stype)) <
                                     std::abs(b.state.stepSize.value(stype));
                            })
        ->state.stepSize.value(stype);
  }

  /// Release the step-size for all components
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param [in] stype The step size type to be released
  void releaseStepSize(State& state, ConstrainedStep::Type stype) const {
    for (auto& component : state.components) {
      SingleStepper::releaseStepSize(component.state, stype);
    }
  }

  /// Output the Step Size of all components into one std::string
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  std::string outputStepSize(const State& state) const {
    std::stringstream ss;
    for (const auto& component : state.components) {
      ss << component.state.stepSize.toString() << " || ";
    }

    return ss.str();
  }

  /// Create and return the bound state at the current position
  ///
  /// @brief This transports (if necessary) the covariance
  /// to the surface and creates a bound state. It does not check
  /// if the transported state is at the surface, this needs to
  /// be guaranteed by the propagator.
  /// @note This is done by combining the gaussian mixture on the specified
  /// surface. If the conversion to bound states of some components
  /// fails, these components are ignored unless all components fail. In this
  /// case an error code is returned.
  ///
  /// @param [in] state State that will be presented as @c BoundState
  /// @param [in] surface The surface to which we bind the state
  /// @param [in] transportCov Flag steering covariance transport
  /// @param [in] freeToBoundCorrection Flag steering non-linear correction during global to local correction
  ///
  /// @return A bound state:
  ///   - the parameters at the surface
  ///   - the stepwise jacobian towards it (from last bound)
  ///   - and the path length (from start - for ordering)
  Result<BoundState> boundState(
      State& state, const Surface& surface, bool transportCov = true,
      const FreeToBoundCorrection& freeToBoundCorrection =
          FreeToBoundCorrection(false)) const;

  /// @brief If necessary fill additional members needed for curvilinearState
  ///
  /// Compute path length derivatives in case they have not been computed
  /// yet, which is the case if no step has been executed yet.
  ///
  /// @param [in, out] state The stepping state (thread-local cache)
  /// @return true if nothing is missing after this call, false otherwise.
  bool prepareCurvilinearState(State& state) const {
    (void)state;
    return true;
  }

  /// Create and return a curvilinear state at the current position
  ///
  /// @brief This transports (if necessary) the covariance
  /// to the current position and creates a curvilinear state.
  /// @note This is done as a simple average over the free representation
  /// and covariance of the components.
  ///
  /// @param [in] state State that will be presented as @c CurvilinearState
  /// @param [in] transportCov Flag steering covariance transport
  ///
  /// @return A curvilinear state:
  ///   - the curvilinear parameters at given position
  ///   - the stepweise jacobian towards it (from last bound)
  ///   - and the path length (from start - for ordering)
  CurvilinearState curvilinearState(State& state,
                                    bool transportCov = true) const;

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state
  ///
  /// @param [in,out] state State of the stepper
  void transportCovarianceToCurvilinear(State& state) const {
    for (auto& component : state.components) {
      SingleStepper::transportCovarianceToCurvilinear(component.state);
    }
  }

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current position,
  /// or direction of the state
  ///
  /// @tparam surface_t the Surface type
  ///
  /// @param [in,out] state State of the stepper
  /// @param [in] surface is the surface to which the covariance is forwarded
  /// @param [in] freeToBoundCorrection Flag steering non-linear correction during global to local correction
  /// to
  /// @note no check is done if the position is actually on the surface
  void transportCovarianceToBound(
      State& state, const Surface& surface,
      const FreeToBoundCorrection& freeToBoundCorrection =
          FreeToBoundCorrection(false)) const {
    for (auto& component : state.components) {
      SingleStepper::transportCovarianceToBound(component.state, surface,
                                                freeToBoundCorrection);
    }
  }

  /// Perform a Runge-Kutta track parameter propagation step
  ///
  /// @param [in,out] state The state of the stepper
  /// @param propDir is the direction of propagation
  /// @param material is the material properties
  /// @return the result of the step
  ///
  /// The state contains the desired step size. It can be negative during
  /// backwards track propagation, and since we're using an adaptive
  /// algorithm, it can be modified by the stepper class during propagation.
  Result<double> step(State& state, Direction propDir,
                      const IVolumeMaterial* material) const;
};

}  // namespace Acts

#include "MultiEigenStepperLoop.ipp"
