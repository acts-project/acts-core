// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/MultiEigenStepperLoop.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/GaussianSumFitter.hpp"
#include "Acts/TrackFitting/GsfMixtureReduction.hpp"
#include "Acts/TrackFitting/GsfOptions.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/TrackFitting/RefittingCalibrator.hpp"
#include "ActsExamples/TrackFitting/TrackFitterFunction.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace Acts {
class MagneticFieldProvider;
class SourceLink;
class Surface;
class TrackingGeometry;
}  // namespace Acts

using namespace ActsExamples;

namespace {

using MultiStepper =
    Acts::MultiEigenStepperLoop<Acts::EigenStepperDefaultExtension,
                                Acts::MaxWeightReducerLoop>;
using Propagator = Acts::Propagator<MultiStepper, Acts::Navigator>;
using DirectPropagator = Acts::Propagator<MultiStepper, Acts::DirectNavigator>;

using Fitter = Acts::GaussianSumFitter<Propagator, BetheHeitlerApprox,
                                       Acts::VectorMultiTrajectory>;
using DirectFitter =
    Acts::GaussianSumFitter<DirectPropagator, BetheHeitlerApprox,
                            Acts::VectorMultiTrajectory>;
using TrackContainer =
    Acts::TrackContainer<Acts::VectorTrackContainer,
                         Acts::VectorMultiTrajectory, std::shared_ptr>;

struct GsfFitterFunctionImpl final : public ActsExamples::TrackFitterFunction {
  Fitter fitter;
  DirectFitter directFitter;

  Acts::GainMatrixUpdater updater;

  std::size_t maxComponents = 0;
  double weightCutoff = 0;
  const double momentumCutoff = 0;  // 500_MeV;
  bool abortOnError = false;
  bool disableAllMaterialHandling = false;
  MixtureReductionAlgorithm reductionAlg =
      MixtureReductionAlgorithm::KLDistance;
  Acts::ComponentMergeMethod mergeMethod =
      Acts::ComponentMergeMethod::eMaxWeight;
  double reverseFilteringCovarianceScaling = 1.0;

  IndexSourceLink::SurfaceAccessor m_slSurfaceAccessor;

  GsfFitterFunctionImpl(Fitter&& f, DirectFitter&& df,
                        const Acts::TrackingGeometry& trkGeo)
      : fitter(std::move(f)),
        directFitter(std::move(df)),
        m_slSurfaceAccessor{trkGeo} {}

  template <typename calibrator_t>
  auto makeGsfOptions(const GeneralFitterOptions& options,
                      const calibrator_t& calibrator) const {
    Acts::GsfExtensions<Acts::VectorMultiTrajectory> extensions;
    extensions.updater.connect<
        &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
        &updater);

    Acts::GsfOptions<Acts::VectorMultiTrajectory> gsfOptions{
        options.geoContext, options.magFieldContext,
        options.calibrationContext};
    gsfOptions.extensions = extensions;
    gsfOptions.propagatorPlainOptions = options.propOptions;
    gsfOptions.referenceSurface = options.referenceSurface;
    gsfOptions.maxComponents = maxComponents;
    gsfOptions.weightCutoff = weightCutoff;
    gsfOptions.abortOnError = abortOnError;
    gsfOptions.disableAllMaterialHandling = disableAllMaterialHandling;
    gsfOptions.componentMergeMethod = mergeMethod;
    gsfOptions.reverseFilteringCovarianceScaling =
        reverseFilteringCovarianceScaling;

    gsfOptions.extensions.calibrator.connect<&calibrator_t::calibrate>(
        &calibrator);

    if (options.doRefit) {
      gsfOptions.extensions.surfaceAccessor
          .connect<&RefittingCalibrator::accessSurface>();
    } else {
      gsfOptions.extensions.surfaceAccessor
          .connect<&IndexSourceLink::SurfaceAccessor::operator()>(
              &m_slSurfaceAccessor);
    }
    switch (reductionAlg) {
      case MixtureReductionAlgorithm::weightCut: {
        gsfOptions.extensions.mixtureReducer
            .connect<&Acts::reduceMixtureLargestWeights>();
      } break;
      case MixtureReductionAlgorithm::KLDistance: {
        gsfOptions.extensions.mixtureReducer
            .connect<&Acts::reduceMixtureWithKLDistance>();
      } break;
    }

    return gsfOptions;
  }

  TrackFitterResult operator()(const std::vector<Acts::SourceLink>& sourceLinks,
                               const TrackParameters& initialParameters,
                               const GeneralFitterOptions& options,
                               const MeasurementCalibratorAdapter& calibrator,
                               TrackContainer& tracks) const override {
    const auto gsfOptions = makeGsfOptions(options, calibrator);

    using namespace Acts::GsfConstants;
    if (!tracks.hasColumn(Acts::hashString(kFinalMultiComponentStateColumn))) {
      std::string key(kFinalMultiComponentStateColumn);
      tracks.template addColumn<FinalMultiComponentState>(key);
    }

    if (!tracks.hasColumn(Acts::hashString(kFwdMaxMaterialXOverX0))) {
      tracks.template addColumn<double>(std::string(kFwdMaxMaterialXOverX0));
    }

    if (!tracks.hasColumn(Acts::hashString(kFwdSumMaterialXOverX0))) {
      tracks.template addColumn<double>(std::string(kFwdSumMaterialXOverX0));
    }

    return fitter.fit(sourceLinks.begin(), sourceLinks.end(), initialParameters,
                      gsfOptions, tracks);
  }

  TrackFitterResult operator()(
      const std::vector<Acts::SourceLink>& sourceLinks,
      const TrackParameters& initialParameters,
      const GeneralFitterOptions& options,
      const RefittingCalibrator& calibrator,
      const std::vector<const Acts::Surface*>& surfaceSequence,
      TrackContainer& tracks) const override {
    const auto gsfOptions = makeGsfOptions(options, calibrator);

    using namespace Acts::GsfConstants;
    if (!tracks.hasColumn(Acts::hashString(kFinalMultiComponentStateColumn))) {
      std::string key(kFinalMultiComponentStateColumn);
      tracks.template addColumn<FinalMultiComponentState>(key);
    }

    return directFitter.fit(sourceLinks.begin(), sourceLinks.end(),
                            initialParameters, gsfOptions, surfaceSequence,
                            tracks);
  }
};

}  // namespace

std::shared_ptr<TrackFitterFunction> ActsExamples::makeGsfFitterFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    BetheHeitlerApprox betheHeitlerApprox, std::size_t maxComponents,
    double weightCutoff, Acts::ComponentMergeMethod componentMergeMethod,
    MixtureReductionAlgorithm mixtureReductionAlgorithm,
    double reverseFilteringCovarianceScaling, const Acts::Logger& logger) {
  // Standard fitter
  MultiStepper stepper(magneticField, logger.cloneWithSuffix("Step"));
  const auto& geo = *trackingGeometry;
  Acts::Navigator::Config cfg{std::move(trackingGeometry)};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg, logger.cloneWithSuffix("Navigator"));
  Propagator propagator(std::move(stepper), std::move(navigator),
                        logger.cloneWithSuffix("Propagator"));
  Fitter trackFitter(std::move(propagator),
                     BetheHeitlerApprox(betheHeitlerApprox),
                     logger.cloneWithSuffix("GSF"));

  // Direct fitter
  MultiStepper directStepper(std::move(magneticField),
                             logger.cloneWithSuffix("Step"));
  Acts::DirectNavigator directNavigator{
      logger.cloneWithSuffix("DirectNavigator")};
  DirectPropagator directPropagator(std::move(directStepper),
                                    std::move(directNavigator),
                                    logger.cloneWithSuffix("DirectPropagator"));
  DirectFitter directTrackFitter(std::move(directPropagator),
                                 BetheHeitlerApprox(betheHeitlerApprox),
                                 logger.cloneWithSuffix("DirectGSF"));

  // build the fitter functions. owns the fitter object.
  auto fitterFunction = std::make_shared<GsfFitterFunctionImpl>(
      std::move(trackFitter), std::move(directTrackFitter), geo);
  fitterFunction->maxComponents = maxComponents;
  fitterFunction->weightCutoff = weightCutoff;
  fitterFunction->mergeMethod = componentMergeMethod;
  fitterFunction->reductionAlg = mixtureReductionAlgorithm;
  fitterFunction->reverseFilteringCovarianceScaling =
      reverseFilteringCovarianceScaling;

  return fitterFunction;
}
