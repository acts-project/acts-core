#!/usr/bin/env python3

from pathlib import Path

import acts
import acts.examples

from truth_tracking_kalman import runTruthTrackingKalman

u = acts.UnitConstants


def runRefittingKf(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    digiConfigFile: Path,
    outputDir: Path,
    multipleScattering: bool = True,
    energyLoss: bool = True,
    reverseFilteringMomThreshold=0 * u.GeV,
    reverseFilteringCovarianceScaling=1.0,
    s: acts.examples.Sequencer = None,
):
    s = runTruthTrackingKalman(
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        outputDir=outputDir,
        reverseFilteringMomThreshold=reverseFilteringMomThreshold,
        reverseFilteringCovarianceScaling=reverseFilteringCovarianceScaling,
        s=s,
    )

    kalmanOptions = {
        "multipleScattering": multipleScattering,
        "energyLoss": energyLoss,
        "reverseFilteringMomThreshold": reverseFilteringMomThreshold,
        "reverseFilteringCovarianceScaling": reverseFilteringCovarianceScaling,
        "freeToBoundCorrection": acts.examples.FreeToBoundCorrection(False),
        "level": acts.logging.INFO,
    }

    s.addAlgorithm(
        acts.examples.RefittingAlgorithm(
            level=acts.logging.INFO,
            inputTracks="kf_tracks",
            outputTracks="kf_refit_tracks",
            initialVarInflation=6 * [100.0],
            fit=acts.examples.makeKalmanFitterFunction(
                trackingGeometry, field, **kalmanOptions
            ),
        )
    )

    s.addAlgorithm(
        acts.examples.TrackTruthMatcher(
            level=acts.logging.INFO,
            inputTracks="kf_refit_tracks",
            inputParticles="particles_selected",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputTrackParticleMatching="refit_track_particle_matching",
            outputParticleTrackMatching="refit_particle_track_matching",
        )
    )

    s.addWriter(
        acts.examples.RootTrackStatesWriter(
            level=acts.logging.INFO,
            inputTracks="kf_refit_tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="refit_track_particle_matching",
            inputSimHits="simhits",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(outputDir / "trackstates_kf_refit.root"),
        )
    )

    s.addWriter(
        acts.examples.RootTrackSummaryWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="refit_track_particle_matching",
            filePath=str(outputDir / "tracksummary_kf_refit.root"),
        )
    )

    s.addWriter(
        acts.examples.TrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "performance_kf_refit.root"),
        )
    )

    return s


if __name__ == "__main__":
    srcdir = Path(__file__).resolve().parent.parent.parent.parent
    outputDir = Path.cwd()

    # ODD
    from acts.examples.odd import getOpenDataDetector

    detector = getOpenDataDetector()
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()
    digiConfigFile = (
        srcdir / "thirdparty/OpenDataDetector/config/odd-digi-smearing-config.json"
    )

    ## GenericDetector
    # detector = acts.examples.GenericDetector()
    # trackingGeometry = detector.trackingGeometry()
    # digiConfigFile = (
    #     srcdir
    #     / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
    # )

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runRefittingKf(
        trackingGeometry=trackingGeometry,
        field=field,
        digiConfigFile=digiConfigFile,
        outputDir=Path.cwd(),
    ).run()
