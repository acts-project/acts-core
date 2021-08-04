#!/usr/bin/env python3
from Examples.Scripts.Python.common import getOpenDataDetector
import os
from pathlib import Path

import acts
import acts.examples
import acts.examples.dd4hep

u = acts.UnitConstants


def runTruthTracking(
    trackingGeometry,
    field,
    outputDir: Path,
    digiConfigFile: Path,
    directNavigation=False,
    s: acts.examples.Sequencer = None,
):

    # Preliminaries
    rnd = acts.examples.RandomNumbers()

    # Sequencer
    s = s or acts.examples.Sequencer(
        events=10, numThreads=-1, logLevel=acts.logging.INFO
    )

    # Input
    vtxGen = acts.examples.GaussianVertexGenerator()
    vtxGen.stddev = acts.Vector4(0, 0, 0, 0)

    ptclGen = acts.examples.ParametricParticleGenerator(
        p=(1 * u.GeV, 10 * u.GeV), eta=(-2, 2), numParticles=2
    )

    g = acts.examples.EventGenerator.Generator()
    g.multiplicity = acts.examples.FixedMultiplicityGenerator()
    g.vertex = vtxGen
    g.particles = ptclGen

    evGen = acts.examples.EventGenerator(
        level=acts.logging.INFO,
        generators=[g],
        outputParticles="particles_input",
        randomNumbers=rnd,
    )
    s.addReader(evGen)

    s.addWriter(
        acts.examples.RootParticleWriter(
            level=acts.logging.INFO,
            inputParticles=evGen.config.outputParticles,
            filePath=str(outputDir / "particles.root"),
        )
    )

    # Selector
    selector = acts.examples.ParticleSelector(
        level=acts.logging.INFO,
        inputParticles=evGen.config.outputParticles,
        outputParticles="particles_selected",
    )
    s.addAlgorithm(selector)

    # Simulation
    simAlg = acts.examples.FatrasAlgorithm(
        level=acts.logging.INFO,
        inputParticles=selector.config.outputParticles,
        outputParticlesInitial="particles_initial",
        outputParticlesFinal="particles_final",
        outputSimHits="simhits",
        randomNumbers=rnd,
        trackingGeometry=trackingGeometry,
        magneticField=field,
        generateHitsOnSensitive=True,
    )
    s.addAlgorithm(simAlg)

    # Digitization
    digiCfg = acts.examples.DigitizationConfig(
        acts.examples.readDigiConfigFromJson(
            str(digiConfigFile),
        ),
        trackingGeometry=trackingGeometry,
        randomNumbers=rnd,
        inputSimHits=simAlg.config.outputSimHits,
    )
    digiAlg = acts.examples.DigitizationAlgorithm(digiCfg, acts.logging.INFO)
    s.addAlgorithm(digiAlg)

    selAlg = acts.examples.TruthSeedSelector(
        level=acts.logging.INFO,
        ptMin=500 * u.MeV,
        nHitsMin=9,
        inputParticles=simAlg.config.outputParticlesInitial,
        inputMeasurementParticlesMap=digiCfg.outputMeasurementParticlesMap,
        outputParticles="particles_seed_selected",
    )
    s.addAlgorithm(selAlg)

    inputParticles = selAlg.config.outputParticles

    smearAlg = acts.examples.ParticleSmearing(
        level=acts.logging.INFO,
        inputParticles=inputParticles,
        outputTrackParameters="smearedparameters",
        randomNumbers=rnd,
        # Gaussian sigmas to smear particle parameters
        # sigmaD0=20 * u.um,
        # sigmaD0PtA=30 * u.um,
        # sigmaD0PtB=0.3 / u.GeV,
        # sigmaZ0=20 * u.um,
        # sigmaZ0PtA=30 * u.um,
        # sigmaZ0PtB=0.3 / u.GeV,
        # sigmaPhi=1 * u.degree,
        # sigmaTheta=1 * u.degree,
        # sigmaPRel=0.01,
        # sigmaT0=1 * u.ns,
        # initialVarInflation=[1, 1, 1, 1, 1, 1],
    )
    s.addAlgorithm(smearAlg)

    truthTrkFndAlg = acts.examples.TruthTrackFinder(
        level=acts.logging.INFO,
        inputParticles=inputParticles,
        inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
        outputProtoTracks="prototracks",
    )
    s.addAlgorithm(truthTrkFndAlg)

    if directNavigation:
        srfSortAlg = acts.examples.SurfaceSortingAlgorithm(
            level=acts.logging.INFO,
            inputProtoTracks=truthTrkFndAlg.config.outputProtoTracks,
            inputSimulatedHits=simAlg.config.outputSimHits,
            inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
            outputProtoTracks="sortedprototracks",
        )
        s.addAlgorithm(srfSortAlg)
        inputProtoTracks = srfSortAlg.config.outputProtoTracks
    else:
        inputProtoTracks = truthTrkFndAlg.config.outputProtoTracks

    fitAlg = acts.examples.TrackFittingAlgorithm(
        level=acts.logging.INFO,
        inputMeasurements=digiAlg.config.outputMeasurements,
        inputSourceLinks=digiAlg.config.outputSourceLinks,
        inputProtoTracks=inputProtoTracks,
        inputInitialTrackParameters=smearAlg.config.outputTrackParameters,
        outputTrajectories="trajectories",
        directNavigation=directNavigation,
        multipleScattering=True,
        energyLoss=True,
        pickTrack=-1,
        trackingGeometry=trackingGeometry,
        dFit=acts.examples.TrackFittingAlgorithm.makeTrackFitterFunction(field),
        fit=acts.examples.TrackFittingAlgorithm.makeTrackFitterFunction(
            trackingGeometry, field
        ),
    )
    s.addAlgorithm(fitAlg)

    # Output

    s.addWriter(
        acts.examples.RootTrajectoryStatesWriter(
            level=acts.logging.INFO,
            inputTrajectories=fitAlg.config.outputTrajectories,
            inputParticles=inputParticles,
            inputSimHits=simAlg.config.outputSimHits,
            inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
            inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
            filePath=str(outputDir / "trackstates_fitter.root"),
        )
    )

    s.addWriter(
        acts.examples.RootTrajectorySummaryWriter(
            level=acts.logging.INFO,
            inputTrajectories=fitAlg.config.outputTrajectories,
            inputParticles=inputParticles,
            inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
            filePath=str(outputDir / "tracksummary_fitter.root"),
        )
    )

    s.addWriter(
        acts.examples.TrackFinderPerformanceWriter(
            level=acts.logging.INFO,
            inputProtoTracks=truthTrkFndAlg.config.outputProtoTracks,
            inputParticles=inputParticles,
            inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
            filePath=str(outputDir / "performance_track_finder.root"),
        )
    )

    s.addWriter(
        acts.examples.TrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTrajectories=fitAlg.config.outputTrajectories,
            inputParticles=inputParticles,
            inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
            filePath=str(outputDir / "performance_track_fitter.root"),
        )
    )

    return s


if "__main__" == __name__:

    # detector, trackingGeometry, _ = getOpenDataDetector()
    detector, trackingGeometry, decorators = acts.examples.GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runTruthTracking(
        trackingGeometry,
        field,
        digiConfigFile=Path(
            # "thirdparty/OpenDataDetector/config/odd-digi-smearing-config.json",
            "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
        ),
        outputDir=Path.cwd(),
    ).run()
