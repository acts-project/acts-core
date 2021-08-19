#!/usr/bin/env python3
from pathlib import Path
from typing import Optional

from acts.examples import Sequencer, GenericDetector, RootParticleReader

import acts

from acts import UnitConstants as u


def runCKFTracks(
    trackingGeometry,
    decorators,
    geometrySelection: Path,
    digiConfigFile: Path,
    field,
    outputDir: Path,
    truthSmearedSeeded=False,
    truthEstimatedSeeded=False,
    outputCsv=True,
    inputParticlePath: Optional[Path] = None,
    s=None,
):
    s = s or Sequencer(events=100, numThreads=-1)

    logger = acts.logging.getLogger("CKFExample")

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)

    if inputParticlePath is None:
        logger.info("Generating particles using particle gun")

        evGen = acts.examples.EventGenerator(
            level=acts.logging.INFO,
            generators=[
                acts.examples.EventGenerator.Generator(
                    multiplicity=acts.examples.FixedMultiplicityGenerator(n=2),
                    vertex=acts.examples.GaussianVertexGenerator(
                        stddev=acts.Vector4(0, 0, 0, 0), mean=acts.Vector4(0, 0, 0, 0)
                    ),
                    particles=acts.examples.ParametricParticleGenerator(
                        p=(1 * u.GeV, 10 * u.GeV),
                        eta=(-2, 2),
                        phi=(0, 360 * u.degree),
                        randomizeCharge=True,
                        numParticles=4,
                    ),
                )
            ],
            outputParticles="particles_input",
            randomNumbers=rnd,
        )
        s.addReader(evGen)
        inputParticles = evGen.config.outputParticles
    else:
        logger.info("Reading particles from %s", inputParticlePath.resolve())
        assert inputParticlePath.exists()
        inputParticles = "particles_read"
        s.addReader(
            RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath.resolve()),
                particleCollection=inputParticles,
                orderedEvents=False,
            )
        )

    # Selector
    selector = acts.examples.ParticleSelector(
        level=acts.logging.INFO,
        inputParticles=inputParticles,
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

    # Run the sim hits smearing
    digiCfg = acts.examples.DigitizationConfig(
        acts.examples.readDigiConfigFromJson(str(digiConfigFile)),
        trackingGeometry=trackingGeometry,
        randomNumbers=rnd,
        inputSimHits=simAlg.config.outputSimHits,
    )
    digiAlg = acts.examples.DigitizationAlgorithm(digiCfg, acts.logging.INFO)
    s.addAlgorithm(digiAlg)

    # Run the particle selection
    # The pre-selection will select truth particles satisfying provided criteria
    # from all particles read in by particle reader for further processing. It
    # has no impact on the truth hits themselves
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

    # Create starting parameters from either particle smearing or combined seed
    # finding and track parameters estimation
    if truthSmearedSeeded:
        logger.info("Using smeared truth particles for seeding")
        # Run particle smearing
        ptclSmear = acts.examples.ParticleSmearing(
            level=acts.logging.INFO,
            inputParticles=inputParticles,
            outputTrackParameters="smearedparameters",
            randomNumbers=rnd,
            # gaussian sigmas to smear particle parameters
            sigmaD0=20 * u.um,
            sigmaD0PtA=30 * u.um,
            sigmaD0PtB=0.3 / 1 * u.GeV,
            sigmaZ0=20 * u.um,
            sigmaZ0PtA=30 * u.um,
            sigmaZ0PtB=0.3 / 1 * u.GeV,
            sigmaPhi=1 * u.degree,
            sigmaTheta=1 * u.degree,
            sigmaPRel=0.01,
            sigmaT0=1 * u.ns,
            initialVarInflation=[1, 1, 1, 1, 1, 1],
        )
        outputTrackParameters = ptclSmear.config.outputTrackParameters
        s.addAlgorithm(ptclSmear)
    else:
        # Create space points
        spAlg = acts.examples.SpacePointMaker(
            level=acts.logging.INFO,
            inputSourceLinks=digiAlg.config.outputSourceLinks,
            inputMeasurements=digiAlg.config.outputMeasurements,
            outputSpacePoints="spacepoints",
            trackingGeometry=trackingGeometry,
            geometrySelection=acts.examples.readJsonGeometryList(
                str(geometrySelection)
            ),
        )
        s.addAlgorithm(spAlg)

        # Run either: truth track finding or seeding
        if truthEstimatedSeeded:
            logger.info("Using truth track finding from space points for seeding")
            # Use truth tracking
            truthTrackFinder = acts.examples.TruthTrackFinder(
                level=acts.logging.INFO,
                inputParticles=inputParticles,
                inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
                outputProtoTracks="prototracks",
            )
            s.addAlgorithm(truthTrackFinder)
            inputProtoTracks = truthTrackFinder.config.outputProtoTracks
            inputSeeds = ""
        else:
            logger.info("Using seeding")
            # Use seeding
            seeding = acts.examples.SeedingAlgorithm(
                level=acts.logging.INFO,
                inputSpacePoints=[spAlg.config.outputSpacePoints],
                outputSeeds="seeds",
                outputProtoTracks="prototracks",
            )
            s.addAlgorithm(seeding)
            inputProtoTracks = seeding.config.outputProtoTracks
            inputSeeds = seeding.config.outputSeeds

        # Write truth track finding / seeding performance
        trackFinderPerformanceWriter = acts.examples.TrackFinderPerformanceWriter(
            level=acts.logging.INFO,
            inputProtoTracks=inputProtoTracks,
            inputParticles=inputParticles,  # the original selected particles after digitization
            inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
            filePath=str(outputDir / "performance_seeding_trees.root"),
        )
        s.addWriter(trackFinderPerformanceWriter)

        # Estimate track parameters from seeds
        paramEstimation = acts.examples.TrackParamsEstimationAlgorithm(
            level=acts.logging.INFO,
            inputSeeds=inputSeeds,
            inputProtoTracks=inputProtoTracks,
            inputSpacePoints=[spAlg.config.outputSpacePoints],
            inputSourceLinks=digiCfg.outputSourceLinks,
            outputTrackParameters="estimatedparameters",
            outputProtoTracks="prototracks_estimated",
            trackingGeometry=trackingGeometry,
            magneticField=field,
            bFieldMin=0.1 * u.T,
            deltaRMax=100.0 * u.mm,
            deltaRMin=10.0 * u.mm,
            sigmaLoc0=25.0 * u.um,
            sigmaLoc1=100.0 * u.um,
            sigmaPhi=0.02 * u.degree,
            sigmaTheta=0.02 * u.degree,
            sigmaQOverP=0.1 / 1.0 * u.GeV,
            sigmaT0=1400.0 * u.s,
            initialVarInflation=[1, 1, 1, 1, 1, 1],
        )
        s.addAlgorithm(paramEstimation)
        outputTrackParameters = paramEstimation.config.outputTrackParameters

    # Setup the track finding algorithm with CKF
    # It takes all the source links created from truth hit smearing, seeds from
    # truth particle smearing and source link selection config
    trackFinder = acts.examples.TrackFindingAlgorithm(
        level=acts.logging.INFO,
        measurementSelectorCfg=acts.MeasurementSelector.Config(
            [(acts.GeometryIdentifier(), (15.0, 10))]
        ),
        inputMeasurements=digiAlg.config.outputMeasurements,
        inputSourceLinks=digiAlg.config.outputSourceLinks,
        inputInitialTrackParameters=outputTrackParameters,
        outputTrajectories="trajectories",
        findTracks=acts.examples.TrackFindingAlgorithm.makeTrackFinderFunction(
            trackingGeometry, field
        ),
    )
    s.addAlgorithm(trackFinder)

    # write track states from CKF
    trackStatesWriter = acts.examples.RootTrajectoryStatesWriter(
        level=acts.logging.INFO,
        inputTrajectories=trackFinder.config.outputTrajectories,
        # @note The full particles collection is used here to avoid lots of warnings
        # since the unselected CKF track might have a majority particle not in the
        # filtered particle collection. This could be avoided when a seperate track
        # selection algorithm is used.
        inputParticles=selector.config.outputParticles,
        inputSimHits=simAlg.config.outputSimHits,
        inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
        inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
        filePath=str(outputDir / "trackstates_ckf.root"),
        treeName="trackstates",
    )
    s.addWriter(trackStatesWriter)

    # write track summary from CKF
    trackSummaryWriter = acts.examples.RootTrajectorySummaryWriter(
        level=acts.logging.INFO,
        inputTrajectories=trackFinder.config.outputTrajectories,
        # @note The full particles collection is used here to avoid lots of warnings
        # since the unselected CKF track might have a majority particle not in the
        # filtered particle collection. This could be avoided when a seperate track
        # selection algorithm is used.
        inputParticles=selector.config.outputParticles,
        inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
        filePath=str(outputDir / "tracksummary_ckf.root"),
        treeName="tracksummary",
    )
    s.addWriter(trackSummaryWriter)

    # Write CKF performance data
    ckfPerfWriter = acts.examples.CKFPerformanceWriter(
        level=acts.logging.INFO,
        inputParticles=inputParticles,
        inputTrajectories=trackFinder.config.outputTrajectories,
        inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
        # The bottom seed could be the first, second or third hits on the truth track
        nMeasurementsMin=selAlg.config.nHitsMin - 3,
        ptMin=0.4 * u.GeV,
        filePath=str(outputDir / "performance_ckf.root"),
    )
    s.addWriter(ckfPerfWriter)

    if outputCsv:
        csv_dir = outputDir / "csv"
        csv_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Writing CSV files")
        csvMTJWriter = acts.examples.CsvMultiTrajectoryWriter(
            level=acts.logging.INFO,
            inputTrajectories=trackFinder.config.outputTrajectories,
            inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
            outputDir=str(csv_dir),
        )
        s.addWriter(csvMTJWriter)

    return s


if "__main__" == __name__:
    detector, trackingGeometry, decorators = GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    inputParticlePath = Path("particles.root")
    if not inputParticlePath.exists():
        inputParticlePath = None

    runCKFTracks(
        trackingGeometry,
        decorators,
        field=field,
        geometrySelection=Path(
            "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
        ),
        digiConfigFile=Path(
            "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
        ),
        outputCsv=True,
        truthSmearedSeeded=False,
        truthEstimatedSeeded=False,
        inputParticlePath=inputParticlePath,
        outputDir=Path.cwd(),
    ).run()
