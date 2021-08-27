#!/usr/bin/env python3
import os

import acts
import acts.examples


u = acts.UnitConstants


def runSeeding(trackingGeometry, field, outputDir, s=None):

    csv_dir = os.path.join(outputDir, "csv")
    if not os.path.exists(csv_dir):
        os.mkdir(csv_dir)

    # Input
    rnd = acts.examples.RandomNumbers(seed=42)
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
                    phi=(0, 90 * u.degree),
                    randomizeCharge=True,
                    numParticles=4,
                ),
            )
        ],
        outputParticles="particles_input",
        randomNumbers=rnd,
    )

    # Read input from input collection (e.g. Pythia8 output)

    # evGen = acts.examples.RootParticleReader(
    #     level=acts.logging.INFO,
    #     particleCollection="particles_input",
    #     inputDir="output",
    #     inputFile="pythia8_particles.root",
    # )

    # Simulation
    simAlg = acts.examples.FatrasAlgorithm(
        level=acts.logging.INFO,
        inputParticles=evGen.config.outputParticles,
        # inputParticles=evGen.config.particleCollection,
        outputParticlesInitial="particles_initial",
        outputParticlesFinal="particles_final",
        outputSimHits="simhits",
        randomNumbers=rnd,
        trackingGeometry=trackingGeometry,
        magneticField=field,
        generateHitsOnSensitive=True,
    )

    # Digitization
    digiCfg = acts.examples.DigitizationConfig(
        acts.examples.readDigiConfigFromJson(
            "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
        ),
        trackingGeometry=trackingGeometry,
        randomNumbers=rnd,
        inputSimHits=simAlg.config.outputSimHits,
    )
    digiAlg = acts.examples.DigitizationAlgorithm(digiCfg, acts.logging.INFO)

    selAlg = acts.examples.TruthSeedSelector(
        level=acts.logging.INFO,
        ptMin=1 * u.GeV,
        eta=(-2.5, 2.5),
        nHitsMin=9,
        inputParticles=simAlg.config.outputParticlesFinal,
        inputMeasurementParticlesMap=digiCfg.outputMeasurementParticlesMap,
        outputParticles="particles_selected",
    )

    inputParticles = selAlg.config.outputParticles

    spAlg = acts.examples.SpacePointMaker(
        level=acts.logging.INFO,
        inputSourceLinks=digiCfg.outputSourceLinks,
        inputMeasurements=digiCfg.outputMeasurements,
        outputSpacePoints="spacepoints",
        trackingGeometry=trackingGeometry,
        geometrySelection=acts.examples.readJsonGeometryList(
            "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
        ),
    )

    gridConfig = acts.SpacePointGridConfig(
        bFieldInZ=1.99724 * u.T,
        minPt=500 * u.MeV,
        rMax=200 * u.mm,
        zMax=2000 * u.mm,
        zMin=-2000 * u.mm,
        deltaRMax=60 * u.mm,
        cotThetaMax=7.40627,  # 2.7 eta
    )

    seedFilterConfig = acts.SeedFilterConfig(maxSeedsPerSpM=1, deltaRMin=1 * u.mm)

    seedFinderConfig = acts.SeedfinderConfig(
        rMax=gridConfig.rMax,
        deltaRMin=seedFilterConfig.deltaRMin,
        deltaRMax=gridConfig.deltaRMax,
        collisionRegionMin=-250 * u.mm,
        collisionRegionMax=250 * u.mm,
        zMin=gridConfig.zMin,
        zMax=gridConfig.zMax,
        maxSeedsPerSpM=seedFilterConfig.maxSeedsPerSpM,
        cotThetaMax=gridConfig.cotThetaMax,
        sigmaScattering=50,
        radLengthPerSeed=0.1,
        minPt=gridConfig.minPt,
        bFieldInZ=gridConfig.bFieldInZ,
        beamPos=acts.Vector2(0 * u.mm, 0 * u.mm),
        impactMax=3 * u.mm,
    )

    seedingAlg = acts.examples.SeedingAlgorithm(
        level=acts.logging.VERBOSE,
        inputSpacePoints=[spAlg.config.outputSpacePoints],
        outputSeeds="seeds",
        outputProtoTracks="prototracks",
        gridConfig=gridConfig,
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
    )

    parEstimateAlg = acts.examples.TrackParamsEstimationAlgorithm(
        level=acts.logging.VERBOSE,
        inputProtoTracks=seedingAlg.config.outputProtoTracks,
        inputSpacePoints=[spAlg.config.outputSpacePoints],
        inputSourceLinks=digiCfg.outputSourceLinks,
        outputTrackParameters="estimatedparameters",
        outputProtoTracks="prototracks_estimated",
        trackingGeometry=trackingGeometry,
        magneticField=field,
    )

    s = s or acts.examples.Sequencer(
        events=100, numThreads=-1, logLevel=acts.logging.INFO
    )

    s.addReader(evGen)
    s.addAlgorithm(simAlg)
    s.addAlgorithm(digiAlg)
    s.addAlgorithm(selAlg)
    s.addAlgorithm(spAlg)
    s.addAlgorithm(seedingAlg)
    s.addAlgorithm(parEstimateAlg)

    s.addWriter(
        acts.examples.RootParticleWriter(
            level=acts.logging.INFO,
            inputParticles=evGen.config.outputParticles,
            filePath=outputDir + "/evgen_particles.root",
        )
    )
    s.addWriter(
        acts.examples.CsvParticleWriter(
            level=acts.logging.INFO,
            inputParticles=evGen.config.outputParticles,
            outputStem="evgen_particles",
            outputDir=outputDir + "/csv",
        )
    )

    s.addWriter(
        acts.examples.RootParticleWriter(
            level=acts.logging.INFO,
            inputParticles=simAlg.config.outputParticlesFinal,
            filePath=outputDir + "/fatras_particles_final.root",
        )
    )
    s.addWriter(
        acts.examples.CsvParticleWriter(
            level=acts.logging.INFO,
            inputParticles=simAlg.config.outputParticlesFinal,
            outputStem="fatras_particles_final",
            outputDir=outputDir + "/csv",
        )
    )

    s.addWriter(
        acts.examples.RootParticleWriter(
            level=acts.logging.INFO,
            inputParticles=simAlg.config.outputParticlesInitial,
            filePath=outputDir + "/fatras_particles_initial.root",
        )
    )
    s.addWriter(
        acts.examples.CsvParticleWriter(
            level=acts.logging.INFO,
            inputParticles=simAlg.config.outputParticlesInitial,
            outputStem="fatras_particles_initial",
            outputDir=outputDir + "/csv",
        )
    )

    s.addWriter(
        acts.examples.TrackFinderPerformanceWriter(
            level=acts.logging.INFO,
            inputProtoTracks=seedingAlg.config.outputProtoTracks,
            inputParticles=inputParticles,
            inputMeasurementParticlesMap=digiCfg.outputMeasurementParticlesMap,
            filePath=outputDir + "/performance_seeding_trees.root",
        )
    )

    s.addWriter(
        acts.examples.SeedingPerformanceWriter(
            level=acts.logging.DEBUG,
            inputProtoTracks=seedingAlg.config.outputProtoTracks,
            inputParticles=inputParticles,
            inputMeasurementParticlesMap=digiCfg.outputMeasurementParticlesMap,
            filePath=outputDir + "/performance_seeding_hists.root",
        )
    )

    s.addWriter(
        acts.examples.RootTrackParameterWriter(
            level=acts.logging.VERBOSE,
            inputTrackParameters=parEstimateAlg.config.outputTrackParameters,
            inputProtoTracks=parEstimateAlg.config.outputProtoTracks,
            inputParticles=simAlg.config.outputParticlesFinal,
            inputSimHits=simAlg.config.outputSimHits,
            inputMeasurementParticlesMap=digiCfg.outputMeasurementParticlesMap,
            inputMeasurementSimHitsMap=digiCfg.outputMeasurementSimHitsMap,
            filePath=outputDir + "/estimatedparams.root",
            treeName="estimatedparams",
        )
    )

    return s


if "__main__" == __name__:
    from common import getOpenDataDetector

    # detector, trackingGeometry, _ = getOpenDataDetector()
    detector, trackingGeometry, _ = acts.examples.GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runSeeding(trackingGeometry, field, outputDir=os.getcwd()).run()
