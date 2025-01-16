from pathlib import Path
from typing import Optional, Union, List
from enum import Enum
from collections import namedtuple

import acts
import acts.examples

u = acts.UnitConstants

SeedingAlgorithm = Enum(
    "SeedingAlgorithm",
    "Default TruthSmeared TruthEstimated Orthogonal HoughTransform Gbts Hashing",
)

TrackSmearingSigmas = namedtuple(
    "TrackSmearingSigmas",
    [
        "loc0",
        "loc0PtA",
        "loc0PtB",
        "loc1",
        "loc1PtA",
        "loc1PtB",
        "time",
        "phi",
        "theta",
        "ptRel",
    ],
    defaults=[None] * 10,
)

SeedFinderConfigArg = namedtuple(
    "SeedFinderConfig",
    [
        "maxSeedsPerSpM",
        "cotThetaMax",
        "sigmaScattering",
        "radLengthPerSeed",
        "minPt",
        "impactMax",
        "deltaPhiMax",
        "interactionPointCut",
        "deltaZMax",
        "maxPtScattering",
        "zBinEdges",
        "zBinsCustomLooping",
        "rRangeMiddleSP",
        "useVariableMiddleSPRange",
        "binSizeR",
        "seedConfirmation",
        "centralSeedConfirmationRange",
        "forwardSeedConfirmationRange",
        "deltaR",  # (min,max)
        "deltaRBottomSP",  # (min,max)
        "deltaRTopSP",  # (min,max)
        "deltaRMiddleSPRange",  # (min,max)
        "collisionRegion",  # (min,max)
        "r",  # (min,max)
        "z",  # (min,max)
    ],
    defaults=[None] * 18 + [(None, None)] * 7,
)
SeedFinderOptionsArg = namedtuple(
    "SeedFinderOptions", ["beamPos", "bFieldInZ"], defaults=[(None, None), None]
)

SeedFilterConfigArg = namedtuple(
    "SeedFilterConfig",
    [
        "impactWeightFactor",
        "zOriginWeightFactor",
        "compatSeedWeight",
        "compatSeedLimit",
        "numSeedIncrement",
        "seedWeightIncrement",
        "seedConfirmation",
        "maxSeedsPerSpMConf",
        "maxQualitySeedsPerSpMConf",
        "useDeltaRorTopRadius",
        "deltaRMin",
    ],
    defaults=[None] * 11,
)

SpacePointGridConfigArg = namedtuple(
    "SeedGridConfig",
    [
        "rMax",
        "zBinEdges",
        "phiBinDeflectionCoverage",
        "impactMax",
        "deltaRMax",
        "maxPhiBins",
        "phi",  # (min,max)
    ],
    defaults=[None] * 6 + [(None, None)] * 1,
)

SeedingAlgorithmConfigArg = namedtuple(
    "SeedingAlgorithmConfig",
    [
        "allowSeparateRMax",
        "zBinNeighborsTop",
        "zBinNeighborsBottom",
        "numPhiNeighbors",
        "useExtraCuts",
    ],
    defaults=[None] * 5,
)

TruthEstimatedSeedingAlgorithmConfigArg = namedtuple(
    "TruthSeederConfig",
    [
        "deltaR",  # (min,max)
    ],
    defaults=[(None, None)],
)

TrackSelectorConfig = namedtuple(
    "TrackSelectorConfig",
    [
        "loc0",
        "loc1",
        "time",
        "eta",
        "absEta",
        "pt",
        "phi",
        "nMeasurementsMin",
        "maxHoles",
        "maxOutliers",
        "maxHolesAndOutliers",
        "maxSharedHits",
        "maxChi2",
        "nMeasurementsGroupMin",
        "requireReferenceSurface",
    ],
    defaults=[(None, None)] * 7 + [None] * 8,
)


def trackSelectorDefaultKWArgs(c):
    """
    Encapsulate this boilerplate code into a function so different uses do not get out of sync
    """
    return acts.examples.defaultKWArgs(
        loc0Min=c.loc0[0],
        loc0Max=c.loc0[1],
        loc1Min=c.loc1[0],
        loc1Max=c.loc1[1],
        timeMin=c.time[0],
        timeMax=c.time[1],
        phiMin=c.phi[0],
        phiMax=c.phi[1],
        etaMin=c.eta[0],
        etaMax=c.eta[1],
        absEtaMin=c.absEta[0],
        absEtaMax=c.absEta[1],
        ptMin=c.pt[0],
        ptMax=c.pt[1],
        minMeasurements=c.nMeasurementsMin,
        maxHoles=c.maxHoles,
        maxOutliers=c.maxOutliers,
        maxHolesAndOutliers=c.maxHolesAndOutliers,
        maxSharedHits=c.maxSharedHits,
        maxChi2=c.maxChi2,
        measurementCounter=c.nMeasurementsGroupMin,
        requireReferenceSurface=c.requireReferenceSurface,
    )


CkfConfig = namedtuple(
    "CkfConfig",
    [
        "chi2CutOffMeasurement",
        "chi2CutOffOutlier",
        "numMeasurementsCutOff",
        "maxSteps",
        "seedDeduplication",
        "stayOnSeed",
        "pixelVolumes",
        "stripVolumes",
        "maxPixelHoles",
        "maxStripHoles",
        "trimTracks",
        "constrainToVolumes",
        "endOfWorldVolumes",
    ],
    defaults=[
        15.0,
        25.0,
        10,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    ],
)

AmbiguityResolutionConfig = namedtuple(
    "AmbiguityResolutionConfig",
    ["maximumSharedHits", "nMeasurementsMin", "maximumIterations"],
    defaults=[None] * 3,
)

ScoreBasedAmbiguityResolutionConfig = namedtuple(
    "ScoreBasedAmbiguityResolutionConfig",
    [
        "minScore",
        "minScoreSharedTracks",
        "maxShared",
        "maxSharedTracksPerMeasurement",
        "pTMax",
        "pTMin",
        "phiMax",
        "phiMin",
        "etaMax",
        "etaMin",
        "useAmbiguityFunction",
    ],
    defaults=[None] * 11,
)

AmbiguityResolutionMLConfig = namedtuple(
    "AmbiguityResolutionMLConfig",
    ["maximumSharedHits", "nMeasurementsMin", "maximumIterations"],
    defaults=[None] * 3,
)

SeedFilterMLDBScanConfig = namedtuple(
    "SeedFilterMLDBScanConfig",
    ["epsilonDBScan", "minPointsDBScan", "minSeedScore"],
    defaults=[None] * 3,
)

HashingTrainingConfigArg = namedtuple(
    "HashingTrainingConfig",
    ["annoySeed", "f"],
    defaults=[None] * 2,
)

HashingAlgorithmConfigArg = namedtuple(
    "HashingAlgorithmConfig",
    ["bucketSize", "zBins", "phiBins"],
    defaults=[None] * 3,
)


class VertexFinder(Enum):
    Truth = (1,)
    AMVF = (2,)
    Iterative = (3,)


@acts.examples.NamedTypeArgs(
    seedingAlgorithm=SeedingAlgorithm,
    trackSmearingSigmas=TrackSmearingSigmas,
    seedFinderConfigArg=SeedFinderConfigArg,
    seedFinderOptionsArg=SeedFinderOptionsArg,
    seedFilterConfigArg=SeedFilterConfigArg,
    spacePointGridConfigArg=SpacePointGridConfigArg,
    seedingAlgorithmConfigArg=SeedingAlgorithmConfigArg,
    hashingTrainingConfigArg=HashingTrainingConfigArg,
    hashingAlgorithmConfigArg=HashingAlgorithmConfigArg,
    truthEstimatedSeedingAlgorithmConfigArg=TruthEstimatedSeedingAlgorithmConfigArg,
    logLevel=acts.logging.Level,
)
def addSeeding(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    geoSelectionConfigFile: Optional[Union[Path, str]] = None,
    layerMappingConfigFile: Optional[Union[Path, str]] = None,
    connector_inputConfigFile: Optional[Union[Path, str]] = None,
    seedingAlgorithm: SeedingAlgorithm = SeedingAlgorithm.Default,
    trackSmearingSigmas: TrackSmearingSigmas = TrackSmearingSigmas(),
    initialSigmas: Optional[list] = None,
    initialSigmaPtRel: Optional[float] = None,
    initialVarInflation: Optional[list] = None,
    seedFinderConfigArg: SeedFinderConfigArg = SeedFinderConfigArg(),
    seedFinderOptionsArg: SeedFinderOptionsArg = SeedFinderOptionsArg(),
    seedFilterConfigArg: SeedFilterConfigArg = SeedFilterConfigArg(),
    spacePointGridConfigArg: SpacePointGridConfigArg = SpacePointGridConfigArg(),
    seedingAlgorithmConfigArg: SeedingAlgorithmConfigArg = SeedingAlgorithmConfigArg(),
    houghTransformConfig: acts.examples.HoughTransformSeeder.Config = acts.examples.HoughTransformSeeder.Config(),
    hashingTrainingConfigArg: Optional[
        HashingTrainingConfigArg
    ] = HashingTrainingConfigArg(),
    hashingAlgorithmConfigArg: Optional[
        HashingAlgorithmConfigArg
    ] = HashingAlgorithmConfigArg(),
    truthEstimatedSeedingAlgorithmConfigArg: TruthEstimatedSeedingAlgorithmConfigArg = TruthEstimatedSeedingAlgorithmConfigArg(),
    particleHypothesis: Optional[
        acts.ParticleHypothesis
    ] = acts.ParticleHypothesis.pion,
    inputParticles: str = "particles",
    selectedParticles: str = "particles_selected",
    outputDirRoot: Optional[Union[Path, str]] = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
    rnd: Optional[acts.examples.RandomNumbers] = None,
) -> None:
    """This function steers the seeding
    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addSeeding)
    trackingGeometry : tracking geometry
    field : magnetic field
    geoSelectionConfigFile : Path|str, path, None
        Json file for space point geometry selection. Not required for SeedingAlgorithm.TruthSmeared.
    seedingAlgorithm : SeedingAlgorithm, Default
        seeding algorithm to use: one of Default (no truth information used), TruthSmeared, TruthEstimated
    trackSmearingSigmas : TrackSmearingSigmas(loc0, loc0PtA, loc0PtB, loc1, loc1PtA, loc1PtB, time, phi, theta, ptRel)
        TrackSmearing configuration.
        Defaults specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TrackParameterSmearing.hpp
    initialSigmas : list
        Sets the initial covariance matrix diagonal. This is ignored in case of TruthSmearing.
        Defaults specified in Examples/Algorithms/TrackFinding/include/ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp
    initialVarInflation : list
        List of 6 scale factors to inflate the initial covariance matrix
        Defaults (all 1) specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TrackParameterSmearing.hpp
    seedFinderConfigArg : SeedFinderConfigArg(maxSeedsPerSpM, cotThetaMax, sigmaScattering, radLengthPerSeed, minPt, impactMax, deltaPhiMax, interactionPointCut, deltaZMax, maxPtScattering, zBinEdges, zBinsCustomLooping, rRangeMiddleSP, useVariableMiddleSPRange, binSizeR, seedConfirmation, centralSeedConfirmationRange, forwardSeedConfirmationRange, deltaR, deltaRBottomSP, deltaRTopSP, deltaRMiddleSPRange, collisionRegion, r, z)
        SeedFinderConfig settings. deltaR, deltaRBottomSP, deltaRTopSP, deltaRMiddleSPRange, collisionRegion, r, z.
        Defaults specified in Core/include/Acts/Seeding/SeedFinderConfig.hpp
    seedFinderOptionsArg :  SeedFinderOptionsArg(bFieldInZ, beamPos)
        Defaults specified in Core/include/Acts/Seeding/SeedFinderConfig.hpp
    seedFilterConfigArg : SeedFilterConfigArg(compatSeedWeight, compatSeedLimit, numSeedIncrement, seedWeightIncrement, seedConfirmation, maxSeedsPerSpMConf, maxQualitySeedsPerSpMConf, useDeltaRorTopRadius)
                                Defaults specified in Core/include/Acts/Seeding/SeedFilterConfig.hpp
    spacePointGridConfigArg : SpacePointGridConfigArg(rMax, zBinEdges, phiBinDeflectionCoverage, phi, maxPhiBins, impactMax)
                                SpacePointGridConfigArg settings. phi is specified as a tuple of (min,max).
        Defaults specified in Core/include/Acts/Seeding/SpacePointGrid.hpp
    seedingAlgorithmConfigArg : SeedingAlgorithmConfigArg(allowSeparateRMax, zBinNeighborsTop, zBinNeighborsBottom, numPhiNeighbors, useExtraCuts)
                                Defaults specified in Examples/Algorithms/TrackFinding/include/ActsExamples/TrackFinding/SeedingAlgorithm.hpp
    hashingTrainingConfigArg : HashingTrainingConfigArg(annoySeed, f)
                                Defaults specified in Plugins/Hashing/include/Acts/Plugins/Hashing/HashingTrainingConfig.hpp
    hashingAlgorithmConfigArg : HashingAlgorithmConfigArg(bucketSize, zBins, phiBins)
                                Defaults specified in Plugins/Hashing/include/Acts/Plugins/Hashing/HashingAlgorithmConfig.hpp
    truthEstimatedSeedingAlgorithmConfigArg : TruthEstimatedSeedingAlgorithmConfigArg(deltaR)
        Currently only deltaR=(min,max) range specified here.
    particleHypothesis : Optional[acts.ParticleHypothesis]
        The hypothesis used for track finding. Defaults to pion.
    inputParticles : str, "particles"
        input particles name in the WhiteBoard
    selectedParticles : str, "particles_selected"
        selected particles name in the WhiteBoard
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    logLevel : acts.logging.Level, None
        logging level to override setting given in `s`
    rnd : RandomNumbers, None
        random number generator. Only used by SeedingAlgorithm.TruthSmeared.
    """

    logLevel = acts.examples.defaultLogging(s, logLevel)()
    logger = acts.logging.getLogger("addSeeding")
    logger.setLevel(logLevel)

    # Create starting parameters from either particle smearing or combined seed
    # finding and track parameters estimation
    if seedingAlgorithm == SeedingAlgorithm.TruthSmeared:
        logger.info("Using smeared truth particles for seeding")
        addTruthSmearedSeeding(
            s=s,
            rnd=rnd,
            selectedParticles=selectedParticles,
            trackSmearingSigmas=trackSmearingSigmas,
            initialSigmas=initialSigmas,
            initialSigmaPtRel=initialSigmaPtRel,
            initialVarInflation=initialVarInflation,
            particleHypothesis=particleHypothesis,
            logLevel=logLevel,
        )
    else:
        spacePoints = addSpacePointsMaking(
            s, trackingGeometry, geoSelectionConfigFile, logLevel
        )
        # Run either: truth track finding or seeding
        if seedingAlgorithm == SeedingAlgorithm.TruthEstimated:
            logger.info("Using truth track finding from space points for seeding")
            seeds = addTruthEstimatedSeeding(
                s,
                spacePoints,
                selectedParticles,
                truthEstimatedSeedingAlgorithmConfigArg,
                logLevel,
            )
        elif seedingAlgorithm == SeedingAlgorithm.Default:
            logger.info("Using default seeding")
            seeds = addStandardSeeding(
                s,
                spacePoints,
                seedingAlgorithmConfigArg,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                spacePointGridConfigArg,
                logLevel,
            )
        elif seedingAlgorithm == SeedingAlgorithm.Orthogonal:
            logger.info("Using orthogonal seeding")
            seeds = addOrthogonalSeeding(
                s,
                spacePoints,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                logLevel,
            )
        elif seedingAlgorithm == SeedingAlgorithm.HoughTransform:
            logger.info("Using Hough Transform seeding")
            houghTransformConfig.inputSpacePoints = [spacePoints]
            houghTransformConfig.inputMeasurements = "measurements"
            houghTransformConfig.outputProtoTracks = "prototracks"
            houghTransformConfig.outputSeeds = "seeds"
            houghTransformConfig.trackingGeometry = trackingGeometry
            seeds = addHoughTransformSeeding(s, houghTransformConfig, logLevel)
        elif seedingAlgorithm == SeedingAlgorithm.Gbts:
            logger.info("Using Gbts seeding")
            # output of algs changed, only one output now
            seeds = addGbtsSeeding(
                s,
                spacePoints,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                trackingGeometry,
                logLevel,
                layerMappingConfigFile,
                geoSelectionConfigFile,
                connector_inputConfigFile,
            )
        elif seedingAlgorithm == SeedingAlgorithm.Hashing:
            logger.info("Using Hashing seeding")
            seeds, buckets = addHashingSeeding(
                s,
                spacePoints,
                seedingAlgorithmConfigArg,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                spacePointGridConfigArg,
                hashingTrainingConfigArg,
                hashingAlgorithmConfigArg,
                logLevel,
            )
        else:
            logger.fatal("unknown seedingAlgorithm %s", seedingAlgorithm)

        parEstimateAlg = acts.examples.TrackParamsEstimationAlgorithm(
            level=logLevel,
            inputSeeds=seeds,
            outputTrackParameters="estimatedparameters",
            outputSeeds="estimatedseeds",
            trackingGeometry=trackingGeometry,
            magneticField=field,
            **acts.examples.defaultKWArgs(
                initialSigmas=initialSigmas,
                initialSigmaPtRel=initialSigmaPtRel,
                initialVarInflation=initialVarInflation,
                particleHypothesis=particleHypothesis,
            ),
        )
        s.addAlgorithm(parEstimateAlg)

        prototracks = "seed-prototracks"
        s.addAlgorithm(
            acts.examples.SeedsToPrototracks(
                level=logLevel,
                inputSeeds=seeds,
                outputProtoTracks=prototracks,
            )
        )

        if outputDirRoot is not None:
            addSeedPerformanceWriters(
                s,
                outputDirRoot,
                seeds,
                prototracks,
                selectedParticles,
                inputParticles,
                parEstimateAlg.config.outputTrackParameters,
                logLevel,
            )

        if outputDirCsv is not None:
            outputDirCsv = Path(outputDirCsv)

            if not outputDirCsv.exists():
                outputDirCsv.mkdir()

            csvSeedWriter = acts.examples.CsvSeedWriter(
                level=logLevel,
                inputTrackParameters=parEstimateAlg.config.outputTrackParameters,
                inputSimSeeds=seeds,
                inputSimHits="simhits",
                inputMeasurementParticlesMap="measurement_particles_map",
                inputMeasurementSimHitsMap="measurement_simhits_map",
                outputDir=str(outputDirCsv),
                fileName=str(f"seed.csv"),
            )
            s.addWriter(csvSeedWriter)

            if seedingAlgorithm == SeedingAlgorithm.Hashing:
                s.addWriter(
                    acts.examples.CsvSpacePointsBucketWriter(
                        level=logLevel,
                        inputBuckets=buckets,
                        outputDir=str(outputDirCsv),
                    )
                )

    return s


def addTruthSmearedSeeding(
    s: acts.examples.Sequencer,
    rnd: Optional[acts.examples.RandomNumbers],
    selectedParticles: str,
    trackSmearingSigmas: TrackSmearingSigmas,
    initialSigmas: Optional[List[float]],
    initialSigmaPtRel: Optional[float],
    initialVarInflation: Optional[List[float]],
    particleHypothesis: Optional[acts.ParticleHypothesis],
    logLevel: acts.logging.Level = None,
):
    """adds algorithm that would mimic detector response uncertainties for truth seeding
    For parameters description see addSeeding
    """

    rnd = rnd or acts.examples.RandomNumbers(seed=42)

    trkParamExtractor = acts.examples.ParticleTrackParamExtractor(
        level=logLevel,
        inputParticles=selectedParticles,
        outputTrackParameters="trueparameters",
    )
    s.addAlgorithm(trkParamExtractor)

    # Smearing track parameters
    trkSmear = acts.examples.TrackParameterSmearing(
        level=logLevel,
        inputTrackParameters=trkParamExtractor.config.outputTrackParameters,
        outputTrackParameters="estimatedparameters",
        randomNumbers=rnd,
        # gaussian sigmas to smear particle parameters
        **acts.examples.defaultKWArgs(
            sigmaLoc0=trackSmearingSigmas.loc0,
            sigmaLoc0PtA=trackSmearingSigmas.loc0PtA,
            sigmaLoc0PtB=trackSmearingSigmas.loc0PtB,
            sigmaLoc1=trackSmearingSigmas.loc1,
            sigmaLoc1PtA=trackSmearingSigmas.loc1PtA,
            sigmaLoc1PtB=trackSmearingSigmas.loc1PtB,
            sigmaTime=trackSmearingSigmas.time,
            sigmaPhi=trackSmearingSigmas.phi,
            sigmaTheta=trackSmearingSigmas.theta,
            sigmaPtRel=trackSmearingSigmas.ptRel,
            initialSigmas=initialSigmas,
            initialSigmaPtRel=initialSigmaPtRel,
            initialVarInflation=initialVarInflation,
            particleHypothesis=particleHypothesis,
        ),
    )
    s.addAlgorithm(trkSmear)

    truthTrkFndAlg = acts.examples.TruthTrackFinder(
        level=logLevel,
        inputParticles=selectedParticles,
        inputParticleMeasurementsMap="particle_measurements_map",
        outputProtoTracks="truth_particle_tracks",
    )
    s.addAlgorithm(truthTrkFndAlg)


def addTruthEstimatedSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    inputParticles: str,
    TruthEstimatedSeedingAlgorithmConfigArg: TruthEstimatedSeedingAlgorithmConfigArg,
    logLevel: acts.logging.Level = None,
):
    """adds truth seeding
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()

    truthSeeding = acts.examples.TruthSeedingAlgorithm(
        level=logLevel,
        inputParticles=inputParticles,
        inputParticleMeasurementsMap="particle_measurements_map",
        inputSpacePoints=[spacePoints],
        outputParticles="truth_seeded_particles",
        outputProtoTracks="truth_particle_tracks",
        outputSeeds="seeds",
        **acts.examples.defaultKWArgs(
            deltaRMin=TruthEstimatedSeedingAlgorithmConfigArg.deltaR[0],
            deltaRMax=TruthEstimatedSeedingAlgorithmConfigArg.deltaR[1],
        ),
    )
    sequence.addAlgorithm(truthSeeding)

    return truthSeeding.config.outputSeeds


def addSpacePointsMaking(
    sequence: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    geoSelectionConfigFile: Union[Path, str],
    logLevel: acts.logging.Level = None,
):
    """adds space points making
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    spAlg = acts.examples.SpacePointMaker(
        level=logLevel,
        inputMeasurements="measurements",
        outputSpacePoints="spacepoints",
        trackingGeometry=trackingGeometry,
        geometrySelection=acts.examples.readJsonGeometryList(
            str(geoSelectionConfigFile)
        ),
    )
    sequence.addAlgorithm(spAlg)
    return spAlg.config.outputSpacePoints


def addStandardSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    seedingAlgorithmConfigArg: SeedingAlgorithmConfigArg,
    seedFinderConfigArg: SeedFinderConfigArg,
    seedFinderOptionsArg: SeedFinderOptionsArg,
    seedFilterConfigArg: SeedFilterConfigArg,
    spacePointGridConfigArg: SpacePointGridConfigArg,
    logLevel: acts.logging.Level = None,
):
    """adds standard seeding
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()

    seedFinderConfig = acts.SeedFinderConfig(
        **acts.examples.defaultKWArgs(
            rMin=seedFinderConfigArg.r[0],
            rMax=seedFinderConfigArg.r[1],
            deltaRMin=seedFinderConfigArg.deltaR[0],
            deltaRMax=seedFinderConfigArg.deltaR[1],
            deltaRMinTopSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRTopSP[0] is None
                else seedFinderConfigArg.deltaRTopSP[0]
            ),
            deltaRMaxTopSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRTopSP[1] is None
                else seedFinderConfigArg.deltaRTopSP[1]
            ),
            deltaRMinBottomSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRBottomSP[0] is None
                else seedFinderConfigArg.deltaRBottomSP[0]
            ),
            deltaRMaxBottomSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRBottomSP[1] is None
                else seedFinderConfigArg.deltaRBottomSP[1]
            ),
            deltaRMiddleMinSPRange=seedFinderConfigArg.deltaRMiddleSPRange[0],
            deltaRMiddleMaxSPRange=seedFinderConfigArg.deltaRMiddleSPRange[1],
            collisionRegionMin=seedFinderConfigArg.collisionRegion[0],
            collisionRegionMax=seedFinderConfigArg.collisionRegion[1],
            zMin=seedFinderConfigArg.z[0],
            zMax=seedFinderConfigArg.z[1],
            maxSeedsPerSpM=seedFinderConfigArg.maxSeedsPerSpM,
            cotThetaMax=seedFinderConfigArg.cotThetaMax,
            sigmaScattering=seedFinderConfigArg.sigmaScattering,
            radLengthPerSeed=seedFinderConfigArg.radLengthPerSeed,
            minPt=seedFinderConfigArg.minPt,
            impactMax=seedFinderConfigArg.impactMax,
            interactionPointCut=seedFinderConfigArg.interactionPointCut,
            deltaZMax=seedFinderConfigArg.deltaZMax,
            maxPtScattering=seedFinderConfigArg.maxPtScattering,
            zBinEdges=seedFinderConfigArg.zBinEdges,
            zBinsCustomLooping=seedFinderConfigArg.zBinsCustomLooping,
            rRangeMiddleSP=seedFinderConfigArg.rRangeMiddleSP,
            useVariableMiddleSPRange=seedFinderConfigArg.useVariableMiddleSPRange,
            binSizeR=seedFinderConfigArg.binSizeR,
            seedConfirmation=seedFinderConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfigArg.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfigArg.forwardSeedConfirmationRange,
        ),
    )
    seedFinderOptions = acts.SeedFinderOptions(
        **acts.examples.defaultKWArgs(
            beamPos=(
                acts.Vector2(0.0, 0.0)
                if seedFinderOptionsArg.beamPos == (None, None)
                else acts.Vector2(
                    seedFinderOptionsArg.beamPos[0], seedFinderOptionsArg.beamPos[1]
                )
            ),
            bFieldInZ=seedFinderOptionsArg.bFieldInZ,
        )
    )
    seedFilterConfig = acts.SeedFilterConfig(
        **acts.examples.defaultKWArgs(
            maxSeedsPerSpM=seedFinderConfig.maxSeedsPerSpM,
            deltaRMin=(
                seedFinderConfig.deltaRMin
                if seedFilterConfigArg.deltaRMin is None
                else seedFilterConfigArg.deltaRMin
            ),
            impactWeightFactor=seedFilterConfigArg.impactWeightFactor,
            zOriginWeightFactor=seedFilterConfigArg.zOriginWeightFactor,
            compatSeedWeight=seedFilterConfigArg.compatSeedWeight,
            compatSeedLimit=seedFilterConfigArg.compatSeedLimit,
            numSeedIncrement=seedFilterConfigArg.numSeedIncrement,
            seedWeightIncrement=seedFilterConfigArg.seedWeightIncrement,
            seedConfirmation=seedFilterConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfig.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfig.forwardSeedConfirmationRange,
            maxSeedsPerSpMConf=seedFilterConfigArg.maxSeedsPerSpMConf,
            maxQualitySeedsPerSpMConf=seedFilterConfigArg.maxQualitySeedsPerSpMConf,
            useDeltaRorTopRadius=seedFilterConfigArg.useDeltaRorTopRadius,
        )
    )

    gridConfig = acts.SpacePointGridConfig(
        **acts.examples.defaultKWArgs(
            minPt=seedFinderConfig.minPt,
            rMax=(
                seedFinderConfig.rMax
                if spacePointGridConfigArg.rMax is None
                else spacePointGridConfigArg.rMax
            ),
            zMax=seedFinderConfig.zMax,
            zMin=seedFinderConfig.zMin,
            deltaRMax=(
                seedFinderConfig.deltaRMax
                if spacePointGridConfigArg.deltaRMax is None
                else spacePointGridConfigArg.deltaRMax
            ),
            cotThetaMax=seedFinderConfig.cotThetaMax,
            phiMin=spacePointGridConfigArg.phi[0],
            phiMax=spacePointGridConfigArg.phi[1],
            maxPhiBins=spacePointGridConfigArg.maxPhiBins,
            impactMax=spacePointGridConfigArg.impactMax,
            zBinEdges=spacePointGridConfigArg.zBinEdges,
            phiBinDeflectionCoverage=spacePointGridConfigArg.phiBinDeflectionCoverage,
        )
    )

    gridOptions = acts.SpacePointGridOptions(
        **acts.examples.defaultKWArgs(
            bFieldInZ=seedFinderOptions.bFieldInZ,
        )
    )

    seedingAlg = acts.examples.SeedingAlgorithm(
        level=logLevel,
        inputSpacePoints=[spacePoints],
        outputSeeds="seeds",
        **acts.examples.defaultKWArgs(
            allowSeparateRMax=seedingAlgorithmConfigArg.allowSeparateRMax,
            zBinNeighborsTop=seedingAlgorithmConfigArg.zBinNeighborsTop,
            zBinNeighborsBottom=seedingAlgorithmConfigArg.zBinNeighborsBottom,
            numPhiNeighbors=seedingAlgorithmConfigArg.numPhiNeighbors,
            useExtraCuts=seedingAlgorithmConfigArg.useExtraCuts,
        ),
        gridConfig=gridConfig,
        gridOptions=gridOptions,
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
        seedFinderOptions=seedFinderOptions,
    )
    sequence.addAlgorithm(seedingAlg)

    return seedingAlg.config.outputSeeds


def addOrthogonalSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    seedFinderConfigArg: SeedFinderConfigArg,
    seedFinderOptionsArg: SeedFinderOptionsArg,
    seedFilterConfigArg: SeedFilterConfigArg,
    logLevel: acts.logging.Level = None,
):
    """adds orthogonal seeding algorithm
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    seedFinderConfig = acts.SeedFinderOrthogonalConfig(
        **acts.examples.defaultKWArgs(
            rMin=seedFinderConfigArg.r[0],
            rMax=seedFinderConfigArg.r[1],
            deltaRMinTopSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRTopSP[0] is None
                else seedFinderConfigArg.deltaRTopSP[0]
            ),
            deltaRMaxTopSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRTopSP[1] is None
                else seedFinderConfigArg.deltaRTopSP[1]
            ),
            deltaRMinBottomSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRBottomSP[0] is None
                else seedFinderConfigArg.deltaRBottomSP[0]
            ),
            deltaRMaxBottomSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRBottomSP[1] is None
                else seedFinderConfigArg.deltaRBottomSP[1]
            ),
            collisionRegionMin=seedFinderConfigArg.collisionRegion[0],
            collisionRegionMax=seedFinderConfigArg.collisionRegion[1],
            zMin=seedFinderConfigArg.z[0],
            zMax=seedFinderConfigArg.z[1],
            maxSeedsPerSpM=seedFinderConfigArg.maxSeedsPerSpM,
            cotThetaMax=seedFinderConfigArg.cotThetaMax,
            sigmaScattering=seedFinderConfigArg.sigmaScattering,
            radLengthPerSeed=seedFinderConfigArg.radLengthPerSeed,
            minPt=seedFinderConfigArg.minPt,
            impactMax=seedFinderConfigArg.impactMax,
            deltaPhiMax=seedFinderConfigArg.deltaPhiMax,
            interactionPointCut=seedFinderConfigArg.interactionPointCut,
            deltaZMax=seedFinderConfigArg.deltaZMax,
            maxPtScattering=seedFinderConfigArg.maxPtScattering,
            rRangeMiddleSP=seedFinderConfigArg.rRangeMiddleSP,
            useVariableMiddleSPRange=seedFinderConfigArg.useVariableMiddleSPRange,
            seedConfirmation=seedFinderConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfigArg.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfigArg.forwardSeedConfirmationRange,
        ),
    )
    seedFinderOptions = acts.SeedFinderOptions(
        **acts.examples.defaultKWArgs(
            beamPos=(
                acts.Vector2(0.0, 0.0)
                if seedFinderOptionsArg.beamPos == (None, None)
                else acts.Vector2(
                    seedFinderOptionsArg.beamPos[0], seedFinderOptionsArg.beamPos[1]
                )
            ),
            bFieldInZ=seedFinderOptionsArg.bFieldInZ,
        )
    )
    seedFilterConfig = acts.SeedFilterConfig(
        **acts.examples.defaultKWArgs(
            maxSeedsPerSpM=seedFinderConfig.maxSeedsPerSpM,
            deltaRMin=(
                seedFinderConfigArg.deltaR[0]
                if seedFilterConfigArg.deltaRMin is None
                else seedFilterConfigArg.deltaRMin
            ),
            impactWeightFactor=seedFilterConfigArg.impactWeightFactor,
            zOriginWeightFactor=seedFilterConfigArg.zOriginWeightFactor,
            compatSeedWeight=seedFilterConfigArg.compatSeedWeight,
            compatSeedLimit=seedFilterConfigArg.compatSeedLimit,
            numSeedIncrement=seedFilterConfigArg.numSeedIncrement,
            seedWeightIncrement=seedFilterConfigArg.seedWeightIncrement,
            seedConfirmation=seedFilterConfigArg.seedConfirmation,
            maxSeedsPerSpMConf=seedFilterConfigArg.maxSeedsPerSpMConf,
            maxQualitySeedsPerSpMConf=seedFilterConfigArg.maxQualitySeedsPerSpMConf,
            useDeltaRorTopRadius=seedFilterConfigArg.useDeltaRorTopRadius,
        )
    )
    seedingAlg = acts.examples.SeedingOrthogonalAlgorithm(
        level=logLevel,
        inputSpacePoints=[spacePoints],
        outputSeeds="seeds",
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
        seedFinderOptions=seedFinderOptions,
    )
    sequence.addAlgorithm(seedingAlg)

    return seedingAlg.config.outputSeeds


def addHashingSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    seedingAlgorithmConfigArg: SeedingAlgorithmConfigArg,
    seedFinderConfigArg: SeedFinderConfigArg,
    seedFinderOptionsArg: SeedFinderOptionsArg,
    seedFilterConfigArg: SeedFilterConfigArg,
    spacePointGridConfigArg: SpacePointGridConfigArg,
    hashingTrainingConfigArg: HashingTrainingConfigArg,
    hashingAlgorithmConfigArg: HashingAlgorithmConfigArg,
    logLevel: acts.logging.Level = None,
):
    """adds Hashing seeding
    For parameters description see addSeeding docstring
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    from acts.examples.hashing import SeedingAlgorithmHashing

    # Same configuration than the standard seeding
    seedFinderConfig = acts.SeedFinderConfig(
        **acts.examples.defaultKWArgs(
            rMin=seedFinderConfigArg.r[0],
            rMax=seedFinderConfigArg.r[1],
            deltaRMin=seedFinderConfigArg.deltaR[0],
            deltaRMax=seedFinderConfigArg.deltaR[1],
            deltaRMinTopSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRTopSP[0] is None
                else seedFinderConfigArg.deltaRTopSP[0]
            ),
            deltaRMaxTopSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRTopSP[1] is None
                else seedFinderConfigArg.deltaRTopSP[1]
            ),
            deltaRMinBottomSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRBottomSP[0] is None
                else seedFinderConfigArg.deltaRBottomSP[0]
            ),
            deltaRMaxBottomSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRBottomSP[1] is None
                else seedFinderConfigArg.deltaRBottomSP[1]
            ),
            deltaRMiddleMinSPRange=seedFinderConfigArg.deltaRMiddleSPRange[0],
            deltaRMiddleMaxSPRange=seedFinderConfigArg.deltaRMiddleSPRange[1],
            collisionRegionMin=seedFinderConfigArg.collisionRegion[0],
            collisionRegionMax=seedFinderConfigArg.collisionRegion[1],
            zMin=seedFinderConfigArg.z[0],
            zMax=seedFinderConfigArg.z[1],
            maxSeedsPerSpM=seedFinderConfigArg.maxSeedsPerSpM,
            cotThetaMax=seedFinderConfigArg.cotThetaMax,
            sigmaScattering=seedFinderConfigArg.sigmaScattering,
            radLengthPerSeed=seedFinderConfigArg.radLengthPerSeed,
            minPt=seedFinderConfigArg.minPt,
            impactMax=seedFinderConfigArg.impactMax,
            interactionPointCut=seedFinderConfigArg.interactionPointCut,
            deltaZMax=seedFinderConfigArg.deltaZMax,
            maxPtScattering=seedFinderConfigArg.maxPtScattering,
            zBinEdges=seedFinderConfigArg.zBinEdges,
            zBinsCustomLooping=seedFinderConfigArg.zBinsCustomLooping,
            rRangeMiddleSP=seedFinderConfigArg.rRangeMiddleSP,
            useVariableMiddleSPRange=seedFinderConfigArg.useVariableMiddleSPRange,
            binSizeR=seedFinderConfigArg.binSizeR,
            seedConfirmation=seedFinderConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfigArg.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfigArg.forwardSeedConfirmationRange,
        ),
    )
    seedFinderOptions = acts.SeedFinderOptions(
        **acts.examples.defaultKWArgs(
            beamPos=(
                acts.Vector2(0.0, 0.0)
                if seedFinderOptionsArg.beamPos == (None, None)
                else acts.Vector2(
                    seedFinderOptionsArg.beamPos[0], seedFinderOptionsArg.beamPos[1]
                )
            ),
            bFieldInZ=seedFinderOptionsArg.bFieldInZ,
        )
    )
    seedFilterConfig = acts.SeedFilterConfig(
        **acts.examples.defaultKWArgs(
            maxSeedsPerSpM=seedFinderConfig.maxSeedsPerSpM,
            deltaRMin=(
                seedFinderConfig.deltaRMin
                if seedFilterConfigArg.deltaRMin is None
                else seedFilterConfigArg.deltaRMin
            ),
            impactWeightFactor=seedFilterConfigArg.impactWeightFactor,
            zOriginWeightFactor=seedFilterConfigArg.zOriginWeightFactor,
            compatSeedWeight=seedFilterConfigArg.compatSeedWeight,
            compatSeedLimit=seedFilterConfigArg.compatSeedLimit,
            numSeedIncrement=seedFilterConfigArg.numSeedIncrement,
            seedWeightIncrement=seedFilterConfigArg.seedWeightIncrement,
            seedConfirmation=seedFilterConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfig.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfig.forwardSeedConfirmationRange,
            maxSeedsPerSpMConf=seedFilterConfigArg.maxSeedsPerSpMConf,
            maxQualitySeedsPerSpMConf=seedFilterConfigArg.maxQualitySeedsPerSpMConf,
            useDeltaRorTopRadius=seedFilterConfigArg.useDeltaRorTopRadius,
        )
    )

    gridConfig = acts.SpacePointGridConfig(
        **acts.examples.defaultKWArgs(
            minPt=seedFinderConfig.minPt,
            rMax=(
                seedFinderConfig.rMax
                if spacePointGridConfigArg.rMax is None
                else spacePointGridConfigArg.rMax
            ),
            zMax=seedFinderConfig.zMax,
            zMin=seedFinderConfig.zMin,
            deltaRMax=(
                seedFinderConfig.deltaRMax
                if spacePointGridConfigArg.deltaRMax is None
                else spacePointGridConfigArg.deltaRMax
            ),
            cotThetaMax=seedFinderConfig.cotThetaMax,
            phiMin=spacePointGridConfigArg.phi[0],
            phiMax=spacePointGridConfigArg.phi[1],
            maxPhiBins=spacePointGridConfigArg.maxPhiBins,
            impactMax=spacePointGridConfigArg.impactMax,
            zBinEdges=spacePointGridConfigArg.zBinEdges,
            phiBinDeflectionCoverage=spacePointGridConfigArg.phiBinDeflectionCoverage,
        )
    )

    gridOptions = acts.SpacePointGridOptions(
        **acts.examples.defaultKWArgs(
            bFieldInZ=seedFinderOptions.bFieldInZ,
        )
    )

    # Hashing configuration
    hashingTrainingConfig = acts.hashing.HashingTrainingConfig(
        **acts.examples.defaultKWArgs(
            annoySeed=hashingTrainingConfigArg.annoySeed,
            f=hashingTrainingConfigArg.f,
        ),
    )

    hashingConfig = acts.hashing.HashingAlgorithmConfig(
        **acts.examples.defaultKWArgs(
            bucketSize=hashingAlgorithmConfigArg.bucketSize,
            zBins=hashingAlgorithmConfigArg.zBins,
            phiBins=hashingAlgorithmConfigArg.phiBins,
        ),
    )

    # Seeding algorithm
    seedingAlg = SeedingAlgorithmHashing(
        level=logLevel,
        inputSpacePoints=[spacePoints],
        outputSeeds="seeds",
        outputBuckets="buckets",
        **acts.examples.defaultKWArgs(
            allowSeparateRMax=seedingAlgorithmConfigArg.allowSeparateRMax,
            zBinNeighborsTop=seedingAlgorithmConfigArg.zBinNeighborsTop,
            zBinNeighborsBottom=seedingAlgorithmConfigArg.zBinNeighborsBottom,
            numPhiNeighbors=seedingAlgorithmConfigArg.numPhiNeighbors,
        ),
        gridConfig=gridConfig,
        gridOptions=gridOptions,
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
        seedFinderOptions=seedFinderOptions,
        hashingConfig=hashingConfig,
        hashingTrainingConfig=hashingTrainingConfig,
    )
    sequence.addAlgorithm(seedingAlg)

    return seedingAlg.config.outputSeeds, seedingAlg.config.outputBuckets


def addHoughTransformSeeding(
    sequence: acts.examples.Sequencer,
    config: acts.examples.HoughTransformSeeder.Config,
    logLevel: acts.logging.Level = None,
):
    """
    Configures HoughTransform (HT) for seeding, instead of extra proxy config objects it takes
    directly the HT example algorithm config.
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    ht = acts.examples.HoughTransformSeeder(config=config, level=logLevel)
    sequence.addAlgorithm(ht)
    # potentially HT can be extended to also produce seeds, but it is not implemented yet
    # configuration option (outputSeeds) exists
    return ht.config.outputSeeds


def addGbtsSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    seedFinderConfigArg: SeedFinderConfigArg,
    seedFinderOptionsArg: SeedFinderOptionsArg,
    seedFilterConfigArg: SeedFilterConfigArg,
    trackingGeometry: acts.TrackingGeometry,
    logLevel: acts.logging.Level = None,
    layerMappingConfigFile: Union[Path, str] = None,
    geoSelectionConfigFile: Union[Path, str] = None,
    connector_inputConfigFile: Union[Path, str] = None,
):
    """Gbts seeding"""

    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    layerMappingFile = str(layerMappingConfigFile)  # turn path into string
    connector_inputFile = str(connector_inputConfigFile)
    seedFinderConfig = acts.SeedFinderGbtsConfig(
        **acts.examples.defaultKWArgs(
            sigmaScattering=seedFinderConfigArg.sigmaScattering,
            maxSeedsPerSpM=seedFinderConfigArg.maxSeedsPerSpM,
            minPt=seedFinderConfigArg.minPt,
            connector_input_file=connector_inputFile,
            m_useClusterWidth=False,
        ),
    )
    seedFinderOptions = acts.SeedFinderOptions(
        **acts.examples.defaultKWArgs(
            beamPos=(
                acts.Vector2(0.0, 0.0)
                if seedFinderOptionsArg.beamPos == (None, None)
                else acts.Vector2(
                    seedFinderOptionsArg.beamPos[0], seedFinderOptionsArg.beamPos[1]
                )
            ),
            bFieldInZ=seedFinderOptionsArg.bFieldInZ,
        )
    )
    seedFilterConfig = acts.SeedFilterConfig(
        **acts.examples.defaultKWArgs(
            maxSeedsPerSpM=seedFinderConfig.maxSeedsPerSpM,
            deltaRMin=(
                seedFinderConfigArg.deltaR[0]
                if seedFilterConfigArg.deltaRMin is None
                else seedFilterConfigArg.deltaRMin
            ),
            impactWeightFactor=seedFilterConfigArg.impactWeightFactor,
            zOriginWeightFactor=seedFilterConfigArg.zOriginWeightFactor,
            compatSeedWeight=seedFilterConfigArg.compatSeedWeight,
            compatSeedLimit=seedFilterConfigArg.compatSeedLimit,
            numSeedIncrement=seedFilterConfigArg.numSeedIncrement,
            seedWeightIncrement=seedFilterConfigArg.seedWeightIncrement,
            seedConfirmation=seedFilterConfigArg.seedConfirmation,
            # curvatureSortingInFilter=seedFilterConfigArg.curvatureSortingInFilter,
            maxSeedsPerSpMConf=seedFilterConfigArg.maxSeedsPerSpMConf,
            maxQualitySeedsPerSpMConf=seedFilterConfigArg.maxQualitySeedsPerSpMConf,
            useDeltaRorTopRadius=seedFilterConfigArg.useDeltaRorTopRadius,
        )
    )

    seedingAlg = acts.examples.GbtsSeedingAlgorithm(
        level=logLevel,
        inputSpacePoints=[spacePoints],
        outputSeeds="seeds",
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
        seedFinderOptions=seedFinderOptions,
        layerMappingFile=layerMappingFile,
        geometrySelection=acts.examples.readJsonGeometryList(
            str(geoSelectionConfigFile)
        ),
        trackingGeometry=trackingGeometry,
        fill_module_csv=False,
        inputClusters="clusters",
    )

    sequence.addAlgorithm(seedingAlg)
    return seedingAlg.config.outputSeeds


def addSeedPerformanceWriters(
    sequence: acts.examples.Sequencer,
    outputDirRoot: Union[Path, str],
    seeds: str,
    prototracks: str,
    selectedParticles: str,
    inputParticles: str,
    outputTrackParameters: str,
    logLevel: acts.logging.Level = None,
):
    """Writes seeding related performance output"""
    customLogLevel = acts.examples.defaultLogging(sequence, logLevel)
    outputDirRoot = Path(outputDirRoot)
    if not outputDirRoot.exists():
        outputDirRoot.mkdir()

    sequence.addWriter(
        acts.examples.SeedingPerformanceWriter(
            level=customLogLevel(minLevel=acts.logging.DEBUG),
            inputSeeds=seeds,
            inputParticles=selectedParticles,
            inputMeasurementParticlesMap="measurement_particles_map",
            filePath=str(outputDirRoot / "performance_seeding.root"),
        )
    )

    sequence.addWriter(
        acts.examples.RootTrackParameterWriter(
            level=customLogLevel(),
            inputTrackParameters=outputTrackParameters,
            inputProtoTracks=prototracks,
            inputParticles=inputParticles,
            inputSimHits="simhits",
            inputMeasurementParticlesMap="measurement_particles_map",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(outputDirRoot / "estimatedparams.root"),
            treeName="estimatedparams",
        )
    )


acts.examples.NamedTypeArgs(
    config=SeedFilterMLDBScanConfig,
)


def addSeedFilterML(
    s,
    config: SeedFilterMLDBScanConfig = SeedFilterMLDBScanConfig(),
    onnxModelFile: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
) -> None:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)()
    from acts.examples.onnx import SeedFilterMLAlgorithm

    inputParticles = "particles"
    selectedParticles = "particles_selected"
    seeds = "seeds"
    estParams = "estimatedparameters"

    filterML = SeedFilterMLAlgorithm(
        level=customLogLevel,
        inputTrackParameters="estimatedparameters",
        inputSimSeeds="seeds",
        inputSeedFilterNN=onnxModelFile,
        outputTrackParameters="filtered-parameters",
        outputSimSeeds="filtered-seeds",
        **acts.examples.defaultKWArgs(
            epsilonDBScan=config.epsilonDBScan,
            minPointsDBScan=config.minPointsDBScan,
            minSeedScore=config.minSeedScore,
        ),
    )
    s.addAlgorithm(filterML)
    s.addWhiteboardAlias(seeds, "filtered-seeds")
    s.addWhiteboardAlias("estimatedparameters", "filtered-parameters")

    prototracks = "seed-prototracks-ML"
    s.addAlgorithm(
        acts.examples.SeedsToPrototracks(
            level=customLogLevel,
            inputSeeds=seeds,
            outputProtoTracks=prototracks,
        )
    )

    if outputDirRoot is not None:
        addSeedPerformanceWriters(
            s,
            outputDirRoot,
            seeds,
            prototracks,
            selectedParticles,
            inputParticles,
            estParams,
            customLogLevel,
        )

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)

        if not outputDirCsv.exists():
            outputDirCsv.mkdir()

        csvSeedWriter = acts.examples.CsvSeedWriter(
            level=customLogLevel,
            inputTrackParameters=estParams,
            inputSimSeeds=seeds,
            inputSimHits="simhits",
            inputMeasurementParticlesMap="measurement_particles_map",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            outputDir=str(outputDirCsv),
            fileName=str(f"seed.csv"),
        )
        s.addWriter(csvSeedWriter)

    return s


def addKalmanTracks(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    reverseFilteringMomThreshold: float = 0 * u.GeV,
    inputProtoTracks: str = "truth_particle_tracks",
    multipleScattering: bool = True,
    energyLoss: bool = True,
    clusters: str = None,
    calibrator: acts.examples.MeasurementCalibrator = acts.examples.makePassThroughCalibrator(),
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    kalmanOptions = {
        "multipleScattering": multipleScattering,
        "energyLoss": energyLoss,
        "reverseFilteringMomThreshold": reverseFilteringMomThreshold,
        "freeToBoundCorrection": acts.examples.FreeToBoundCorrection(False),
        "level": customLogLevel(),
    }

    fitAlg = acts.examples.TrackFittingAlgorithm(
        level=customLogLevel(),
        inputMeasurements="measurements",
        inputProtoTracks=inputProtoTracks,
        inputInitialTrackParameters="estimatedparameters",
        inputClusters=clusters if clusters is not None else "",
        outputTracks="kf_tracks",
        pickTrack=-1,
        fit=acts.examples.makeKalmanFitterFunction(
            trackingGeometry, field, **kalmanOptions
        ),
        calibrator=calibrator,
    )
    s.addAlgorithm(fitAlg)
    s.addWhiteboardAlias("tracks", fitAlg.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=fitAlg.config.outputTracks,
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="kf_track_particle_matching",
        outputParticleTrackMatching="kf_particle_track_matching",
        doubleMatching=True,
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    return s


def addTruthTrackingGsf(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    inputProtoTracks: str = "truth_particle_tracks",
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    # NOTE we specify clampToRange as True to silence warnings in the test about
    # queries to the loss distribution outside the specified range, since no dedicated
    # approximation for the ODD is done yet.
    bha = acts.examples.AtlasBetheHeitlerApprox.makeDefault(clampToRange=True)

    gsfOptions = {
        "betheHeitlerApprox": bha,
        "maxComponents": 12,
        "componentMergeMethod": acts.examples.ComponentMergeMethod.maxWeight,
        "mixtureReductionAlgorithm": acts.examples.MixtureReductionAlgorithm.KLDistance,
        "weightCutoff": 1.0e-4,
        "level": customLogLevel(),
    }

    gsfAlg = acts.examples.TrackFittingAlgorithm(
        level=customLogLevel(),
        inputMeasurements="measurements",
        inputProtoTracks=inputProtoTracks,
        inputInitialTrackParameters="estimatedparameters",
        outputTracks="gsf_tracks",
        pickTrack=-1,
        fit=acts.examples.makeGsfFitterFunction(trackingGeometry, field, **gsfOptions),
        calibrator=acts.examples.makePassThroughCalibrator(),
    )
    s.addAlgorithm(gsfAlg)
    s.addWhiteboardAlias("tracks", gsfAlg.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=gsfAlg.config.outputTracks,
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="gsf_track_particle_matching",
        outputParticleTrackMatching="gsf_particle_track_matching",
        doubleMatching=True,
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    return s


@acts.examples.NamedTypeArgs(
    trackSelectorConfig=TrackSelectorConfig,
    ckfConfig=CkfConfig,
)
def addCKFTracks(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    trackSelectorConfig: Optional[
        Union[TrackSelectorConfig, List[TrackSelectorConfig]]
    ] = None,
    ckfConfig: CkfConfig = CkfConfig(),
    twoWay: bool = True,
    reverseSearch: bool = False,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeTrackSummary: bool = True,
    writeTrackStates: bool = False,
    writePerformance: bool = True,
    writeCovMat=False,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """This function steers the seeding

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addSeeding)
    trackingGeometry : tracking geometry
    field : magnetic field
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    trackSelectorConfig : TrackSelectorConfig(loc0, loc1, time, eta, absEta, pt, phi, minMeasurements)
        TrackSelector configuration. Each range is specified as a tuple of (min,max).
        Specify as a list(TrackSelectorConfig) for eta-dependent cuts, with binning specified by absEta[1].
        Defaults of no cuts specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TrackSelector.hpp
    writeTrackSummary : bool, True
        write tracksummary_ckf.root ntuple?
    writeTrackStates : bool, False
        write trackstates_ckf.root ntuple? This can be quite large.
    writePerformance : bool, True
        write performance_fitting_ckf.root and performance_finding_ckf.root ntuples?
    writeCovMat : bool, False
        write covaraiance matrices to tracksummary_ckf.root ntuple?
    """

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    tslist = (
        []
        if trackSelectorConfig is None
        else (
            [trackSelectorConfig]
            if type(trackSelectorConfig) is TrackSelectorConfig
            else trackSelectorConfig
        )
    )

    overwriteArgs = dict() if len(tslist) == 1 else dict(absEtaMax=None)
    cutSets = [
        acts.TrackSelector.Config(**(trackSelectorDefaultKWArgs(c) | overwriteArgs))
        for c in tslist
    ]
    if len(tslist) == 0:
        trkSelCfg = None
    elif len(tslist) == 1:
        trkSelCfg = cutSets[0]
    else:
        trkSelCfg = acts.TrackSelector.EtaBinnedConfig(
            cutSets=cutSets,
            absEtaEdges=[cutSets[0].absEtaMin] + [c.absEta[1] for c in tslist],
        )

    # Setup the track finding algorithm with CKF
    # It takes all the source links created from truth hit smearing, seeds from
    # truth particle smearing and source link selection config
    trackFinder = acts.examples.TrackFindingAlgorithm(
        level=customLogLevel(),
        measurementSelectorCfg=acts.MeasurementSelector.Config(
            [
                (
                    acts.GeometryIdentifier(),
                    (
                        [],
                        [ckfConfig.chi2CutOffMeasurement],
                        [ckfConfig.chi2CutOffOutlier],
                        [ckfConfig.numMeasurementsCutOff],
                    ),
                )
            ]
        ),
        inputMeasurements="measurements",
        inputInitialTrackParameters="estimatedparameters",
        inputSeeds=(
            "estimatedseeds"
            if ckfConfig.seedDeduplication or ckfConfig.stayOnSeed
            else ""
        ),
        outputTracks="ckf_tracks",
        findTracks=acts.examples.TrackFindingAlgorithm.makeTrackFinderFunction(
            trackingGeometry, field, customLogLevel()
        ),
        **acts.examples.defaultKWArgs(
            trackingGeometry=trackingGeometry,
            magneticField=field,
            trackSelectorCfg=trkSelCfg,
            maxSteps=ckfConfig.maxSteps,
            twoWay=twoWay,
            reverseSearch=reverseSearch,
            seedDeduplication=ckfConfig.seedDeduplication,
            stayOnSeed=ckfConfig.stayOnSeed,
            pixelVolumeIds=ckfConfig.pixelVolumes,
            stripVolumeIds=ckfConfig.stripVolumes,
            maxPixelHoles=ckfConfig.maxPixelHoles,
            maxStripHoles=ckfConfig.maxStripHoles,
            trimTracks=ckfConfig.trimTracks,
            constrainToVolumeIds=ckfConfig.constrainToVolumes,
            endOfWorldVolumeIds=ckfConfig.endOfWorldVolumes,
        ),
    )
    s.addAlgorithm(trackFinder)
    s.addWhiteboardAlias("tracks", trackFinder.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=trackFinder.config.outputTracks,
        inputParticles="particles_selected",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="ckf_track_particle_matching",
        outputParticleTrackMatching="ckf_particle_track_matching",
        doubleMatching=True,
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    addTrackWriters(
        s,
        name="ckf",
        tracks=trackFinder.config.outputTracks,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeSummary=writeTrackSummary,
        writeStates=writeTrackStates,
        writeFitterPerformance=writePerformance,
        writeFinderPerformance=writePerformance,
        writeCovMat=writeCovMat,
        logLevel=logLevel,
    )

    return s


def addGx2fTracks(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    inputProtoTracks: str = "truth_particle_tracks",
    multipleScattering: bool = False,
    energyLoss: bool = False,
    nUpdateMax: int = 5,
    relChi2changeCutOff: float = 1e-7,
    clusters: str = None,
    calibrator: acts.examples.MeasurementCalibrator = acts.examples.makePassThroughCalibrator(),
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    gx2fOptions = {
        "multipleScattering": multipleScattering,
        "energyLoss": energyLoss,
        "freeToBoundCorrection": acts.examples.FreeToBoundCorrection(False),
        "nUpdateMax": nUpdateMax,
        "relChi2changeCutOff": relChi2changeCutOff,
        "level": customLogLevel(),
    }

    fitAlg = acts.examples.TrackFittingAlgorithm(
        level=customLogLevel(),
        inputMeasurements="measurements",
        inputProtoTracks=inputProtoTracks,
        inputInitialTrackParameters="estimatedparameters",
        inputClusters=clusters if clusters is not None else "",
        outputTracks="gx2f_tracks",
        pickTrack=-1,
        fit=acts.examples.makeGlobalChiSquareFitterFunction(
            trackingGeometry, field, **gx2fOptions
        ),
        calibrator=calibrator,
    )
    s.addAlgorithm(fitAlg)
    s.addWhiteboardAlias("tracks", fitAlg.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=fitAlg.config.outputTracks,
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="gx2f_track_particle_matching",
        outputParticleTrackMatching="gx2f_particle_track_matching",
        doubleMatching=True,
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    return s


def addTrackWriters(
    s: acts.examples.Sequencer,
    name: str,
    tracks: str = "tracks",
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeSummary: bool = True,
    writeStates: bool = False,
    writeFitterPerformance: bool = False,
    writeFinderPerformance: bool = False,
    logLevel: Optional[acts.logging.Level] = None,
    writeCovMat=False,
):
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()

        if writeSummary:
            trackSummaryWriter = acts.examples.RootTrackSummaryWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching="track_particle_matching",
                filePath=str(outputDirRoot / f"tracksummary_{name}.root"),
                treeName="tracksummary",
                writeCovMat=writeCovMat,
            )
            s.addWriter(trackSummaryWriter)

        if writeStates:
            trackStatesWriter = acts.examples.RootTrackStatesWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching="track_particle_matching",
                inputSimHits="simhits",
                inputMeasurementSimHitsMap="measurement_simhits_map",
                filePath=str(outputDirRoot / f"trackstates_{name}.root"),
                treeName="trackstates",
            )
            s.addWriter(trackStatesWriter)

        if writeFitterPerformance:
            trackFitterPerformanceWriter = acts.examples.TrackFitterPerformanceWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching="track_particle_matching",
                filePath=str(outputDirRoot / f"performance_fitting_{name}.root"),
            )
            s.addWriter(trackFitterPerformanceWriter)

        if writeFinderPerformance:
            trackFinderPerfWriter = acts.examples.TrackFinderPerformanceWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching="track_particle_matching",
                inputParticleTrackMatching="particle_track_matching",
                filePath=str(outputDirRoot / f"performance_finding_{name}.root"),
            )
            s.addWriter(trackFinderPerfWriter)

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()

        if writeSummary:
            csvWriter = acts.examples.CsvTrackWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputMeasurementParticlesMap="measurement_particles_map",
                outputDir=str(outputDirCsv),
                fileName=str(f"tracks_{name}.csv"),
            )
            s.addWriter(csvWriter)


@acts.examples.NamedTypeArgs(
    trackSelectorConfig=TrackSelectorConfig,
)
def addTrackSelection(
    s: acts.examples.Sequencer,
    trackSelectorConfig: TrackSelectorConfig,
    inputTracks: str,
    outputTracks: str,
    logLevel: Optional[acts.logging.Level] = None,
) -> acts.examples.TrackSelectorAlgorithm:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    # single cut config for implicit single bin eta configuration
    selectorConfig = acts.TrackSelector.Config(
        **trackSelectorDefaultKWArgs(trackSelectorConfig)
    )

    trackSelector = acts.examples.TrackSelectorAlgorithm(
        level=customLogLevel(),
        inputTracks=inputTracks,
        outputTracks=outputTracks,
        selectorConfig=selectorConfig,
    )

    s.addAlgorithm(trackSelector)

    return trackSelector


ExaTrkXBackend = Enum("ExaTrkXBackend", "Torch Onnx")


def addExaTrkX(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    geometrySelection: Union[Path, str],
    modelDir: Union[Path, str],
    outputDirRoot: Optional[Union[Path, str]] = None,
    backend: Optional[ExaTrkXBackend] = ExaTrkXBackend.Torch,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    # Create space points
    s.addAlgorithm(
        acts.examples.SpacePointMaker(
            level=customLogLevel(),
            inputMeasurements="measurements",
            outputSpacePoints="spacepoints",
            trackingGeometry=trackingGeometry,
            geometrySelection=acts.examples.readJsonGeometryList(
                str(geometrySelection)
            ),
        )
    )

    metricLearningConfig = {
        "level": customLogLevel(),
        "embeddingDim": 8,
        "rVal": 1.6,
        "knnVal": 100,
    }

    filterConfig = {
        "level": customLogLevel(),
        "cut": 0.01,
    }

    gnnConfig = {
        "level": customLogLevel(),
        "cut": 0.5,
    }

    if backend == ExaTrkXBackend.Torch:
        metricLearningConfig["modelPath"] = str(modelDir / "embed.pt")
        metricLearningConfig["selectedFeatures"] = [0, 1, 2]
        filterConfig["modelPath"] = str(modelDir / "filter.pt")
        filterConfig["selectedFeatures"] = [0, 1, 2]
        gnnConfig["modelPath"] = str(modelDir / "gnn.pt")
        gnnConfig["undirected"] = True
        gnnConfig["selectedFeatures"] = [0, 1, 2]

        graphConstructor = acts.examples.TorchMetricLearning(**metricLearningConfig)
        edgeClassifiers = [
            acts.examples.TorchEdgeClassifier(**filterConfig),
            acts.examples.TorchEdgeClassifier(**gnnConfig),
        ]
        trackBuilder = acts.examples.BoostTrackBuilding(customLogLevel())
    elif backend == ExaTrkXBackend.Onnx:
        metricLearningConfig["modelPath"] = str(modelDir / "embedding.onnx")
        metricLearningConfig["spacepointFeatures"] = 3
        filterConfig["modelPath"] = str(modelDir / "filtering.onnx")
        gnnConfig["modelPath"] = str(modelDir / "gnn.onnx")

        graphConstructor = acts.examples.OnnxMetricLearning(**metricLearningConfig)
        edgeClassifiers = [
            acts.examples.OnnxEdgeClassifier(**filterConfig),
            acts.examples.OnnxEdgeClassifier(**gnnConfig),
        ]
        trackBuilder = acts.examples.BoostTrackBuilding(customLogLevel())

    findingAlg = acts.examples.TrackFindingAlgorithmExaTrkX(
        level=customLogLevel(),
        inputSpacePoints="spacepoints",
        outputProtoTracks="exatrkx_prototracks",
        graphConstructor=graphConstructor,
        edgeClassifiers=edgeClassifiers,
        trackBuilder=trackBuilder,
    )
    s.addAlgorithm(findingAlg)
    s.addWhiteboardAlias("prototracks", findingAlg.config.outputProtoTracks)

    s.addAlgorithm(
        acts.examples.PrototracksToTracks(
            level=customLogLevel(),
            inputProtoTracks="prototracks",
            inputMeasurements="measurements",
            outputTracks="tracks",
        )
    )

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks="tracks",
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="exatrkx_track_particle_matching",
        outputParticleTrackMatching="exatrkx_particle_track_matching",
        doubleMatching=True,
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    if outputDirRoot is not None:
        s.addWriter(
            acts.examples.TrackFinderNTupleWriter(
                level=customLogLevel(),
                inputTracks="tracks",
                inputParticles="particles",
                inputParticleMeasurementsMap="particle_measurements_map",
                inputTrackParticleMatching=matchAlg.config.outputTrackParticleMatching,
                filePath=str(Path(outputDirRoot) / "performance_track_finding.root"),
            )
        )

    return s


@acts.examples.NamedTypeArgs(
    config=AmbiguityResolutionConfig,
)
def addAmbiguityResolution(
    s,
    config: AmbiguityResolutionConfig = AmbiguityResolutionConfig(),
    tracks: str = "tracks",
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeTrackSummary: bool = True,
    writeTrackStates: bool = False,
    writePerformance: bool = True,
    writeCovMat=False,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    from acts.examples import GreedyAmbiguityResolutionAlgorithm

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    alg = GreedyAmbiguityResolutionAlgorithm(
        level=customLogLevel(),
        inputTracks=tracks,
        outputTracks="ambi_tracks",
        **acts.examples.defaultKWArgs(
            maximumSharedHits=config.maximumSharedHits,
            nMeasurementsMin=config.nMeasurementsMin,
            maximumIterations=config.maximumIterations,
        ),
    )
    s.addAlgorithm(alg)
    s.addWhiteboardAlias("tracks", alg.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=alg.config.outputTracks,
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="ambi_track_particle_matching",
        outputParticleTrackMatching="ambi_particle_track_matching",
        doubleMatching=True,
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    addTrackWriters(
        s,
        name="ambi",
        tracks=alg.config.outputTracks,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeSummary=writeTrackSummary,
        writeStates=writeTrackStates,
        writeFitterPerformance=writePerformance,
        writeFinderPerformance=writePerformance,
        writeCovMat=writeCovMat,
        logLevel=logLevel,
    )

    return s


@acts.examples.NamedTypeArgs(
    config=ScoreBasedAmbiguityResolutionConfig,
)
def addScoreBasedAmbiguityResolution(
    s,
    config: ScoreBasedAmbiguityResolutionConfig = ScoreBasedAmbiguityResolutionConfig(),
    tracks: str = "tracks",
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    ambiVolumeFile: Optional[Union[Path, str]] = None,
    writeTrackSummary: bool = True,
    writeTrackStates: bool = False,
    writePerformance: bool = True,
    writeCovMat=False,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    from acts.examples import ScoreBasedAmbiguityResolutionAlgorithm

    customLogLevel = acts.examples.defaultLogging(s, acts.logging.INFO)

    algScoreBased = ScoreBasedAmbiguityResolutionAlgorithm(
        level=customLogLevel(),
        inputTracks=tracks,
        configFile=ambiVolumeFile,
        outputTracks="ambiTracksScoreBased",
        **acts.examples.defaultKWArgs(
            minScore=config.minScore,
            minScoreSharedTracks=config.minScoreSharedTracks,
            maxShared=config.maxShared,
            maxSharedTracksPerMeasurement=config.maxSharedTracksPerMeasurement,
            phiMax=config.phiMax,
            phiMin=config.phiMin,
            etaMax=config.etaMax,
            etaMin=config.etaMin,
            useAmbiguityFunction=config.useAmbiguityFunction,
        ),
    )
    s.addAlgorithm(algScoreBased)
    s.addWhiteboardAlias("tracks", algScoreBased.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=algScoreBased.config.outputTracks,
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="ambi_scorebased_track_particle_matching",
        outputParticleTrackMatching="ambi_scorebased_particle_track_matching",
        doubleMatching=True,
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    addTrackWriters(
        s,
        name="ambi_scorebased",
        tracks=algScoreBased.config.outputTracks,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeSummary=writeTrackSummary,
        writeStates=writeTrackStates,
        writeFitterPerformance=writePerformance,
        writeFinderPerformance=writePerformance,
        writeCovMat=writeCovMat,
        logLevel=logLevel,
    )

    return s


@acts.examples.NamedTypeArgs(
    config=AmbiguityResolutionMLConfig,
)
def addAmbiguityResolutionML(
    s,
    config: AmbiguityResolutionMLConfig = AmbiguityResolutionMLConfig(),
    tracks: str = "tracks",
    onnxModelFile: Optional[Union[Path, str]] = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeTrackSummary: bool = True,
    writeTrackStates: bool = False,
    writePerformance: bool = True,
    writeCovMat=False,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    from acts.examples.onnx import AmbiguityResolutionMLAlgorithm
    from acts.examples import GreedyAmbiguityResolutionAlgorithm

    customLogLevel = acts.examples.defaultLogging(s, logLevel)
    algML = AmbiguityResolutionMLAlgorithm(
        level=customLogLevel(),
        inputTracks=tracks,
        inputDuplicateNN=onnxModelFile,
        outputTracks="ambiTracksML",
        **acts.examples.defaultKWArgs(
            nMeasurementsMin=config.nMeasurementsMin,
        ),
    )

    algGreedy = GreedyAmbiguityResolutionAlgorithm(
        level=customLogLevel(),
        inputTracks=algML.config.outputTracks,
        outputTracks="ambiTracksMLGreedy",
        **acts.examples.defaultKWArgs(
            maximumSharedHits=config.maximumSharedHits,
            nMeasurementsMin=config.nMeasurementsMin,
            maximumIterations=config.maximumIterations,
        ),
    )

    s.addAlgorithm(algML)
    s.addAlgorithm(algGreedy)

    s.addWhiteboardAlias("tracks", algGreedy.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=algGreedy.config.outputTracks,
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="ambiML_track_particle_matching",
        outputParticleTrackMatching="ambiML_particle_track_matching",
        doubleMatching=True,
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    addTrackWriters(
        s,
        name="ambiML",
        tracks=algGreedy.config.outputTracks,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeSummary=writeTrackSummary,
        writeStates=writeTrackStates,
        writeFitterPerformance=writePerformance,
        writeFinderPerformance=writePerformance,
        writeCovMat=writeCovMat,
        logLevel=logLevel,
    )

    return s


@acts.examples.NamedTypeArgs(
    trackSelectorConfig=TrackSelectorConfig,
)
def addVertexFitting(
    s,
    field,
    tracks: Optional[str] = "tracks",
    trackParameters: Optional[str] = None,
    outputProtoVertices: str = "protovertices",
    outputVertices: str = "fittedVertices",
    vertexFinder: VertexFinder = VertexFinder.Truth,
    maxIterations: Optional[int] = None,
    initialVariances: Optional[List[float]] = None,
    useTime: Optional[bool] = False,
    seeder: Optional[acts.VertexSeedFinder] = acts.VertexSeedFinder.GaussianSeeder,
    spatialBinExtent: Optional[float] = None,
    temporalBinExtent: Optional[float] = None,
    trackSelectorConfig: Optional[TrackSelectorConfig] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """This function steers the vertex fitting

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from
        addVertexFitting)
    field : magnetic field
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    vertexFinder : VertexFinder, Truth
        vertexFinder algorithm: one of Truth, AMVF, Iterative
    seeder : enum member
        determines vertex seeder for AMVF, can be acts.seeder.GaussianSeeder or
        acts.seeder.AdaptiveGridSeeder
    useTime : bool, False
        determines whether time information is used in vertex seeder, finder,
        and fitter
        only implemented for the AMVF and the AdaptiveGridSeeder
    spatialBinExtent : float, None
        spatial bin extent for the AdaptiveGridSeeder
    temporalBinExtent : float, None
        temporal bin extent for the AdaptiveGridSeeder
    logLevel : acts.logging.Level, None
        logging level to override setting given in `s`
    """
    from acts.examples import (
        TruthVertexFinder,
        VertexFitterAlgorithm,
        IterativeVertexFinderAlgorithm,
        AdaptiveMultiVertexFinderAlgorithm,
        VertexNTupleWriter,
    )

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    if tracks is not None and trackSelectorConfig is not None:
        trackSelector = addTrackSelection(
            s,
            trackSelectorConfig,
            inputTracks=tracks,
            outputTracks="selectedTracksVertexing",
            logLevel=customLogLevel(),
        )
        tracks = trackSelector.config.outputTracks

    if trackParameters is None:
        converter = acts.examples.TracksToParameters(
            level=customLogLevel(),
            inputTracks=tracks,
            outputTrackParameters="selectedTracksParametersVertexing",
        )
        s.addAlgorithm(converter)
        trackParameters = converter.config.outputTrackParameters

    tracks = tracks if tracks is not None else ""
    inputParticles = "particles"
    selectedParticles = "particles_selected"
    inputTruthVertices = "vertices_truth"

    if vertexFinder == VertexFinder.Truth:
        findVertices = TruthVertexFinder(
            level=customLogLevel(),
            inputTracks=tracks,
            inputParticles=selectedParticles,
            inputMeasurementParticlesMap="measurement_particles_map",
            outputProtoVertices=outputProtoVertices,
            excludeSecondaries=True,
        )
        s.addAlgorithm(findVertices)
        fitVertices = VertexFitterAlgorithm(
            level=customLogLevel(),
            inputTrackParameters=trackParameters,
            inputProtoVertices=findVertices.config.outputProtoVertices,
            outputVertices=outputVertices,
            bField=field,
        )
        s.addAlgorithm(fitVertices)
    elif vertexFinder == VertexFinder.Iterative:
        findVertices = IterativeVertexFinderAlgorithm(
            level=customLogLevel(),
            inputTrackParameters=trackParameters,
            outputProtoVertices=outputProtoVertices,
            outputVertices=outputVertices,
            bField=field,
        )
        s.addAlgorithm(findVertices)
    elif vertexFinder == VertexFinder.AMVF:
        findVertices = AdaptiveMultiVertexFinderAlgorithm(
            level=customLogLevel(),
            inputTrackParameters=trackParameters,
            inputTruthParticles=selectedParticles,
            inputTruthVertices=inputTruthVertices,
            outputProtoVertices=outputProtoVertices,
            outputVertices=outputVertices,
            bField=field,
            seedFinder=seeder,
            **acts.examples.defaultKWArgs(
                maxIterations=maxIterations,
                initialVariances=initialVariances,
                useTime=useTime,
                spatialBinExtent=spatialBinExtent,
                temporalBinExtent=temporalBinExtent,
            ),
        )
        s.addAlgorithm(findVertices)
    else:
        raise RuntimeError("Invalid finder argument")

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()
        s.addWriter(
            VertexNTupleWriter(
                level=customLogLevel(),
                inputVertices=outputVertices,
                inputTracks=tracks,
                inputTruthVertices=inputTruthVertices,
                inputParticles=inputParticles,
                inputSelectedParticles=selectedParticles,
                inputTrackParticleMatching="track_particle_matching",
                bField=field,
                treeName="vertexing",
                filePath=str(outputDirRoot / "performance_vertexing.root"),
            )
        )

    return s


def addSingleSeedVertexFinding(
    s,
    outputDirRoot: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
    inputSpacePoints: Optional[str] = "spacepoints",
    outputVertices: Optional[str] = "fittedSeedVertices",
) -> None:
    from acts.examples import (
        SingleSeedVertexFinderAlgorithm,
        VertexNTupleWriter,
    )

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    findSingleSeedVertex = SingleSeedVertexFinderAlgorithm(
        level=customLogLevel(),
        inputSpacepoints=inputSpacePoints,
        outputVertices=outputVertices,
    )
    s.addAlgorithm(findSingleSeedVertex)

    inputParticles = "particles"
    selectedParticles = "particles_selected"

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()

        s.addWriter(
            VertexNTupleWriter(
                level=customLogLevel(),
                inputAllTruthParticles=inputParticles,
                inputSelectedTruthParticles=selectedParticles,
                useTracks=False,
                inputVertices=outputVertices,
                treeName="seedvertexing",
                filePath=str(outputDirRoot / "performance_seedvertexing.root"),
            )
        )

    return s
