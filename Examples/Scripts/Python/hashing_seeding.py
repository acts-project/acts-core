#!/usr/bin/env python3

import argparse
import pathlib

import acts
import acts.examples
from acts.examples.simulation import (
    addPythia8,
    ParticleSelectorConfig,
    addGenParticleSelection,
    addFatras,
    addDigitization,
    addDigiParticleSelection,
)

from acts.examples.reconstruction import (
    addSeeding,
    addCKFTracks,
    TrackSelectorConfig,
    SeedingAlgorithm,
)

from pathlib import Path

from typing import Optional
from enum import Enum

DetectorName = Enum("DetectorName", "ODD generic")

HashingMetric = Enum("HashingMetric", "dphi dR")


class Config:
    def __init__(
        self,
        mu: int = None,
        bucketSize: int = 100,
        maxSeedsPerSpM: int = 1000,
        detector: DetectorName = DetectorName.generic,
        seedingAlgorithm: SeedingAlgorithm = SeedingAlgorithm.Hashing,
        metric: str = HashingMetric.dphi,
        annoySeed: int = 123456789,
        zBins: int = 100_000,
        phiBins: int = 0,
    ):
        self.mu = mu
        self.bucketSize = bucketSize
        self.maxSeedsPerSpM = maxSeedsPerSpM
        self.detector = detector
        self.seedingAlgorithm = (
            seedingAlgorithm
            if type(seedingAlgorithm) == SeedingAlgorithm
            else SeedingAlgorithm[seedingAlgorithm]
        )
        self.metric = metric if type(metric) == HashingMetric else HashingMetric[metric]
        self.annoySeed = annoySeed
        self.zBins = zBins
        self.phiBins = phiBins

        if seedingAlgorithm == SeedingAlgorithm.Default:
            self.bucketSize = 0
            self.metric = HashingMetric.dphi
            self.annoySeed = 123456789
            self.zBins = 0
            self.phiBins = 0

    def __str__(self):
        return (
            f"mu={self.mu}, bucketSize={self.bucketSize}, maxSeedsPerSpM={self.maxSeedsPerSpM}, "
            f"detector={self.detector}, seedingAlgorithm={self.seedingAlgorithm}, metric={self.metric}, "
            f"annoySeed={self.annoySeed}, zBins={self.zBins}, phiBins={self.phiBins}"
        )

    def __repr__(self):
        return self.__str__()

    def save(self, path: pathlib.Path):
        with open(path, "w") as f:
            f.write(str(self))

    @property
    def doHashing(self):
        return self.bucketSize > 0

    @property
    def directory(self):
        outDir = f"detector_{extractEnumName(self.detector)}"
        outDir += "_output"
        doHashing = self.bucketSize > 0
        if doHashing:
            outDir += "_hashing"

        outDir += f"_mu_{self.mu}"
        if doHashing:
            outDir += f"_bucket_{self.bucketSize}"

        outDir += f"_maxSeedsPerSpM_{self.maxSeedsPerSpM}"

        outDir += f"_seedFinderConfig_{'TrackML' if self.detector == DetectorName.generic else 'ODD'}"

        outDir += f"_seedingAlgorithm_{extractEnumName(self.seedingAlgorithm)}Seeding"
        if doHashing:
            if self.metric != "angular":
                outDir += f"_metric_{extractEnumName(self.metric)}"
            outDir += f"_annoySeed_{self.annoySeed}"
            if self.zBins != 0:
                outDir += f"_zBins_{self.zBins}"
            if self.phiBins != 0:
                outDir += f"_phiBins_{self.phiBins}"

        return outDir

    def getDetectorInfo(self):
        actsExamplesDir = getActsExamplesDirectory()

        if self.detector == DetectorName.ODD:
            from acts.examples.odd import (
                getOpenDataDetector,
                getOpenDataDetectorDirectory,
            )

            geoDir = getOpenDataDetectorDirectory()

            oddMaterialMap = geoDir / "data/odd-material-maps.root"
            oddDigiConfig = geoDir / "config/odd-digi-smearing-config.json"
            oddSeedingSel = geoDir / "config/odd-seeding-config.json"
            oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

            detector = getOpenDataDetector(odd_dir=geoDir, mdecorator=oddMaterialDeco)
            trackingGeometry = detector.trackingGeometry()

            digiConfig = oddDigiConfig

            geoSelectionConfigFile = oddSeedingSel

        elif self.detector == DetectorName.generic:
            print("Create detector and tracking geometry")

            detector = acts.examples.GenericDetector()
            trackingGeometry = detector.trackingGeometry()
            digiConfig = (
                actsExamplesDir
                / "Algorithms/Digitization/share/default-smearing-config-generic.json"
            )
            geoSelectionConfigFile = (
                actsExamplesDir
                / "Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
            )
        else:
            exit("Detector not supported")

        return detector, trackingGeometry, digiConfig, geoSelectionConfigFile


def extractEnumName(enumvar):
    return str(enumvar).split(".")[-1]


def getActsExamplesDirectory():
    return Path(__file__).parent.parent.parent


def runHashingSeeding(
    nevents: int,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    outputDir: pathlib.Path,
    saveFiles: bool,
    npileup: int,
    seedingAlgorithm: SeedingAlgorithm,
    maxSeedsPerSpM: int,
    digiConfig: pathlib.Path,
    geoSelectionConfigFile: pathlib.Path,
    eta=4,
    rnd: acts.examples.RandomNumbers = acts.examples.RandomNumbers(seed=42),
    config: Config = None,
    s=None,
):
    from acts.examples.reconstruction import (
        SeedFinderConfigArg,
        SeedFinderOptionsArg,
        HashingTrainingConfigArg,
        HashingAlgorithmConfigArg,
    )

    u = acts.UnitConstants

    outputDir = Path(outputDir)

    s = s or acts.examples.Sequencer(
        events=nevents,
        numThreads=1,
        outputDir=str(outputDir),
        trackFpes=False,
    )

    addPythia8(
        s,
        hardProcess=["Top:qqbar2ttbar=on"],
        npileup=npileup,
        vtxGen=acts.examples.GaussianVertexGenerator(
            stddev=acts.Vector4(0, 0, 50 * u.mm, 0),
            mean=acts.Vector4(0, 0, 0, 0),
        ),
        rnd=rnd,
    )

    addGenParticleSelection(
        s,
        ParticleSelectorConfig(
            rho=(0.0, 24 * u.mm),
            absZ=(0.0, 1.0 * u.m),
        ),
    )

    addFatras(
        s,
        trackingGeometry,
        field,
        enableInteractions=True,
        # outputDirRoot=outputDir,  # RootParticle ERROR when setting the outputDirRoot
        outputDirCsv=outputDir if saveFiles else None,
        rnd=rnd,
    )

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfig,
        outputDirRoot=outputDir,
        outputDirCsv=outputDir if saveFiles else None,
        rnd=rnd,
    )

    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(1.0 * u.GeV, None),
            eta=(-eta, eta),
            measurements=(9, None),
            removeNeutral=True,
        ),
    )

    import numpy as np

    cotThetaMax = 1 / (np.tan(2 * np.arctan(np.exp(-eta))))  # =1/tan(2×atan(e^(-eta)))
    seedFinderConfigArg = SeedFinderConfigArg(
        r=(None, 200 * u.mm),  # rMin=default, 33mm
        deltaR=(1 * u.mm, 60 * u.mm),
        collisionRegion=(-250 * u.mm, 250 * u.mm),
        z=(-2000 * u.mm, 2000 * u.mm),
        maxSeedsPerSpM=maxSeedsPerSpM,
        sigmaScattering=5,
        radLengthPerSeed=0.1,
        minPt=500 * u.MeV,
        impactMax=3 * u.mm,
        cotThetaMax=cotThetaMax,
    )

    seedFinderOptionsArg: SeedFinderOptionsArg = SeedFinderOptionsArg(bFieldInZ=2 * u.T)
    hashingTrainingConfigArg: HashingTrainingConfigArg = HashingTrainingConfigArg(
        annoySeed=config.annoySeed, f=2 if config.metric == HashingMetric.dR else 1
    )
    hashingAlgorithmConfigArg: HashingAlgorithmConfigArg = HashingAlgorithmConfigArg(
        bucketSize=config.bucketSize,
        zBins=config.zBins,
        phiBins=config.phiBins,
    )
    initialSigmas: Optional[list] = [
        1 * u.mm,
        1 * u.mm,
        1 * u.degree,
        1 * u.degree,
        0.1 / u.GeV,
        1 * u.ns,
    ]

    initialVarInflation: Optional[list] = [1.0] * 6

    addSeeding(
        s,
        trackingGeometry,
        field,
        seedFinderConfigArg,
        seedFinderOptionsArg,
        hashingTrainingConfigArg,
        hashingAlgorithmConfigArg,
        seedingAlgorithm=seedingAlgorithm,
        geoSelectionConfigFile=geoSelectionConfigFile,
        initialSigmas=initialSigmas,
        initialVarInflation=initialVarInflation,
        outputDirRoot=outputDir,
        outputDirCsv=outputDir if saveFiles else None,
    )

    return s


if __name__ == "__main__":
    eta = 4
    # eta = 2.5

    parser = argparse.ArgumentParser(
        description="Example script to run seed finding with hashing"
    )
    parser.add_argument("--mu", type=int, default=50)
    parser.add_argument("--bucketSize", type=int, default=100)
    parser.add_argument(
        "--nevents",
        type=int,
        default=20,
    )
    parser.add_argument("--maxSeedsPerSpM", type=int, default=1000)
    parser.add_argument("--seedingAlgorithm", type=str, default="Hashing")
    parser.add_argument("--saveFiles", type=bool, default=True)
    parser.add_argument("--annoySeed", type=int, default=123456789)
    parser.add_argument("--zBins", type=int, default=100000)
    parser.add_argument("--phiBins", type=int, default=0)
    parser.add_argument("--metric", type=str, default="dphi")
    args = parser.parse_args()

    print(args)

    mu = args.mu
    bucketSize = args.bucketSize
    nevents = args.nevents
    saveFiles = args.saveFiles
    annoySeed = args.annoySeed
    zBins = args.zBins
    phiBins = args.phiBins
    maxSeedsPerSpM = args.maxSeedsPerSpM

    metric = HashingMetric[args.metric]

    seedingAlgorithm = SeedingAlgorithm[args.seedingAlgorithm]

    detector = DetectorName.generic

    u = acts.UnitConstants

    config = Config(
        mu=mu,
        bucketSize=bucketSize,
        maxSeedsPerSpM=maxSeedsPerSpM,
        detector=detector,
        seedingAlgorithm=seedingAlgorithm,
        metric=metric,
        annoySeed=annoySeed,
        zBins=zBins,
        phiBins=phiBins,
    )

    doHashing = config.doHashing
    bucketSize = config.bucketSize
    npileup = config.mu
    maxSeedsPerSpM = config.maxSeedsPerSpM

    outputDir = pathlib.Path.cwd() / config.directory

    if not outputDir.exists():
        outputDir.mkdir(parents=True)

    config.save(outputDir / "config_file.txt")

    field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
    rnd = acts.examples.RandomNumbers(seed=42)

    _, trackingGeometry, digiConfig, geoSelectionConfigFile = config.getDetectorInfo()

    s = runHashingSeeding(
        nevents,
        trackingGeometry,
        field,
        outputDir,
        saveFiles,
        npileup,
        seedingAlgorithm,
        maxSeedsPerSpM,
        digiConfig,
        geoSelectionConfigFile,
        rnd=rnd,
        eta=eta,
        config=config,
    )

    logLevel = acts.logging.VERBOSE
    rootSpacepointsWriter = acts.examples.RootSpacepointWriter(
        level=logLevel,
        inputSpacepoints="spacepoints",
        filePath=str(outputDir / "spacepoints.root"),
    )
    s.addWriter(rootSpacepointsWriter)

    rootSeedsWriter = acts.examples.RootSeedWriter(
        level=logLevel,
        inputSeeds="seeds",
        filePath=str(outputDir / "seeds.root"),
    )
    s.addWriter(rootSeedsWriter)

    # TrackFinding ERROR no intersection found; TrackExtrapolationError:2
    addCKFTracks(
        s,
        trackingGeometry,
        field,
        TrackSelectorConfig(
            pt=(1.0 * u.GeV, None),
            absEta=(None, eta),
            nMeasurementsMin=6,
        ),
        outputDirRoot=outputDir,
        writeTrajectories=False,
        twoWay=False,
    )

    s.run()
