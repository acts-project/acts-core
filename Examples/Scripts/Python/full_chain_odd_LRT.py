#!/usr/bin/env python3

import os
import argparse
import pathlib
import math

import acts
import acts.examples
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    ParticleSelectorConfig,
    addPythia8,
    ParticleSelectorConfig,
    addGenParticleSelection,
    addFatras,
    addGeant4,
    addSimParticleSelection,
    addDigitization,
    addDigiParticleSelection,
)
from acts.examples.reconstruction import (
    addSeeding,
    CkfConfig,
    addCKFTracks,
    TrackSelectorConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addAmbiguityResolutionML,
    AmbiguityResolutionMLConfig,
    addScoreBasedAmbiguityResolution,
    ScoreBasedAmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
    addSeedFilterML,
    SeedFilterMLDBScanConfig,
)
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

u = acts.UnitConstants


parser = argparse.ArgumentParser(description="Full chain with the OpenDataDetector")
parser.add_argument(
    "--output",
    "-o",
    help="Output directory",
    type=pathlib.Path,
    default=pathlib.Path.cwd() / "odd_output",
)
parser.add_argument("--events", "-n", help="Number of events", type=int, default=100)
parser.add_argument("--skip", "-s", help="Number of events", type=int, default=0)
parser.add_argument("--edm4hep", help="Use edm4hep inputs", type=pathlib.Path)
parser.add_argument(
    "--geant4", help="Use Geant4 instead of fatras", action="store_true"
)
parser.add_argument(
    "--ttbar",
    help="Use Pythia8 (ttbar, pile-up 200) instead of particle gun",
    action="store_true",
)
parser.add_argument(
    "--ttbar-pu",
    help="Number of pile-up events for ttbar",
    type=int,
    default=200,
)
parser.add_argument(
    "--gun-particles",
    help="Multiplicity (no. of particles) of the particle gun",
    type=int,
    default=4,
)
parser.add_argument(
    "--gun-multiplicity",
    help="Multiplicity (no. of vertices) of the particle gun",
    type=int,
    default=200,
)
parser.add_argument(
    "--gun-eta-range",
    nargs=2,
    help="Eta range of the particle gun",
    type=float,
    default=[-3.0, 3.0],
)
parser.add_argument(
    "--gun-pt-range",
    nargs=2,
    help="Pt range of the particle gun (GeV)",
    type=float,
    default=[1.0 * u.GeV, 10.0 * u.GeV],
)
parser.add_argument(
    "--digi-config", help="Digitization configuration file", type=pathlib.Path
)
parser.add_argument(
    "--material-config", help="Material map configuration file", type=pathlib.Path
)
parser.add_argument(
    "--ambi-solver",
    help="Set which ambiguity solver to use, default is the classical one",
    type=str,
    choices=["greedy", "scoring", "ML"],
    default="greedy",
)
parser.add_argument(
    "--ambi-config",
    help="Set the configuration file for the Score Based ambiguity resolution",
    type=pathlib.Path,
    default=pathlib.Path.cwd() / "ambi_config.json",
)

parser.add_argument(
    "--MLSeedFilter",
    help="Use the Ml seed filter to select seed after the seeding step",
    action="store_true",
)
parser.add_argument(
    "--reco",
    help="Switch reco on/off",
    default=True,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--output-root",
    help="Switch root output on/off",
    default=True,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--output-csv",
    help="Switch csv output on/off",
    default=True,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--output-obj",
    help="Switch obj output on/off",
    default=True,
    action=argparse.BooleanOptionalAction,
)

args = parser.parse_args()

outputDir = args.output
ambi_ML = args.ambi_solver == "ML"
ambi_scoring = args.ambi_solver == "scoring"
ambi_config = args.ambi_config
seedFilter_ML = args.MLSeedFilter
geoDir = getOpenDataDetectorDirectory()
# acts.examples.dump_args_calls(locals())  # show python binding calls

oddMaterialMap = (
    args.material_config
    if args.material_config
    else geoDir / "data/odd-material-maps.root"
)

oddDigiConfig = (
    args.digi_config
    if args.digi_config
    else geoDir / "config/odd-digi-smearing-config.json"
)

oddSeedingSel = geoDir / "config/odd-seeding-config.json"
oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

detector = getOpenDataDetector(odd_dir=geoDir, mdecorator=oddMaterialDeco)
trackingGeometry = detector.trackingGeometry()
decorators = detector.contextDecorators()
field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(
    events=args.events,
    skip=args.skip,
    numThreads=1 if args.geant4 else -1,
    outputDir=str(outputDir),
)

if args.edm4hep:
    import acts.examples.edm4hep

    edm4hepReader = acts.examples.edm4hep.EDM4hepReader(
        inputPath=str(args.edm4hep),
        inputSimHits=[
            "PixelBarrelReadout",
            "PixelEndcapReadout",
            "ShortStripBarrelReadout",
            "ShortStripEndcapReadout",
            "LongStripBarrelReadout",
            "LongStripEndcapReadout",
        ],
        outputParticlesGenerator="particles_generated",
        outputParticlesSimulation="particles_simulated",
        outputSimHits="simhits",
        graphvizOutput="graphviz",
        dd4hepDetector=detector,
        trackingGeometry=trackingGeometry,
        sortSimHitsInTime=True,
        level=acts.logging.INFO,
    )
    s.addReader(edm4hepReader)

    s.addWhiteboardAlias("particles", edm4hepReader.config.outputParticlesGenerator)

    addSimParticleSelection(
        s,
        ParticleSelectorConfig(
            rho=(0.0, 24 * u.mm),
            absZ=(0.0, 1.0 * u.m),
            eta=(-3.0, 3.0),
            pt=(150 * u.MeV, None),
            removeNeutral=True,
        ),
    )
else:
    if not args.ttbar:
        addParticleGun(
            s,
            MomentumConfig(
                args.gun_pt_range[0] * u.GeV,
                args.gun_pt_range[1] * u.GeV,
                transverse=True,
            ),
            EtaConfig(args.gun_eta_range[0], args.gun_eta_range[1]),
            PhiConfig(0.0, 360.0 * u.degree),
            ParticleConfig(
                args.gun_particles, acts.PdgParticle.eMuon, randomizeCharge=True
            ),
            vtxGen=acts.examples.GaussianDisplacedVertexPositionGenerator(
                rMean=2,
                rStdDev=0.0125 * u.mm,
                zMean=2,
                zStdDev=55.5 * u.mm,
                tMean=0,
                tStdDev=1.0 * u.ns,
            ),
            multiplicity=args.gun_multiplicity,
            rnd=rnd,
        )
    else:
        addPythia8(
            s,
            hardProcess=["Top:qqbar2ttbar=on"],
            npileup=args.ttbar_pu,
            vtxGen=acts.examples.GaussianDisplacedVertexPositionGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns
                ),
            ),
            rnd=rnd,
            outputDirRoot=outputDir if args.output_root else None,
            outputDirCsv=outputDir if args.output_csv else None,
        )

        addGenParticleSelection(
            s,
            ParticleSelectorConfig(
                rho=(0.0, 24 * u.mm),
                absZ=(0.0, 1.0 * u.m),
                eta=(-3.0, 3.0),
                pt=(150 * u.MeV, None),
            ),
        )

    if args.geant4:
        if s.config.numThreads != 1:
            raise ValueError("Geant 4 simulation does not support multi-threading")

        # Pythia can sometime simulate particles outside the world volume, a cut on the Z of the track help mitigate this effect
        # Older version of G4 might not work, this as has been tested on version `geant4-11-00-patch-03`
        # For more detail see issue #1578
        addGeant4(
            s,
            detector,
            trackingGeometry,
            field,
            outputDirRoot=outputDir if args.output_root else None,
            outputDirCsv=outputDir if args.output_csv else None,
            outputDirObj=outputDir if args.output_obj else None,
            rnd=rnd,
            killVolume=trackingGeometry.highestTrackingVolume,
            killAfterTime=25 * u.ns,
        )
    else:
        addFatras(
            s,
            trackingGeometry,
            field,
            enableInteractions=True,
            outputDirRoot=outputDir if args.output_root else None,
            outputDirCsv=outputDir if args.output_csv else None,
            outputDirObj=outputDir if args.output_obj else None,
            rnd=rnd,
        )

addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=oddDigiConfig,
    outputDirRoot=outputDir if args.output_root else None,
    outputDirCsv=outputDir if args.output_csv else None,
    rnd=rnd,
)

addDigiParticleSelection(
    s,
    ParticleSelectorConfig(
        pt=(1.0 * u.GeV, None),
        eta=(-3.0, 3.0),
        measurements=(9, None),
        removeNeutral=True,
    ),
)

if args.reco:
    addSeeding(
        s,
        trackingGeometry,
        field,
        initialSigmas=[
            1 * u.mm,
            1 * u.mm,
            1 * u.degree,
            1 * u.degree,
            0.1 * u.e / u.GeV,
            1 * u.ns,
        ],
        initialSigmaPtRel=0.1,
        initialVarInflation=[1.0] * 6,
        geoSelectionConfigFile=oddSeedingSel,
        outputDirRoot=outputDir if args.output_root else None,
        outputDirCsv=outputDir if args.output_csv else None,
    )

    if seedFilter_ML:
        addSeedFilterML(
            s,
            SeedFilterMLDBScanConfig(
                epsilonDBScan=0.03, minPointsDBScan=2, minSeedScore=0.1
            ),
            onnxModelFile=os.path.dirname(__file__)
            + "/MLAmbiguityResolution/seedDuplicateClassifier.onnx",
            outputDirRoot=outputDir if args.output_root else None,
            outputDirCsv=outputDir if args.output_csv else None,
        )

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        TrackSelectorConfig(
            pt=(1.0 * u.GeV if args.ttbar else 0.0, None),
            absEta=(None, 3.0),
            loc0=(-4.0 * u.mm, 4.0 * u.mm),
            nMeasurementsMin=7,
            maxHoles=2,
            maxOutliers=2,
        ),
        CkfConfig(
            chi2CutOffMeasurement=15.0,
            chi2CutOffOutlier=25.0,
            numMeasurementsCutOff=10,
            seedDeduplication=True,
            stayOnSeed=True,
            pixelVolumes=[16, 17, 18],
            stripVolumes=[23, 24, 25],
            maxPixelHoles=1,
            maxStripHoles=2,
            constrainToVolumes=[
                2,  # beam pipe
                32,
                4,  # beam pip gap
                16,
                17,
                18,  # pixel
                20,  # PST
                23,
                24,
                25,  # short strip
                26,
                8,  # long strip gap
                28,
                29,
                30,  # long strip
            ],
        ),
        outputDirRoot=outputDir if args.output_root else None,
        outputDirCsv=outputDir if args.output_csv else None,
        writeCovMat=True,
    )

    if ambi_ML:
        addAmbiguityResolutionML(
            s,
            AmbiguityResolutionMLConfig(
                maximumSharedHits=3, maximumIterations=1000000, nMeasurementsMin=7
            ),
            outputDirRoot=outputDir if args.output_root else None,
            outputDirCsv=outputDir if args.output_csv else None,
            onnxModelFile=os.path.dirname(__file__)
            + "/MLAmbiguityResolution/duplicateClassifier.onnx",
        )

    elif ambi_scoring:
        addScoreBasedAmbiguityResolution(
            s,
            ScoreBasedAmbiguityResolutionConfig(
                minScore=0,
                minScoreSharedTracks=1,
                maxShared=2,
                maxSharedTracksPerMeasurement=2,
                pTMax=1400,
                pTMin=0.5,
                phiMax=math.pi,
                phiMin=-math.pi,
                etaMax=4,
                etaMin=-4,
                useAmbiguityFunction=False,
            ),
            outputDirRoot=outputDir if args.output_root else None,
            outputDirCsv=outputDir if args.output_csv else None,
            ambiVolumeFile=ambi_config,
            writeCovMat=True,
        )
    else:
        addAmbiguityResolution(
            s,
            AmbiguityResolutionConfig(
                maximumSharedHits=3, maximumIterations=1000000, nMeasurementsMin=7
            ),
            outputDirRoot=outputDir if args.output_root else None,
            outputDirCsv=outputDir if args.output_csv else None,
            writeCovMat=True,
        )

    addVertexFitting(
        s,
        field,
        vertexFinder=VertexFinder.AMVF,
        outputDirRoot=outputDir if args.output_root else None,
    )

s.run()
