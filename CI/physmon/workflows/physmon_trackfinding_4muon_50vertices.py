#!/usr/bin/env python3

import tempfile
from pathlib import Path
import shutil

import acts
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addFatras,
    addDigitization,
    ParticleSelectorConfig,
    addDigiParticleSelection,
)
from acts.examples.reconstruction import (
    addSeeding,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedingAlgorithm,
    CkfConfig,
    addCKFTracks,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
    TrackSelectorConfig,
)

from physmon_common import makeSetup

u = acts.UnitConstants

setup = makeSetup()


with tempfile.TemporaryDirectory() as temp:
    s = acts.examples.Sequencer(
        events=3,
        numThreads=-1,
        logLevel=acts.logging.INFO,
    )

    tp = Path(temp)

    for d in setup.decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)

    addParticleGun(
        s,
        MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
        EtaConfig(-3.0, 3.0, uniform=True),
        PhiConfig(0.0, 360.0 * u.degree),
        ParticleConfig(10, acts.PdgParticle.eMuon, randomizeCharge=True),
        vtxGen=acts.examples.GaussianVertexGenerator(
            mean=acts.Vector4(0, 0, 0, 0),
            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 1.0 * u.ns),
        ),
        multiplicity=50,
        rnd=rnd,
    )

    addFatras(
        s,
        setup.trackingGeometry,
        setup.field,
        rnd=rnd,
    )

    addDigitization(
        s,
        setup.trackingGeometry,
        setup.field,
        digiConfigFile=setup.digiConfig,
        rnd=rnd,
    )

    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(0.9 * u.GeV, None),
            measurements=(9, None),
            removeNeutral=True,
        ),
    )

    addSeeding(
        s,
        setup.trackingGeometry,
        setup.field,
        SeedFinderConfigArg(
            r=(33 * u.mm, 200 * u.mm),
            deltaR=(1 * u.mm, 60 * u.mm),
            collisionRegion=(-250 * u.mm, 250 * u.mm),
            z=(-2000 * u.mm, 2000 * u.mm),
            maxSeedsPerSpM=1,
            sigmaScattering=5,
            radLengthPerSeed=0.1,
            minPt=0.5 * u.GeV,
            impactMax=3 * u.mm,
        ),
        SeedFinderOptionsArg(bFieldInZ=2 * u.T, beamPos=(0.0, 0.0)),
        seedingAlgorithm=SeedingAlgorithm.Default,
        initialSigmas=[
            1 * u.mm,
            1 * u.mm,
            1 * u.degree,
            1 * u.degree,
            0.1 * u.e / u.GeV,
            1 * u.ns,
        ],
        initialSigmaPtRel=0.01,
        initialVarInflation=[1.0] * 6,
        geoSelectionConfigFile=setup.geoSel,
        outputDirRoot=tp,
    )

    addCKFTracks(
        s,
        setup.trackingGeometry,
        setup.field,
        TrackSelectorConfig(
            pt=(0.9 * u.GeV, None),
            loc0=(-4.0 * u.mm, 4.0 * u.mm),
            nMeasurementsMin=6,
            maxHoles=2,
            maxOutliers=2,
        ),
        CkfConfig(
            chi2CutOffMeasurement=15.0,
            chi2CutOffOutlier=25.0,
            numMeasurementsCutOff=10,
            seedDeduplication=True,
            stayOnSeed=True,
        ),
        outputDirRoot=tp,
    )

    addAmbiguityResolution(
        s,
        AmbiguityResolutionConfig(
            maximumSharedHits=3,
            maximumIterations=100000,
            nMeasurementsMin=6,
        ),
        outputDirRoot=tp,
    )

    s.addAlgorithm(
        acts.examples.TracksToParameters(
            level=acts.logging.INFO,
            inputTracks="tracks",
            outputTrackParameters="trackParameters",
        )
    )

    # Choosing a seeder only has an effect on VertexFinder.AMVF. For
    # VertexFinder.IVF we always use acts.VertexSeedFinder.GaussianSeeder
    # (Python binding is not implemented).
    # Setting useTime also only has an effect on VertexFinder.AMVF due to
    # the same reason.
    addVertexFitting(
        s,
        setup.field,
        trackParameters="trackParameters",
        outputProtoVertices="ivf_notime_protovertices",
        outputVertices="ivf_notime_fittedVertices",
        vertexFinder=VertexFinder.Iterative,
        outputDirRoot=tp / "ivf_notime",
    )

    addVertexFitting(
        s,
        setup.field,
        trackParameters="trackParameters",
        outputProtoVertices="amvf_gauss_notime_protovertices",
        outputVertices="amvf_gauss_notime_fittedVertices",
        seeder=acts.VertexSeedFinder.GaussianSeeder,
        useTime=False,  # Time seeding not implemented for the Gaussian seeder
        vertexFinder=VertexFinder.AMVF,
        outputDirRoot=tp / "amvf_gauss_notime",
    )

    addVertexFitting(
        s,
        setup.field,
        trackParameters="trackParameters",
        outputProtoVertices="amvf_grid_time_protovertices",
        outputVertices="amvf_grid_time_fittedVertices",
        seeder=acts.VertexSeedFinder.AdaptiveGridSeeder,
        useTime=True,
        vertexFinder=VertexFinder.AMVF,
        outputDirRoot=tp / "amvf_grid_time",
    )

    s.run()

    shutil.move(
        tp / "performance_finding_ambi.root",
        tp / "performance_finding_ckf_ambi.root",
    )
    shutil.move(
        tp / "performance_fitting_ambi.root",
        tp / "performance_fitting_ckf_ambi.root",
    )
    for vertexing in ["ivf_notime", "amvf_gauss_notime", "amvf_grid_time"]:
        shutil.move(
            tp / f"{vertexing}/performance_vertexing.root",
            tp / f"performance_vertexing_{vertexing}.root",
        )

    for file in [
        "performance_seeding.root",
        "tracksummary_ckf.root",
        "performance_finding_ckf.root",
        "performance_fitting_ckf.root",
        "performance_finding_ckf_ambi.root",
        "performance_fitting_ckf_ambi.root",
        "performance_vertexing_ivf_notime.root",
        "performance_vertexing_amvf_gauss_notime.root",
        "performance_vertexing_amvf_grid_time.root",
    ]:
        perf_file = tp / file
        assert perf_file.exists(), f"Performance file not found {perf_file}"
        shutil.copy(perf_file, setup.outdir / file)
