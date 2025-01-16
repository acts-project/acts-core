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
    TrackSmearingSigmas,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedingAlgorithm,
    TruthEstimatedSeedingAlgorithmConfigArg,
    CkfConfig,
    addCKFTracks,
    TrackSelectorConfig,
)

from physmon_common import makeSetup

u = acts.UnitConstants

setup = makeSetup()


def run_ckf_tracking(label, seeding):
    with tempfile.TemporaryDirectory() as temp:
        s = acts.examples.Sequencer(
            events=10000,
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
            ParticleConfig(1, acts.PdgParticle.eMuon, randomizeCharge=True),
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 1.0 * u.ns
                ),
            ),
            multiplicity=1,
            rnd=rnd,
        )

        addFatras(
            s,
            setup.trackingGeometry,
            setup.field,
            enableInteractions=True,
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
            TrackSmearingSigmas(  # only used by SeedingAlgorithm.TruthSmeared
                # zero eveything so the CKF has a chance to find the measurements
                loc0=0,
                loc0PtA=0,
                loc0PtB=0,
                loc1=0,
                loc1PtA=0,
                loc1PtB=0,
                time=0,
                phi=0,
                theta=0,
                ptRel=0,
            ),
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
            SeedFinderOptionsArg(bFieldInZ=2 * u.T),
            TruthEstimatedSeedingAlgorithmConfigArg(deltaR=(10.0 * u.mm, None)),
            seedingAlgorithm=seeding,
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
            rnd=rnd,  # only used by SeedingAlgorithm.TruthSmeared
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
                seedDeduplication=(
                    True if seeding != SeedingAlgorithm.TruthSmeared else False
                ),
                stayOnSeed=True if seeding != SeedingAlgorithm.TruthSmeared else False,
            ),
            outputDirRoot=tp,
        )

        s.run()

        for file in (
            ["performance_seeding.root"]
            if seeding != SeedingAlgorithm.TruthSmeared
            else []
        ) + [
            "tracksummary_ckf.root",
            "performance_finding_ckf.root",
            "performance_fitting_ckf.root",
        ]:
            perf_file = tp / file
            assert perf_file.exists(), f"Performance file not found {perf_file}"
            (setup.outdir / label).mkdir(parents=True, exist_ok=True)
            shutil.copy(perf_file, setup.outdir / label / file)


for label, seeding in [
    ("truth_smeared", SeedingAlgorithm.TruthSmeared),
    ("truth_estimated", SeedingAlgorithm.TruthEstimated),
    ("seeded", SeedingAlgorithm.Default),
    ("orthogonal", SeedingAlgorithm.Orthogonal),
]:
    run_ckf_tracking(label, seeding)
