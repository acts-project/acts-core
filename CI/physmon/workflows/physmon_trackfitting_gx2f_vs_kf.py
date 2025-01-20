#!/usr/bin/env python3

import tempfile

from pathlib import Path

import acts
import acts.examples

from acts.examples.simulation import (
    addParticleGun,
    ParticleConfig,
    EtaConfig,
    PhiConfig,
    MomentumConfig,
    ParticleSelectorConfig,
    addFatras,
    addDigitization,
)
from acts.examples.reconstruction import (
    addSeeding,
    SeedingAlgorithm,
    addGx2fTracks,
    addKalmanTracks,
)

from physmon_common import makeSetup

u = acts.UnitConstants

setup = makeSetup()

trackingGeometry = setup.trackingGeometry
field = setup.field
digiConfigFile = setup.digiConfig

with tempfile.TemporaryDirectory() as temp:
    s = acts.examples.Sequencer(
        events=100000,
        numThreads=-1,
        logLevel=acts.logging.INFO,
        skip=0,
        trackFpes=False,
    )

    rnd = acts.examples.RandomNumbers(seed=42)
    outputDirTemp = Path(temp)
    outputDir = setup.outdir

    addParticleGun(
        s,
        ParticleConfig(num=1, pdg=acts.PdgParticle.eMuon, randomizeCharge=True),
        EtaConfig(-3.0, 3.0, uniform=True),
        MomentumConfig(1.0 * u.GeV, 100.0 * u.GeV, transverse=True),
        PhiConfig(0.0, 360.0 * u.degree),
        vtxGen=acts.examples.GaussianVertexGenerator(
            mean=acts.Vector4(0, 0, 0, 0),
            stddev=acts.Vector4(0, 0, 0, 0),
        ),
        multiplicity=1,
        rnd=rnd,
        outputDirRoot=outputDirTemp,
    )

    addFatras(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        enableInteractions=True,
        postSelectParticles=ParticleSelectorConfig(
            pt=(0.9 * u.GeV, None),
            measurements=(7, None),
            removeNeutral=True,
            removeSecondaries=True,
        ),
    )

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        rnd=rnd,
    )

    addSeeding(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        inputParticles="particles_input",
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
        particleHypothesis=acts.ParticleHypothesis.muon,
    )

    addGx2fTracks(
        s,
        trackingGeometry,
        field,
        nUpdateMax=17,
        relChi2changeCutOff=1e-7,
        multipleScattering=True,
    )

    s.addAlgorithm(
        acts.examples.TrackSelectorAlgorithm(
            level=acts.logging.INFO,
            inputTracks="tracks",
            outputTracks="selected-tracks-gx2f",
            selectorConfig=acts.TrackSelector.Config(
                minMeasurements=7,
            ),
        )
    )

    s.addWriter(
        acts.examples.TrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="selected-tracks-gx2f",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "performance_trackfitting_gx2f.root"),
        )
    )

    addKalmanTracks(
        s,
        trackingGeometry,
        field,
        directNavigation=False,
        reverseFilteringMomThreshold=0 * u.GeV,
    )

    s.addAlgorithm(
        acts.examples.TrackSelectorAlgorithm(
            level=acts.logging.INFO,
            inputTracks="tracks",
            outputTracks="selected-tracks-kf",
            selectorConfig=acts.TrackSelector.Config(
                minMeasurements=7,
            ),
        )
    )

    s.addWriter(
        acts.examples.TrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="selected-tracks-kf",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "performance_trackfitting_kf.root"),
        )
    )

    s.run()
    del s
