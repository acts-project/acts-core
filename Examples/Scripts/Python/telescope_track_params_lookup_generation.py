#!/usr/bin/env python3

import argparse

import acts
import acts.examples
from acts.examples.simulation import (
    addParticleGun,
    addFatras,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
)

u = acts.UnitConstants


def estimateLookup(trackingGeometry, numEvents, outputPath):

    # Set up the dipole magnetic field
    field = acts.ConstantBField(acts.Vector3(50 * u.T, 0, 0))

    # Fatras simulation of muons
    rnd = acts.examples.RandomNumbers(seed=42)

    s = acts.examples.Sequencer(
        events=numEvents, numThreads=1, logLevel=acts.logging.INFO
    )

    vertexGen = acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(0, 0, 0, 0), mean=acts.Vector4(0, 9, 0, 0)
    )

    addParticleGun(
        s=s,
        etaConfig=EtaConfig(10.0, 10.0),
        phiConfig=PhiConfig(0, 0),
        momentumConfig=MomentumConfig(0.5 * u.GeV, 10 * u.GeV),
        particleConfig=ParticleConfig(1, acts.PdgParticle.eMuon, False),
        multiplicity=1,
        rnd=rnd,
        vtxGen=vertexGen,
    )

    addFatras(
        s,
        trackingGeometry,
        field,
        inputParticles="particles_generated",
        outputSimHits="sim_hits",
        rnd=rnd,
    )

    # Set up the track lookup grid writer
    jsonWriterConfig = acts.examples.JsonTrackParamsLookupWriter.Config(path=outputPath)
    jsonWriter = acts.examples.JsonTrackParamsLookupWriter(jsonWriterConfig)

    # Set up the track estimation algorithm
    surfaces = list(trackingGeometry.geoIdSurfaceMap().values())
    refSurface = surfaces[0]
    refGeometryId = refSurface.geometryId()

    trackEstConfig = acts.examples.TrackParamsLookupEstimation.Config(
        refLayers={refGeometryId: refSurface},
        bins=(1, 1000),
        inputHits="sim_hits",
        inputParticles="particles_generated",
        trackLookupGridWriters=[jsonWriter],
    )
    trackEstAlg = acts.examples.TrackParamsLookupEstimation(
        trackEstConfig, acts.logging.INFO
    )

    s.addAlgorithm(trackEstAlg)

    s.run()


if __name__ == "__main__":
    p = argparse.ArgumentParser()

    p.add_argument(
        "-n",
        "--events",
        type=int,
        default=100000,
        help="Number of events for lookup estimation",
    )
    p.add_argument(
        "-o",
        "--output",
        type=str,
        default="lookup.json",
        help="Output lookup file name",
    )

    args = p.parse_args()

    # Initialize the geometry
    detector = acts.examples.TelescopeDetector(
        bounds=[4, 10],
        positions=[30, 60, 90],
        stereos=[0, 0, 0],
        binValue=2,
        surfaceType=0,
    )
    trackingGeometry = detector.trackingGeometry()

    # Estimate the lookup
    estimateLookup(trackingGeometry, args.events, args.output)
