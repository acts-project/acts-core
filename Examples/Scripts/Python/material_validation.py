#!/usr/bin/env python3

import os
import argparse

import acts
from acts.examples import Sequencer, RootMaterialTrackWriter
from acts.examples.odd import getOpenDataDetector
from acts.examples.simulation import addParticleGun, EtaConfig, ParticleConfig


def runMaterialValidation(
    nevents,
    ntracks,
    trackingGeometry,
    decorators,
    field,
    outputDir,
    outputName="propagation-material",
    s=None,
):
    # Create a sequencer
    s = s or Sequencer(events=nevents, numThreads=-1)

    for decorator in decorators:
        s.addContextDecorator(decorator)

    nav = acts.Navigator(trackingGeometry=trackingGeometry)
    stepper = acts.StraightLineStepper()
    # stepper = acts.EigenStepper(field)

    prop = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))

    rnd = acts.examples.RandomNumbers(seed=42)

    addParticleGun(
        s,
        ParticleConfig(num=ntracks, pdg=acts.PdgParticle.eMuon, randomizeCharge=True),
        EtaConfig(-4.0, 4.0),
        rnd=rnd,
    )

    trkParamExtractor = acts.examples.ParticleTrackParamExtractor(
        level=acts.logging.INFO,
        inputParticles="particles_generated",
        outputTrackParameters="params_particles_generated",
    )
    s.addAlgorithm(trkParamExtractor)

    alg = acts.examples.PropagationAlgorithm(
        propagatorImpl=prop,
        level=acts.logging.INFO,
        sterileLogger=True,
        recordMaterialInteractions=True,
        inputTrackParameters="params_particles_generated",
        outputSummaryCollection="propagation_summary",
        outputMaterialCollection="material_tracks",
    )
    s.addAlgorithm(alg)

    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks=alg.config.outputMaterialCollection,
            filePath=os.path.join(outputDir, (outputName + ".root")),
            storeSurface=True,
            storeVolume=True,
        )
    )

    return s


if "__main__" == __name__:
    p = argparse.ArgumentParser()

    p.add_argument(
        "-n", "--events", type=int, default=1000, help="Number of events to process"
    )
    p.add_argument(
        "-t", "--tracks", type=int, default=1000, help="Number of tracks per event"
    )
    p.add_argument(
        "-m", "--map", type=str, help="Input file (optional) for the material map"
    )
    p.add_argument(
        "-o",
        "--output",
        type=str,
        default="propagation-material",
        help="Output file name",
    )

    args = p.parse_args()

    materialDecorator = (
        acts.IMaterialDecorator.fromFile(args.map) if args.map != None else None
    )

    detector = getOpenDataDetector(materialDecorator)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    field = acts.ConstantBField(acts.Vector3(0, 0, 0 * acts.UnitConstants.T))

    runMaterialValidation(
        args.events,
        args.tracks,
        trackingGeometry,
        decorators,
        field,
        outputDir=os.getcwd(),
        outputName=args.output,
    ).run()
