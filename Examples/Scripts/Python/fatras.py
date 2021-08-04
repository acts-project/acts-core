#!/usr/bin/env python3
import os

import acts
import acts.examples

u = acts.UnitConstants


def runFatras(trackingGeometry, field, outputDir, s: acts.examples.Sequencer = None):

    # Preliminaries
    rnd = acts.examples.RandomNumbers()

    # Input
    vtxGen = acts.examples.GaussianVertexGenerator()
    vtxGen.stddev = acts.Vector4(0, 0, 0, 0)

    ptclGen = acts.examples.ParametricParticleGenerator(
        p=(1 * u.GeV, 10 * u.GeV), eta=(-2, 2)
    )

    g = acts.examples.EventGenerator.Generator()
    g.multiplicity = acts.examples.FixedMultiplicityGenerator()
    g.vertex = vtxGen
    g.particles = ptclGen

    evGen = acts.examples.EventGenerator(
        level=acts.logging.INFO,
        generators=[g],
        outputParticles="particles_input",
        randomNumbers=rnd,
    )

    # Selector
    selector = acts.examples.ParticleSelector(
        level=acts.logging.INFO,
        inputParticles=evGen.config.outputParticles,
        outputParticles="particles_selected",
    )

    # Simulation
    alg = acts.examples.FatrasAlgorithm(
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

    # Sequencer
    s = s or acts.examples.Sequencer(
        events=100, numThreads=-1, logLevel=acts.logging.INFO
    )

    s.addReader(evGen)
    s.addAlgorithm(selector)
    s.addAlgorithm(alg)

    # Output
    s.addWriter(
        acts.examples.CsvParticleWriter(
            level=acts.logging.INFO,
            outputDir=outputDir + "/csv",
            inputParticles="particles_final",
            outputStem="particles_final",
        )
    )

    s.addWriter(
        acts.examples.RootParticleWriter(
            level=acts.logging.INFO,
            inputParticles="particles_final",
            filePath=outputDir + "/fatras_particles_final.root",
        )
    )

    s.addWriter(
        acts.examples.CsvParticleWriter(
            level=acts.logging.INFO,
            outputDir=outputDir + "/csv",
            inputParticles="particles_initial",
            outputStem="particles_initial",
        )
    )

    s.addWriter(
        acts.examples.RootParticleWriter(
            level=acts.logging.INFO,
            inputParticles="particles_initial",
            filePath=outputDir + "/fatras_particles_initial.root",
        )
    )

    s.addWriter(
        acts.examples.CsvSimHitWriter(
            level=acts.logging.INFO,
            inputSimHits=alg.config.outputSimHits,
            outputDir=outputDir + "/csv",
            outputStem="hits",
        )
    )

    s.addWriter(
        acts.examples.RootSimHitWriter(
            level=acts.logging.INFO,
            inputSimHits=alg.config.outputSimHits,
            filePath=outputDir + "/hits.root",
        )
    )

    return s


if "__main__" == __name__:

    gdc = acts.examples.GenericDetector.Config()
    detector = acts.examples.GenericDetector()
    trackingGeometry, contextDecorators = detector.finalize(gdc, None)

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runFatras(trackingGeometry, field, os.getcwd()).run()
