#!/usr/bin/env python3
import os

from acts.examples import (
    Sequencer,
    RandomNumbers,
    GaussianVertexGenerator,
    EventGenerator,
    FixedMultiplicityGenerator,
    ParametricParticleGenerator,
    CsvParticleWriter,
    ParticlesPrinter,
    RootParticleWriter,
)

import acts
from acts import Vector4, UnitConstants as u


def runParticleGun(outputDir, s=None):
    s = s or Sequencer(events=10, numThreads=-1)

    # Preliminaries
    rnd = RandomNumbers(seed=228)

    # Input
    vtxGen = GaussianVertexGenerator()
    vtxGen.stddev = Vector4(0, 0, 0, 0)

    ptclGen = ParametricParticleGenerator(
        p=(1 * u.GeV, 10 * u.GeV), eta=(-4, 4), numParticles=2
    )

    g = EventGenerator.Generator()
    g.multiplicity = FixedMultiplicityGenerator(n=1)
    g.vertex = vtxGen
    g.particles = ptclGen

    evGen = EventGenerator(
        level=acts.logging.INFO,
        generators=[g],
        outputParticles="particles_input",
        randomNumbers=rnd,
    )

    s.addReader(evGen)

    s.addAlgorithm(
        ParticlesPrinter(
            level=acts.logging.INFO, inputParticles=evGen.config.outputParticles
        )
    )

    csv_dir = os.path.join(outputDir, "csv")
    if not os.path.exists(csv_dir):
        os.mkdir(csv_dir)

    s.addWriter(
        CsvParticleWriter(
            level=acts.logging.INFO,
            inputParticles=evGen.config.outputParticles,
            outputDir=csv_dir,
            outputStem="particles",
        ),
    )

    root_file = os.path.join(outputDir, "particles.root")

    s.addWriter(
        RootParticleWriter(
            level=acts.logging.INFO,
            inputParticles=evGen.config.outputParticles,
            filePath=root_file,
        )
    )

    return s


if "__main__" == __name__:
    runParticleGun(os.getcwd()).run()
