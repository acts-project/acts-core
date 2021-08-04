#!/usr/bin/env python3
import os

import acts
import acts.examples

from common import addPythia8

u = acts.UnitConstants


def runPythia8(
    outputDir,
    outputRoot: bool = True,
    outputCsv: bool = True,
    s: acts.examples.Sequencer = None,
):
    # Preliminaries
    rnd = acts.examples.RandomNumbers()

    # Sequencer
    s = s or acts.examples.Sequencer(
        events=10, numThreads=-1, logLevel=acts.logging.INFO
    )

    evGen = addPythia8(s, rnd)

    if outputRoot:
        s.addWriter(
            acts.examples.RootParticleWriter(
                level=acts.logging.INFO,
                inputParticles=evGen.config.outputParticles,
                filePath=outputDir + "/pythia8_particles.root",
            )
        )

    if outputCsv:
        s.addWriter(
            acts.examples.CsvParticleWriter(
                level=acts.logging.INFO,
                inputParticles=evGen.config.outputParticles,
                outputDir=outputDir + "/csv/",
                outputStem="particles",
            )
        )

    return s


if "__main__" == __name__:
    runPythia8(os.getcwd()).run()
