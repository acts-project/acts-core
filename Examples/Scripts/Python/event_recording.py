#!/usr/bin/env python3

import os

import acts
import acts.examples
import acts.examples.hepmc3
import acts.examples.dd4hep
import acts.examples.geant4
import acts.examples.geant4.hepmc3
from acts.examples.odd import getOpenDataDetector


u = acts.UnitConstants


def runEventRecording(detector, outputDir, s=None):
    hepmc_dir = os.path.join(outputDir, "hepmc3")
    if not os.path.exists(hepmc_dir):
        os.mkdir(hepmc_dir)

    s = s or acts.examples.Sequencer(
        events=int(os.environ.get("NEVENTS", 100)), numThreads=1
    )

    rnd = acts.examples.RandomNumbers(seed=42)
    evGen = acts.examples.EventGenerator(
        level=acts.logging.INFO,
        generators=[
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=2),
                vertex=acts.examples.GaussianVertexGenerator(
                    stddev=acts.Vector4(0, 0, 0, 0), mean=acts.Vector4(0, 0, 0, 0)
                ),
                particles=acts.examples.ParametricParticleGenerator(
                    p=(1 * u.GeV, 10 * u.GeV),
                    eta=(-2, 2),
                    phi=(0, 90 * u.degree),
                    randomizeCharge=True,
                    numParticles=4,
                ),
            )
        ],
        outputParticles="particles_generated",
        outputVertices="vertices_input",
        randomNumbers=rnd,
    )

    s.addReader(evGen)

    erAlgCfg = acts.examples.geant4.hepmc3.EventRecording.Config(
        inputParticles=evGen.config.outputParticles,
        outputHepMcTracks="geant-event",
        seed1=43,
        seed2=44,
        detector=detector,
    )

    erAlg = acts.examples.geant4.hepmc3.EventRecording(
        config=erAlgCfg, level=acts.logging.INFO
    )

    s.addAlgorithm(erAlg)

    s.addWriter(
        acts.examples.hepmc3.HepMC3AsciiWriter(
            level=acts.logging.INFO,
            outputDir=hepmc_dir,
            outputStem="events",
            inputEvents=erAlg.config.outputHepMcTracks,
        )
    )

    return s


if "__main__" == __name__:
    detector = getOpenDataDetector()

    runEventRecording(
        detector=detector,
        outputDir=os.getcwd(),
    ).run()
