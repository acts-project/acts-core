#!/usr/bin/env python3
import os
from pathlib import Path

import acts
import acts.examples

import acts.examples.hepmc3
import acts.examples.dd4hep
import acts.examples.geant4
import acts.examples.geant4.dd4hep
import acts.examples.geant4.hepmc3

u = acts.UnitConstants


def runEventRecording(geoFactory, outputDir, s=None):
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
        outputParticles="particles_input",
        randomNumbers=rnd,
    )

    s.addReader(evGen)

    dd4hepSvc = acts.examples.dd4hep.DD4hepGeometryService(
        xmlFileNames=["thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
    )

    erAlgCfg = acts.examples.geant4.hepmc3.EventRecording.Config(
        inputParticles=evGen.config.outputParticles,
        outputHepMcTracks="geant-event",
        seed1=43,
        seed2=44,
        detectorConstructionFactory=geoFactory,
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
    odd = (
        Path(__file__).parent
        / "../../../thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"
    ).resolve()

    dd4hepSvc = acts.examples.dd4hep.DD4hepGeometryService(xmlFileNames=[str(odd)])

    geoFactory = acts.examples.geant4.dd4hep.DD4hepDetectorConstructionFactory(
        dd4hepSvc
    )

    runEventRecording(geoFactory=geoFactory, outputDir=os.getcwd()).run()
