#!/usr/bin/env python3
import os
import warnings

import acts
import acts.examples

import acts.examples.dd4hep
import acts.examples.geant4
import acts.examples.geant4.dd4hep

u = acts.UnitConstants

_geantino_recording_executed = False


def runGeantinoRecording(geoFactory, outputDir, tracksPerEvent=10000, s=None):
    global _geantino_recording_executed
    if _geantino_recording_executed:
        warnings.warn("Geantino recording already ran in this process. Expect crashes")
    _geantino_recording_executed = True

    s = s or acts.examples.Sequencer(events=10, numThreads=1)

    g4AlgCfg = acts.examples.geant4.GeantinoRecording.Config()
    g4AlgCfg.detectorConstructionFactory = geoFactory
    g4AlgCfg.tracksPerEvent = tracksPerEvent

    g4Alg = acts.examples.geant4.GeantinoRecording(
        level=acts.logging.INFO, config=g4AlgCfg
    )

    s.addAlgorithm(g4Alg)

    s.addWriter(
        acts.examples.RootMaterialTrackWriter(
            prePostStep=True,
            recalculateTotals=True,
            collection=g4Alg.config.outputMaterialTracks,
            filePath=outputDir + "/" + g4Alg.config.outputMaterialTracks + ".root",
            level=acts.logging.INFO,
        )
    )

    return s


if "__main__" == __name__:
    dd4hepSvc = acts.examples.dd4hep.DD4hepGeometryService(
        xmlFileNames=["thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
    )
    geoFactory = acts.examples.geant4.dd4hep.DD4hepDetectorConstructionFactory(
        dd4hepSvc
    )

    runGeantinoRecording(geoFactory=geoFactory, outputDir=os.getcwd()).run()
