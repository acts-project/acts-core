#!/usr/bin/env python3
from pathlib import Path
from typing import Optional
import enum
import warnings

from acts.examples import (
    Sequencer,
    ParticleSelector,
    ParticleSmearing,
    TruthVertexFinder,
    VertexFitterAlgorithm,
    IterativeVertexFinderAlgorithm,
    RootParticleReader,
    AdaptiveMultiVertexFinderAlgorithm,
    RootVertexPerformanceWriter,
    RootTrajectorySummaryReader,
    TrackSelector,
    GenericDetector,
)

import acts

from acts import UnitConstants as u

from common import addPythia8


class VertexFinder(enum.Enum):
    Truth = (1,)
    AMVF = (2,)
    Iterative = (3,)


def runVertexFitting(
    field,
    outputDir: Path,
    outputRoot: bool = True,
    inputParticlePath: Optional[Path] = None,
    inputTrackSummary: Path = None,
    vertexFinder: VertexFinder = VertexFinder.Truth,
    s=None,
):
    s = s or Sequencer(events=100, numThreads=1)

    logger = acts.logging.getLogger("VertexFittingExample")

    rnd = acts.examples.RandomNumbers(seed=42)

    if inputParticlePath is None:
        logger.info("Generating particles using Pythia8")
        evGen = addPythia8(s, rnd)
        inputParticles = evGen.config.outputParticles
    else:
        logger.info("Reading particles from %s", inputParticlePath.resolve())
        assert inputParticlePath.exists()
        inputParticles = "particles_read"
        s.addReader(
            RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath.resolve()),
                particleCollection=inputParticles,
                orderedEvents=False,
            )
        )

    selectedParticles = "particles_selected"
    ptclSelector = ParticleSelector(
        level=acts.logging.INFO,
        inputParticles=inputParticles,
        outputParticles=selectedParticles,
        removeNeutral=True,
        absEtaMax=2.5,
        rhoMax=4.0 * u.mm,
        ptMin=500 * u.MeV,
    )
    s.addAlgorithm(ptclSelector)

    trackParameters = "trackparameters"
    if inputTrackSummary is None:
        logger.info("Using smeared particles")

        ptclSmearing = ParticleSmearing(
            level=acts.logging.INFO,
            inputParticles=selectedParticles,
            outputTrackParameters=trackParameters,
            randomNumbers=rnd,
        )
        s.addAlgorithm(ptclSmearing)
        associatedParticles = selectedParticles
    else:
        logger.info("Reading track summary from %s", inputTrackSummary.resolve())
        assert inputTrackSummary.exists()
        associatedParticles = "associatedTruthParticles"
        trackSummaryReader = RootTrajectorySummaryReader(
            level=acts.logging.VERBOSE,
            outputTracks="fittedTrackParameters",
            outputParticles=associatedParticles,
            filePath=str(inputTrackSummary.resolve()),
            orderedEvents=False,
        )
        s.addReader(trackSummaryReader)

        s.addAlgorithm(
            TrackSelector(
                level=acts.logging.INFO,
                inputTrackParameters=trackSummaryReader.config.outputTracks,
                outputTrackParameters=trackParameters,
                outputTrackIndices="outputTrackIndices",
                removeNeutral=True,
                absEtaMax=2.5,
                loc0Max=4.0 * u.mm,  # rho max
                ptMin=1 * u.GeV,
            )
        )

    logger.info("Using vertex finder: %s", vertexFinder.name)

    outputVertices = "fittedVertices"
    outputTime = ""
    if vertexFinder == VertexFinder.Truth:
        findVertices = TruthVertexFinder(
            level=acts.logging.VERBOSE,
            inputParticles=selectedParticles,
            outputProtoVertices="protovertices",
            excludeSecondaries=True,
        )
        s.addAlgorithm(findVertices)
        fitVertices = VertexFitterAlgorithm(
            level=acts.logging.VERBOSE,
            bField=field,
            inputTrackParameters=trackParameters,
            inputProtoVertices=findVertices.config.outputProtoVertices,
            outputVertices=outputVertices,
        )
        s.addAlgorithm(fitVertices)

    elif vertexFinder == VertexFinder.Iterative:
        findVertices = IterativeVertexFinderAlgorithm(
            level=acts.logging.INFO,
            bField=field,
            inputTrackParameters=trackParameters,
            outputProtoVertices="protovertices",
            outputVertices=outputVertices,
        )
        s.addAlgorithm(findVertices)
    elif vertexFinder == VertexFinder.AMVF:
        outputTime = "outputTime"
        findVertices = AdaptiveMultiVertexFinderAlgorithm(
            level=acts.logging.INFO,
            bField=field,
            inputTrackParameters=trackParameters,
            outputProtoVertices="protovertices",
            outputVertices=outputVertices,
            outputTime=outputTime,
        )

        s.addAlgorithm(findVertices)
    else:
        raise RuntimeError("Invalid finder argument")

    if outputRoot:
        if inputTrackSummary is None:
            warnings.warn(
                "This combination is not necessarily supported. Please get in touch with us"
            )
        s.addWriter(
            RootVertexPerformanceWriter(
                level=acts.logging.INFO,
                inputAllTruthParticles=inputParticles,
                inputSelectedTruthParticles=selectedParticles,
                inputAssociatedTruthParticles=associatedParticles,
                inputFittedTracks=trackParameters,
                inputVertices=outputVertices,
                inputTime=outputTime,
                treeName="vertexing",
                filePath=str(outputDir / "performance_vertexing.root"),
            )
        )

    return s


if "__main__" == __name__:
    detector, trackingGeometry, decorators = GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    inputParticlePath = Path("particles.root")
    if not inputParticlePath.exists():
        inputParticlePath = None
    inputParticlePath = None

    inputTrackSummary = None
    for p in ("tracksummary_fitter.root", "tracksummary_ckf.root"):
        p = Path(p)
        if p.exists():
            inputTrackSummary = p
            break

    runVertexFitting(
        field,
        vertexFinder=VertexFinder.Truth,
        inputParticlePath=inputParticlePath,
        inputTrackSummary=inputTrackSummary,
        outputDir=Path.cwd(),
    ).run()
