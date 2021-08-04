import pytest

import acts

from acts.examples import (
    TutorialVertexFinderAlgorithm,
    AdaptiveMultiVertexFinderAlgorithm,
    VertexFitterAlgorithm,
    IterativeVertexFinderAlgorithm,
    SpacePointMaker,
    TrackFindingAlgorithm,
    SeedingAlgorithm,
    TrackParamsEstimationAlgorithm,
    EventGenerator,
    FatrasAlgorithm,
    MaterialMapping,
    TruthSeedSelector,
    TruthTrackFinder,
    ParticleSelector,
    TruthVertexFinder,
    ParticleSmearing,
    TrackSelector,
    TrackFittingAlgorithm,
    SurfaceSortingAlgorithm,
    ParticlesPrinter,
    HitsPrinter,
    TrackParametersPrinter,
    PropagationAlgorithm,
    DigitizationAlgorithm,
    SmearingAlgorithm,
    PlanarSteppingAlgorithm,
)

from acts.examples.hepmc3 import HepMCProcessExtractor

from acts.examples.geant4.hepmc3 import EventRecording
from acts.examples.geant4 import GeantinoRecording


@pytest.mark.parametrize(
    "alg",
    [
        TutorialVertexFinderAlgorithm,
        AdaptiveMultiVertexFinderAlgorithm,
        VertexFitterAlgorithm,
        IterativeVertexFinderAlgorithm,
        SpacePointMaker,
        TrackFindingAlgorithm,
        SeedingAlgorithm,
        TrackParamsEstimationAlgorithm,
        HepMCProcessExtractor,
        EventGenerator,
        FatrasAlgorithm,
        MaterialMapping,
        TruthSeedSelector,
        TruthTrackFinder,
        ParticleSelector,
        TruthVertexFinder,
        ParticleSmearing,
        TrackSelector,
        TrackFittingAlgorithm,
        SurfaceSortingAlgorithm,
        ParticlesPrinter,
        HitsPrinter,
        TrackParametersPrinter,
        PropagationAlgorithm,
        GeantinoRecording,
        PlanarSteppingAlgorithm,
        EventRecording,
    ],
)
def test_algorithm_interface(alg):
    assert hasattr(alg, "Config")


def test_special_algorithm_interfaces():
    # just assert they exists
    assert DigitizationAlgorithm
    assert SmearingAlgorithm
