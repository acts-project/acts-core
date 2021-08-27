from typing import Type
import inspect

import pytest

from helpers import geant4Enabled, AssertCollectionExistsAlg

import acts

from acts import PlanarModuleStepper
from acts.examples import (
    RootParticleWriter,
    RootParticleReader,
    RootMaterialTrackReader,
    RootTrajectorySummaryReader,
    CsvParticleWriter,
    CsvParticleReader,
    CsvMeasurementWriter,
    CsvMeasurementReader,
    CsvSimHitWriter,
    CsvSimHitReader,
    CsvPlanarClusterWriter,
    CsvPlanarClusterReader,
    PlanarSteppingAlgorithm,
    BareAlgorithm,
    Sequencer,
)


def test_root_particle_reader(tmp_path, conf_const, ptcl_gun):
    # need to write out some particles first
    s = Sequencer(numThreads=1, events=10, logLevel=acts.logging.ERROR)
    evGen = ptcl_gun(s)

    file = tmp_path / "particles.root"
    s.addWriter(
        conf_const(
            RootParticleWriter,
            acts.logging.ERROR,
            inputParticles=evGen.config.outputParticles,
            filePath=str(file),
        )
    )

    s.run()

    del s  # to properly close the root file

    # reset sequencer for reading

    s2 = Sequencer(numThreads=1, logLevel=acts.logging.ERROR)

    s2.addReader(
        conf_const(
            RootParticleReader,
            acts.logging.ERROR,
            particleCollection="input_particles",
            filePath=str(file),
        )
    )

    alg = AssertCollectionExistsAlg("input_particles", "check_alg", acts.logging.ERROR)
    s2.addAlgorithm(alg)

    s2.run()

    assert alg.events_seen == 10


def test_csv_particle_reader(tmp_path, conf_const, ptcl_gun):
    s = Sequencer(numThreads=1, events=10, logLevel=acts.logging.ERROR)
    evGen = ptcl_gun(s)

    out = tmp_path / "csv"

    out.mkdir()

    s.addWriter(
        conf_const(
            CsvParticleWriter,
            acts.logging.ERROR,
            inputParticles=evGen.config.outputParticles,
            outputStem="particle",
            outputDir=str(out),
        )
    )

    s.run()

    # reset the seeder
    s = Sequencer(numThreads=1, logLevel=acts.logging.ERROR)

    s.addReader(
        conf_const(
            CsvParticleReader,
            acts.logging.ERROR,
            inputDir=str(out),
            inputStem="particle",
            outputParticles="input_particles",
        )
    )

    alg = AssertCollectionExistsAlg("input_particles", "check_alg", acts.logging.ERROR)

    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 10


@pytest.mark.parametrize(
    "reader",
    [RootParticleReader, RootTrajectorySummaryReader],
)
@pytest.mark.root
def test_root_reader_interface(reader, conf_const, tmp_path):
    assert hasattr(reader, "Config")

    config = reader.Config

    assert hasattr(config, "filePath")

    kw = {"level": acts.logging.INFO, "filePath": str(tmp_path / "file.root")}

    assert conf_const(reader, **kw)


@pytest.mark.slow
@pytest.mark.root
@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
def test_root_material_track_reader(tmp_path, geantino_recording):

    # recreate sequencer

    s = Sequencer(numThreads=1)

    s.addReader(
        RootMaterialTrackReader(
            level=acts.logging.INFO,
            fileList=[str(geantino_recording / "geant-material-tracks.root")],
        )
    )

    alg = AssertCollectionExistsAlg("material-tracks", "check_alg", acts.logging.ERROR)
    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 200


@pytest.mark.csv
def test_csv_meas_reader(tmp_path, fatras, trk_geo, conf_const):
    s = Sequencer(numThreads=1, events=10)
    evGen, simAlg, digiAlg = fatras(s)

    out = tmp_path / "csv"
    out.mkdir()

    config = CsvMeasurementWriter.Config(
        inputMeasurements=digiAlg.config.outputMeasurements,
        inputClusters=""
        if digiAlg.config.isSimpleSmearer
        else digiAlg.config.outputClusters,
        inputSimHits=simAlg.config.outputSimHits,
        inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
        outputDir=str(out),
    )
    s.addWriter(CsvMeasurementWriter(level=acts.logging.INFO, config=config))
    s.run()

    # read back in
    s = Sequencer(numThreads=1)

    s.addReader(
        conf_const(
            CsvMeasurementReader,
            level=acts.logging.ERROR,
            outputMeasurements="measurements",
            outputMeasurementSimHitsMap="simhitsmap",
            outputSourceLinks="sourcelinks",
            inputDir=str(out),
        )
    )

    algs = [
        AssertCollectionExistsAlg(k, f"check_alg_{k}", acts.logging.ERROR)
        for k in ("measurements", "simhitsmap", "sourcelinks")
    ]
    for alg in algs:
        s.addAlgorithm(alg)

    s.run()

    for alg in algs:
        assert alg.events_seen == 10


@pytest.mark.csv
def test_csv_simhits_reader(tmp_path, fatras, conf_const):
    s = Sequencer(numThreads=1, events=10)
    evGen, simAlg, digiAlg = fatras(s)

    out = tmp_path / "csv"
    out.mkdir()

    s.addWriter(
        CsvSimHitWriter(
            level=acts.logging.INFO,
            inputSimHits=simAlg.config.outputSimHits,
            outputDir=str(out),
            outputStem="hits",
        )
    )

    s.run()

    s = Sequencer(numThreads=1)

    s.addReader(
        conf_const(
            CsvSimHitReader,
            level=acts.logging.INFO,
            inputDir=str(out),
            inputStem="hits",
            outputSimHits="simhits",
        )
    )

    alg = AssertCollectionExistsAlg("simhits", "check_alg", acts.logging.ERROR)
    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 10


@pytest.mark.csv
def test_csv_clusters_reader(tmp_path, fatras, conf_const, trk_geo, rng):
    s = Sequencer(numThreads=1, events=10)  # we're not going to use this one
    evGen, simAlg, _ = fatras(s)
    s = Sequencer(numThreads=1, events=10)
    s.addReader(evGen)
    s.addAlgorithm(simAlg)
    digiAlg = PlanarSteppingAlgorithm(
        level=acts.logging.ERROR,
        inputSimHits=simAlg.config.outputSimHits,
        outputClusters="clusters",
        outputSourceLinks="sourcelinks",
        outputMeasurements="measurements",
        outputMeasurementParticlesMap="meas_ptcl_map",
        outputMeasurementSimHitsMap="meas_sh_map",
        trackingGeometry=trk_geo,
        randomNumbers=rng,
        planarModuleStepper=PlanarModuleStepper(),
    )
    s.addAlgorithm(digiAlg)

    out = tmp_path / "csv"
    out.mkdir()

    s.addWriter(
        CsvPlanarClusterWriter(
            level=acts.logging.ERROR,
            outputDir=str(out),
            inputSimHits=simAlg.config.outputSimHits,
            inputClusters=digiAlg.config.outputClusters,
            trackingGeometry=trk_geo,
        )
    )

    s.run()

    s = Sequencer(numThreads=1)

    s.addReader(
        conf_const(
            CsvPlanarClusterReader,
            level=acts.logging.ERROR,
            outputClusters="clusters",
            inputDir=str(out),
            outputHitIds="hits",
            outputMeasurementParticlesMap="meas_ptcl_map",
            outputSimHits="simhits",
            trackingGeometry=trk_geo,
        )
    )

    algs = [
        AssertCollectionExistsAlg(k, f"check_alg_{k}", acts.logging.ERROR)
        for k in ("clusters", "simhits", "meas_ptcl_map")
    ]
    for alg in algs:
        s.addAlgorithm(alg)

    s.run()

    for alg in algs:
        assert alg.events_seen == 10
