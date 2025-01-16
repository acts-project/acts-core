import pytest
import os
from pathlib import Path
import multiprocessing

from helpers import (
    geant4Enabled,
    edm4hepEnabled,
    AssertCollectionExistsAlg,
)

import acts
from acts import UnitConstants as u
from acts.examples import (
    RootParticleWriter,
    RootParticleReader,
    RootMaterialTrackReader,
    RootTrackSummaryReader,
    CsvParticleWriter,
    CsvParticleReader,
    CsvMeasurementWriter,
    CsvMeasurementReader,
    CsvSimHitWriter,
    CsvSimHitReader,
    Sequencer,
)
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory


@pytest.mark.root
def test_root_particle_reader(tmp_path, conf_const, ptcl_gun):
    # need to write out some particles first
    s = Sequencer(numThreads=1, events=10, logLevel=acts.logging.WARNING)
    evGen = ptcl_gun(s)

    file = tmp_path / "particles.root"
    s.addWriter(
        conf_const(
            RootParticleWriter,
            acts.logging.WARNING,
            inputParticles=evGen.config.outputParticles,
            filePath=str(file),
        )
    )

    s.run()

    # reset sequencer for reading

    s2 = Sequencer(numThreads=1, logLevel=acts.logging.WARNING)

    s2.addReader(
        conf_const(
            RootParticleReader,
            acts.logging.WARNING,
            outputParticles="particles_generated",
            filePath=str(file),
        )
    )

    alg = AssertCollectionExistsAlg(
        "particles_generated", "check_alg", acts.logging.WARNING
    )
    s2.addAlgorithm(alg)

    s2.run()

    assert alg.events_seen == 10


@pytest.mark.csv
def test_csv_particle_reader(tmp_path, conf_const, ptcl_gun):
    s = Sequencer(numThreads=1, events=10, logLevel=acts.logging.WARNING)
    evGen = ptcl_gun(s)

    out = tmp_path / "csv"

    out.mkdir()

    s.addWriter(
        conf_const(
            CsvParticleWriter,
            acts.logging.WARNING,
            inputParticles=evGen.config.outputParticles,
            outputStem="particle",
            outputDir=str(out),
        )
    )

    s.run()

    # reset the seeder
    s = Sequencer(numThreads=1, logLevel=acts.logging.WARNING)

    s.addReader(
        conf_const(
            CsvParticleReader,
            acts.logging.WARNING,
            inputDir=str(out),
            inputStem="particle",
            outputParticles="input_particles",
        )
    )

    alg = AssertCollectionExistsAlg(
        "input_particles", "check_alg", acts.logging.WARNING
    )

    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 10


@pytest.mark.parametrize(
    "reader",
    [RootParticleReader, RootTrackSummaryReader],
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
@pytest.mark.odd
@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
def test_root_material_track_reader(material_recording):
    input_tracks = material_recording / "geant4_material_tracks.root"
    assert input_tracks.exists()

    s = Sequencer(numThreads=1)

    s.addReader(
        RootMaterialTrackReader(
            level=acts.logging.INFO,
            fileList=[str(input_tracks)],
            outputMaterialTracks="material-tracks",
        )
    )

    alg = AssertCollectionExistsAlg(
        "material-tracks", "check_alg", acts.logging.WARNING
    )
    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 2


@pytest.mark.csv
def test_csv_meas_reader(tmp_path, fatras, trk_geo, conf_const):
    s = Sequencer(numThreads=1, events=10)
    evGen, simAlg, digiAlg = fatras(s)

    out = tmp_path / "csv"
    out.mkdir()

    s.addWriter(
        CsvMeasurementWriter(
            level=acts.logging.INFO,
            inputMeasurements=digiAlg.config.outputMeasurements,
            inputClusters=digiAlg.config.outputClusters,
            inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
            outputDir=str(out),
        )
    )

    # Write hits, so we can later construct the measurement-particles-map
    s.addWriter(
        CsvSimHitWriter(
            level=acts.logging.INFO,
            inputSimHits=simAlg.config.outputSimHits,
            outputDir=str(out),
            outputStem="hits",
        )
    )

    s.run()

    # read back in
    s = Sequencer(numThreads=1)

    s.addReader(
        CsvSimHitReader(
            level=acts.logging.INFO,
            outputSimHits=simAlg.config.outputSimHits,
            inputDir=str(out),
            inputStem="hits",
        )
    )

    s.addReader(
        conf_const(
            CsvMeasurementReader,
            level=acts.logging.WARNING,
            outputMeasurements="measurements",
            outputMeasurementSimHitsMap="simhitsmap",
            outputMeasurementParticlesMap="meas_ptcl_map",
            inputSimHits=simAlg.config.outputSimHits,
            inputDir=str(out),
        )
    )

    algs = [
        AssertCollectionExistsAlg(k, f"check_alg_{k}", acts.logging.WARNING)
        for k in ("measurements", "simhitsmap", "meas_ptcl_map")
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

    alg = AssertCollectionExistsAlg("simhits", "check_alg", acts.logging.WARNING)
    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 10


def generate_input_test_edm4hep_simhit_reader(input, output):
    from DDSim.DD4hepSimulation import DD4hepSimulation

    ddsim = DD4hepSimulation()
    if isinstance(ddsim.compactFile, list):
        ddsim.compactFile = [input]
    else:
        ddsim.compactFile = input
    ddsim.enableGun = True
    ddsim.gun.direction = (1, 0, 0)
    ddsim.gun.particle = "pi-"
    ddsim.gun.distribution = "eta"
    ddsim.numberOfEvents = 10
    ddsim.outputFile = output
    ddsim.run()


@pytest.mark.slow
@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_simhit_particle_reader(tmp_path):
    from acts.examples.edm4hep import EDM4hepReader

    tmp_file = str(tmp_path / "output_edm4hep.root")
    odd_xml_file = str(getOpenDataDetectorDirectory() / "xml" / "OpenDataDetector.xml")

    # explicitly ask for "spawn" as CI failures were observed with "fork"
    spawn_context = multiprocessing.get_context("spawn")
    p = spawn_context.Process(
        target=generate_input_test_edm4hep_simhit_reader, args=(odd_xml_file, tmp_file)
    )
    p.start()
    p.join()

    assert os.path.exists(tmp_file)

    s = Sequencer(numThreads=1)

    with getOpenDataDetector() as detector:
        trackingGeometry = detector.trackingGeometry()

        s.addReader(
            EDM4hepReader(
                level=acts.logging.INFO,
                inputPath=tmp_file,
                inputSimHits=[
                    "PixelBarrelReadout",
                    "PixelEndcapReadout",
                    "ShortStripBarrelReadout",
                    "ShortStripEndcapReadout",
                    "LongStripBarrelReadout",
                    "LongStripEndcapReadout",
                ],
                outputParticlesGenerator="particles_generated",
                outputParticlesSimulation="particles_simulated",
                outputSimHits="simhits",
                dd4hepDetector=detector,
                trackingGeometry=trackingGeometry,
            )
        )

        alg = AssertCollectionExistsAlg("simhits", "check_alg", acts.logging.WARNING)
        s.addAlgorithm(alg)

        alg = AssertCollectionExistsAlg(
            "particles_generated", "check_alg", acts.logging.WARNING
        )
        s.addAlgorithm(alg)

        s.run()

    assert alg.events_seen == 10


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_measurement_reader(tmp_path, fatras, conf_const):
    from acts.examples.edm4hep import (
        EDM4hepMeasurementWriter,
        EDM4hepMeasurementReader,
    )

    s = Sequencer(numThreads=1, events=10)
    _, simAlg, digiAlg = fatras(s)

    out = tmp_path / "measurements_edm4hep.root"

    config = EDM4hepMeasurementWriter.Config(
        inputMeasurements=digiAlg.config.outputMeasurements,
        inputClusters=digiAlg.config.outputClusters,
        outputPath=str(out),
    )
    s.addWriter(EDM4hepMeasurementWriter(level=acts.logging.INFO, config=config))
    s.run()

    # read back in
    s = Sequencer(numThreads=1)

    s.addReader(
        conf_const(
            EDM4hepMeasurementReader,
            level=acts.logging.WARNING,
            outputMeasurements="measurements",
            outputMeasurementSimHitsMap="simhitsmap",
            inputPath=str(out),
        )
    )

    algs = [
        AssertCollectionExistsAlg(k, f"check_alg_{k}", acts.logging.WARNING)
        for k in ("measurements", "simhitsmap")
    ]
    for alg in algs:
        s.addAlgorithm(alg)

    s.run()

    for alg in algs:
        assert alg.events_seen == 10


@pytest.mark.edm4hep
@pytest.mark.skipif(not edm4hepEnabled, reason="EDM4hep is not set up")
def test_edm4hep_tracks_reader(tmp_path):
    from acts.examples.edm4hep import EDM4hepTrackWriter, EDM4hepTrackReader

    detector = acts.examples.GenericDetector()
    trackingGeometry = detector.trackingGeometry()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    from truth_tracking_kalman import runTruthTrackingKalman

    s = Sequencer(numThreads=1, events=10)
    runTruthTrackingKalman(
        trackingGeometry,
        field,
        digiConfigFile=Path(
            str(
                Path(__file__).parent.parent.parent.parent
                / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
            )
        ),
        outputDir=tmp_path,
        s=s,
    )

    out = tmp_path / "tracks_edm4hep.root"

    s.addWriter(
        EDM4hepTrackWriter(
            level=acts.logging.VERBOSE,
            inputTracks="kf_tracks",
            outputPath=str(out),
            Bz=2 * u.T,
        )
    )

    s.run()

    s = Sequencer(numThreads=1)
    s.addReader(
        EDM4hepTrackReader(
            level=acts.logging.VERBOSE,
            outputTracks="kf_tracks",
            inputPath=str(out),
            Bz=2 * u.T,
        )
    )

    s.run()
