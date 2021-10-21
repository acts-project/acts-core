from pathlib import Path
import os
import json
import functools
import subprocess

import pytest

from helpers import (
    geant4Enabled,
    rootEnabled,
    dd4hepEnabled,
    hepmc3Enabled,
    AssertCollectionExistsAlg,
    isCI,
)
from helpers.hash_root import assert_root_hash

pytestmark = pytest.mark.skipif(not rootEnabled, reason="ROOT not set up")


import acts
from acts.examples import (
    Sequencer,
    GenericDetector,
    AlignedDetector,
    RootParticleWriter,
)

from common import getOpenDataDetector

u = acts.UnitConstants


@pytest.fixture
def field():
    return acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))


@pytest.fixture
def seq():
    return Sequencer(events=10, numThreads=1)


def assert_csv_output(csv_path, stem):
    __tracebackhide__ = True
    # print(list(csv_path.iterdir()))
    assert len([f for f in csv_path.iterdir() if f.name.endswith(stem + ".csv")]) > 0
    assert all([f.stat().st_size > 100 for f in csv_path.iterdir()])


def assert_entries(root_file, tree_name, exp):
    __tracebackhide__ = True
    import ROOT

    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)

    rf = ROOT.TFile.Open(str(root_file))
    keys = [k.GetName() for k in rf.GetListOfKeys()]
    assert tree_name in keys
    assert rf.Get(tree_name).GetEntries() == exp, f"{root_file}:{tree_name}"


def test_fatras(trk_geo, tmp_path, field):
    from fatras import runFatras

    csv = tmp_path / "csv"
    csv.mkdir()

    nevents = 10

    root_files = [
        (
            "fatras_particles_final.root",
            "particles",
            nevents,
            "09b0d46ad7641fc5af97650cdcbf0e6c85966933a2e1b2bfcf7fbb9c6d90b2e1",
        ),
        (
            "fatras_particles_initial.root",
            "particles",
            nevents,
            "712a41d95ed1fccafccf708afddcd2177793934bb80cd40b1d99b532c7a21031",
        ),
        (
            "hits.root",
            "hits",
            115,
            "f890599a8c94bc80270636c64039f8a0a7d8ddd716a15306d5ca61f15f8da5e8",
        ),
    ]

    assert len(list(csv.iterdir())) == 0
    for rf, _, _, _ in root_files:
        assert not (tmp_path / rf).exists()

    seq = Sequencer(events=nevents)
    runFatras(trk_geo, field, str(tmp_path), s=seq).run()

    del seq

    assert_csv_output(csv, "particles_final")
    assert_csv_output(csv, "particles_initial")
    assert_csv_output(csv, "hits")
    for f, tn, exp_entries, file_hash in root_files:
        rfp = tmp_path / f
        assert rfp.exists()
        assert rfp.stat().st_size > 2 ** 10 * 10

        assert_entries(rfp, tn, exp_entries)
        assert_root_hash(rfp, file_hash)


def test_seeding(tmp_path, trk_geo, field):
    from seeding import runSeeding

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))

    csv = tmp_path / "csv"
    csv.mkdir()

    seq = Sequencer(events=10, numThreads=1)

    root_files = [
        (
            "estimatedparams.root",
            "estimatedparams",
            371,
            "8bbc97cb3d4777c61dd0b483a1c8268fc8411ad182c35bc731e5ed222450deca",
        ),
        (
            "performance_seeding_trees.root",
            "track_finder_tracks",
            371,
            "c900a70a9881ce637a73ff730de4b2f4b6c0446854d666f70c7fc9c44c893666",
        ),
        (
            "performance_seeding_hists.root",
            None,
            0,
            "4943f152bfad6ca6302563b57e495599fad4fea43b8e0b41abe5b8de31b391bc",
        ),
        (
            "evgen_particles.root",
            "particles",
            seq.config.events,
            "4943f152bfad6ca6302563b57e495599fad4fea43b8e0b41abe5b8de31b391bc",
        ),
        (
            "fatras_particles_final.root",
            "particles",
            seq.config.events,
            "934e290545090fe87d177bf40a78cf9375420e854433c2ef5a6d408fb3744d9d",
        ),
        (
            "fatras_particles_initial.root",
            "particles",
            seq.config.events,
            "4943f152bfad6ca6302563b57e495599fad4fea43b8e0b41abe5b8de31b391bc",
        ),
    ]

    for fn, _, _, _ in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    assert len(list(csv.iterdir())) == 0

    runSeeding(trk_geo, field, outputDir=str(tmp_path), s=seq).run()

    del seq

    for fn, tn, exp_entries, file_hash in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 100

        if tn is not None:
            assert_entries(fp, tn, exp_entries)
            assert_root_hash(fp, file_hash)

    assert_csv_output(csv, "evgen_particles")
    assert_csv_output(csv, "evgen_particles")
    assert_csv_output(csv, "fatras_particles_final")
    assert_csv_output(csv, "fatras_particles_initial")


def test_propagation(tmp_path, trk_geo, field, seq):
    from propagation import runPropagation

    obj = tmp_path / "obj"
    obj.mkdir()

    root_files = [
        (
            "propagation_steps.root",
            "propagation_steps",
            10000,
            "2d9a86e1097b8f0a845f33b05830f6b9c78cdeaeeb01685753e8dc9b434c7f9d",
        )
    ]

    for fn, _, _, _ in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    assert len(list(obj.iterdir())) == 0

    runPropagation(trk_geo, field, str(tmp_path), s=seq).run()

    for fn, tn, ee, file_hash in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 2 ** 10 * 50
        assert_entries(fp, tn, ee)
        assert_root_hash(fp, file_hash)

    assert len(list(obj.iterdir())) > 0


@pytest.mark.slow
@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
def test_material_recording(tmp_path, material_recording):

    # Not quite sure why this isn't 200
    root_files = [
        (
            "geant4_material_tracks.root",
            "material-tracks",
            198,
            "aaf899f3584d6a96905dce0251d6cfe060b00008aeca43ea8162d65aa982000b",
        )
    ]

    for fn, tn, ee, file_hash in root_files:
        fp = material_recording / fn
        assert fp.exists()
        assert fp.stat().st_size > 2 ** 10 * 50
        assert_entries(fp, tn, ee)
        assert_root_hash(fp, file_hash)


@pytest.mark.slow
@pytest.mark.skipif(not hepmc3Enabled, reason="HepMC3 plugin not available")
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
def test_event_recording(tmp_path):

    script = (
        Path(__file__).parent.parent.parent.parent
        / "Examples"
        / "Scripts"
        / "Python"
        / "event_recording.py"
    )
    assert script.exists()

    env = os.environ.copy()
    env["NEVENTS"] = "1"
    subprocess.check_call([str(script)], cwd=tmp_path, env=env)

    from acts.examples.hepmc3 import HepMC3AsciiReader

    out_path = tmp_path / "hepmc3"
    # out_path.mkdir()

    assert len([f for f in out_path.iterdir() if f.name.endswith("events.hepmc3")]) > 0
    assert all([f.stat().st_size > 100 for f in out_path.iterdir()])

    s = Sequencer(numThreads=1)

    s.addReader(
        HepMC3AsciiReader(
            level=acts.logging.INFO,
            inputDir=str(out_path),
            inputStem="events",
            outputEvents="hepmc-events",
        )
    )

    alg = AssertCollectionExistsAlg(
        "hepmc-events", name="check_alg", level=acts.logging.INFO
    )
    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 1


def test_particle_gun(tmp_path):
    from particle_gun import runParticleGun

    s = Sequencer(events=20, numThreads=-1)

    csv_dir = tmp_path / "csv"
    root_file = tmp_path / "particles.root"

    assert not csv_dir.exists()
    assert not root_file.exists()

    runParticleGun(str(tmp_path), s=s).run()

    assert csv_dir.exists()
    assert root_file.exists()

    assert len([f for f in csv_dir.iterdir() if f.name.endswith("particles.csv")]) > 0
    assert all([f.stat().st_size > 100 for f in csv_dir.iterdir()])

    assert root_file.stat().st_size > 200
    assert_entries(root_file, "particles", 20)
    assert_root_hash(
        root_file, "78a89f365177423d0834ea6f1bd8afe1488e72b12a25066a20bd9050f5407860"
    )


@pytest.mark.slow
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
def test_material_mapping(material_recording, tmp_path):
    map_file = tmp_path / "material-maps_tracks.root"
    assert not map_file.exists()

    s = Sequencer(numThreads=1)

    detector, trackingGeometry, decorators = getOpenDataDetector()

    from material_mapping import runMaterialMapping

    runMaterialMapping(
        trackingGeometry,
        decorators,
        outputDir=str(tmp_path),
        inputDir=material_recording,
        s=s,
    )

    s.run()

    # MaterialMapping alg only writes on destruct.
    # See https://github.com/acts-project/acts/issues/881
    del s

    mat_file = tmp_path / "material-map.json"

    assert mat_file.exists()
    assert mat_file.stat().st_size > 10

    with mat_file.open() as fh:
        assert json.load(fh)

    assert map_file.exists()
    assert_entries(map_file, "material-tracks", 198)
    assert_root_hash(
        map_file, "7daa335505eb0b0143372fe00c75f62d03a69f94930be53c74a3b8b6e39e7822"
    )

    val_file = tmp_path / "propagation-material.root"
    assert not val_file.exists()

    # test the validation as well

    # we need to destroy the ODD to reload with material
    # del trackingGeometry
    # del detector

    detector, trackingGeometry, decorators = getOpenDataDetector(
        mdecorator=acts.IMaterialDecorator.fromFile(mat_file)
    )

    from material_validation import runMaterialValidation

    s = Sequencer(events=10, numThreads=1)

    field = acts.NullBField()

    runMaterialValidation(
        trackingGeometry, decorators, field, outputDir=str(tmp_path), s=s
    )

    s.run()

    assert val_file.exists()
    assert_entries(val_file, "material-tracks", 10000)
    assert_root_hash(
        val_file, "1fd042f01bc2ff60f740c3179cc7a68a641d7311f0eb65f8937f67cedd254634"
    )


@pytest.mark.parametrize(
    "geoFactory,nobj",
    [
        (GenericDetector.create, 450),
        pytest.param(
            getOpenDataDetector,
            540,
            marks=pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up"),
        ),
        (functools.partial(AlignedDetector.create, iovSize=1), 450),
    ],
)
def test_geometry_example(geoFactory, nobj, tmp_path):
    detector, trackingGeometry, decorators = geoFactory()

    from geometry import runGeometry

    json_dir = tmp_path / "json"
    csv_dir = tmp_path / "csv"
    obj_dir = tmp_path / "obj"

    for d in (json_dir, csv_dir, obj_dir):
        d.mkdir()

    events = 5

    kwargs = dict(
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        events=events,
        outputDir=str(tmp_path),
    )

    runGeometry(outputJson=True, **kwargs)
    runGeometry(outputJson=False, **kwargs)

    assert len(list(obj_dir.iterdir())) == nobj
    assert all(f.stat().st_size > 200 for f in obj_dir.iterdir())

    assert len(list(csv_dir.iterdir())) == 3 * events
    assert all(f.stat().st_size > 200 for f in csv_dir.iterdir())

    detector_files = [csv_dir / f"event{i:>09}-detectors.csv" for i in range(events)]
    for detector_file in detector_files:
        assert detector_file.exists()
        assert detector_file.stat().st_size > 200

    contents = [f.read_text() for f in detector_files]
    ref = contents[0]
    for c in contents[1:]:
        if isinstance(detector, AlignedDetector):
            assert c != ref, "Detector writeout is expected to be different"
        else:
            assert c == ref, "Detector writeout is expected to be identical"

    if not isinstance(detector, AlignedDetector):
        for f in [json_dir / f"event{i:>09}-detector.json" for i in range(events)]:
            assert detector_file.exists()
            with f.open() as fh:
                data = json.load(fh)
                assert data
        material_file = tmp_path / "geometry-map.json"
        assert material_file.exists()
        assert material_file.stat().st_size > 200


def test_digitization_example(trk_geo, tmp_path):
    from digitization import configureDigitization

    s = Sequencer(events=10, numThreads=-1)

    csv_dir = tmp_path / "csv"
    root_file = tmp_path / "measurements.root"

    assert not root_file.exists()
    assert not csv_dir.exists()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    configureDigitization(trk_geo, field, outputDir=tmp_path, s=s)

    s.run()

    assert root_file.exists()
    assert csv_dir.exists()

    assert len(list(csv_dir.iterdir())) == 3 * s.config.events
    assert all(f.stat().st_size > 50 for f in csv_dir.iterdir())
    for tn, nev in (
        (8, 407),
        (9, 0),
        (12, 11),
        (13, 375),
        (14, 2),
        (16, 25),
        (17, 146),
        (18, 9),
    ):
        assert_entries(root_file, f"vol{tn}", nev)

    assert_root_hash(
        root_file, "f2dd54cd8315e4656136571d4802a69293d7feb71f427706410b8f0d2ad76265"
    )


@pytest.mark.xfail(
    reason="Digitization from input currently not reproducible", strict=True
)
def test_digitization_example_input(trk_geo, tmp_path):
    from particle_gun import runParticleGun
    from digitization import configureDigitization

    ptcl_dir = tmp_path / "ptcl"
    ptcl_dir.mkdir()
    pgs = Sequencer(events=20, numThreads=-1)
    runParticleGun(str(ptcl_dir), s=pgs)
    pgs.run()

    s = Sequencer(numThreads=-1)

    csv_dir = tmp_path / "csv"
    root_file = tmp_path / "measurements.root"

    assert not root_file.exists()
    assert not csv_dir.exists()

    assert_root_hash(
        ptcl_dir / "particles.root",
        "78a89f365177423d0834ea6f1bd8afe1488e72b12a25066a20bd9050f5407860",
    )

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    configureDigitization(
        trk_geo,
        field,
        outputDir=tmp_path,
        particlesInput=ptcl_dir / "particles.root",
        s=s,
    )

    s.run()

    assert root_file.exists()
    assert csv_dir.exists()

    assert len(list(csv_dir.iterdir())) == 3 * pgs.config.events
    assert all(f.stat().st_size > 50 for f in csv_dir.iterdir())
    for tn, nev in (
        (7, 0),
        (8, 193),
        (9, 0),
        (12, 1),
        (13, 183),
        (14, 6),
        (16, 3),
        (17, 76),
        (18, 10),
    ):
        assert_entries(root_file, f"vol{tn}", nev)
    assert_root_hash(
        root_file, "ccc92f0ad538d1b62d98f19f947970bcc491843e54d8ffeed16ad2e226b8caee"
    )


def test_digitization_config_example(trk_geo, tmp_path):
    from digitization_config import runDigitizationConfig

    out_file = tmp_path / "output.json"
    assert not out_file.exists()

    input = (
        Path(__file__).parent
        / "../../../Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
    )
    assert input.exists(), input.resolve()

    runDigitizationConfig(trk_geo, input=input, output=out_file)

    assert out_file.exists()

    with out_file.open() as fh:
        data = json.load(fh)
    assert len(data.keys()) == 2
    assert data["acts-geometry-hierarchy-map"]["format-version"] == 0
    assert (
        data["acts-geometry-hierarchy-map"]["value-identifier"]
        == "digitization-configuration"
    )
    assert len(data["entries"]) == 27


def test_ckf_tracks_example_full_seeding(tmp_path):
    # the example as written is only compatible with the generic detector
    detector, trackingGeometry, decorators = GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    events = 10
    s = Sequencer(events=events, numThreads=1)  # Digitization is not thread-safe

    root_files = [
        (
            "performance_ckf.root",
            None,
            None,
            None,
        ),
        (
            "performance_seeding_trees.root",
            "track_finder_tracks",
            368,
            "938bcc9b9425b12c620f5d0efa2c592817dfe92a18c309e97aa9d87412918620",
        ),
        (
            "performance_seeding_trees.root",
            "track_finder_particles",
            80,
            "938bcc9b9425b12c620f5d0efa2c592817dfe92a18c309e97aa9d87412918620",
        ),
        (
            "trackstates_ckf.root",
            "trackstates",
            368,
            "2faceafd4a521ff4030557301723e29c3d870edad052965eb644b824b57e2146",
        ),
        (
            "tracksummary_ckf.root",
            "tracksummary",
            10,
            "530ba8801c465780c8c66a788ae15592bc8ce36f8cfd8126216ea445f73c7e0d",
        ),
    ]

    csv = tmp_path / "csv"

    assert not csv.exists()
    for rf, _, _, _ in root_files:
        assert not (tmp_path / rf).exists()

    from ckf_tracks import runCKFTracks

    runCKFTracks(
        trackingGeometry,
        decorators,
        field=field,
        geometrySelection=Path(
            "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
        ),
        digiConfigFile=Path(
            "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
        ),
        outputCsv=True,
        outputDir=tmp_path,
        truthSmearedSeeded=False,
        truthEstimatedSeeded=False,
        s=s,
    )
    s.run()

    del s  # files are closed in destructors, not great

    assert csv.exists()
    for rf, tn, nume, file_hash in root_files:
        rp = tmp_path / rf
        assert rp.exists()
        if tn is not None and nume is not None:
            assert_entries(rp, tn, nume)
            assert_root_hash(rp, file_hash)

    assert len([f for f in csv.iterdir() if f.name.endswith("CKFtracks.csv")]) == events
    assert all([f.stat().st_size > 300 for f in csv.iterdir()])


def test_ckf_tracks_example_truth_estimate(tmp_path):
    # the example as written is only compatible with the generic detector
    detector, trackingGeometry, decorators = GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    events = 10
    s = Sequencer(events=events, numThreads=1)  # Digitization is not thread-safe

    root_files = [
        ("performance_ckf.root", None, None, None),
        (
            "performance_seeding_trees.root",
            "track_finder_tracks",
            80,
            "5c0cf9e84af64e6814ab1ddf4cbaf4be6008ad8b2371b5b0241085b19d0fc52c",
        ),
        (
            "performance_seeding_trees.root",
            "track_finder_particles",
            80,
            "5c0cf9e84af64e6814ab1ddf4cbaf4be6008ad8b2371b5b0241085b19d0fc52c",
        ),
        (
            "trackstates_ckf.root",
            "trackstates",
            80,
            "f44a901621438c92c6f20f214f5b7528ce5f9e2e72a6ea86fff594afb827afbc",
        ),
        (
            "tracksummary_ckf.root",
            "tracksummary",
            10,
            "cfbb83740cb8abc5180e4304d1dd552c424babb5e7f0caab65c6caa45f2a4313",
        ),
    ]

    csv = tmp_path / "csv"

    assert not csv.exists()
    for rf, _, _, _ in root_files:
        assert not (tmp_path / rf).exists()

    from ckf_tracks import runCKFTracks

    runCKFTracks(
        trackingGeometry,
        decorators,
        field=field,
        geometrySelection=Path(
            "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
        ),
        digiConfigFile=Path(
            "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
        ),
        outputCsv=True,
        outputDir=tmp_path,
        truthSmearedSeeded=False,
        truthEstimatedSeeded=True,
        s=s,
    )
    s.run()

    del s  # files are closed in destructors, not great

    assert csv.exists()
    for rf, tn, nume, file_hash in root_files:
        rp = tmp_path / rf
        assert rp.exists()
        if tn is not None and nume is not None:
            assert_entries(rp, tn, nume)
            assert_root_hash(rp, file_hash)

    assert len([f for f in csv.iterdir() if f.name.endswith("CKFtracks.csv")]) == events
    assert all([f.stat().st_size > 100 for f in csv.iterdir()])


def test_ckf_tracks_example_truth_smeared(tmp_path):
    # the example as written is only compatible with the generic detector
    detector, trackingGeometry, decorators = GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    events = 10
    s = Sequencer(events=events, numThreads=1)  # Digitization is not thread-safe

    root_files = [
        ("performance_ckf.root", None, None, None),
        (
            "trackstates_ckf.root",
            "trackstates",
            80,
            "bfa96b735050dfe77bb93631ce5a4e06e82a3ad8d65c321950d6d549f20f1ffa",
        ),
        (
            "tracksummary_ckf.root",
            "tracksummary",
            10,
            "70568044605d373b49bf941064d5047b1fc74816a390cd76674d13d61c6cb7aa",
        ),
    ]

    csv = tmp_path / "csv"

    assert not csv.exists()
    for rf, _, _, _ in root_files:
        assert not (tmp_path / rf).exists()

    from ckf_tracks import runCKFTracks

    runCKFTracks(
        trackingGeometry,
        decorators,
        field=field,
        geometrySelection=Path(
            "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
        ),
        digiConfigFile=Path(
            "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
        ),
        outputCsv=True,
        outputDir=tmp_path,
        truthSmearedSeeded=True,
        truthEstimatedSeeded=False,
        s=s,
    )
    s.run()

    del s  # files are closed in destructors, not great

    assert csv.exists()
    for rf, tn, nume, file_hash in root_files:
        rp = tmp_path / rf
        assert rp.exists()
        if tn is not None and nume is not None:
            assert_entries(rp, tn, nume)
            assert_root_hash(rp, file_hash)

    assert len([f for f in csv.iterdir() if f.name.endswith("CKFtracks.csv")]) == events
    assert all([f.stat().st_size > 300 for f in csv.iterdir()])
