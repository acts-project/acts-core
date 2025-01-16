import multiprocessing
from pathlib import Path
import sys
import os
import tempfile
import shutil
from typing import Dict
import warnings
import pytest_check as check
from collections import namedtuple


sys.path += [
    str(Path(__file__).parent.parent.parent.parent / "Examples/Scripts/Python/"),
    str(Path(__file__).parent),
]


import helpers
import helpers.hash_root

import pytest

import acts
import acts.examples
from acts.examples.odd import getOpenDataDetector
from acts.examples.simulation import addParticleGun, EtaConfig, ParticleConfig

try:
    import ROOT

    ROOT.gSystem.ResetSignals()
except ImportError:
    pass

try:
    if acts.logging.getFailureThreshold() != acts.logging.WARNING:
        acts.logging.setFailureThreshold(acts.logging.WARNING)
except RuntimeError:
    # Repackage with different error string
    errtype = (
        "negative"
        if acts.logging.getFailureThreshold() < acts.logging.WARNING
        else "positive"
    )
    warnings.warn(
        "Runtime log failure threshold could not be set. "
        "Compile-time value is probably set via CMake, i.e. "
        f"`ACTS_LOG_FAILURE_THRESHOLD={acts.logging.getFailureThreshold().name}` is set, "
        "or `ACTS_ENABLE_LOG_FAILURE_THRESHOLD=OFF`. "
        f"The pytest test-suite can produce false-{errtype} results in this configuration"
    )


u = acts.UnitConstants


class RootHashAssertionError(AssertionError):
    def __init__(
        self, file: Path, key: str, exp_hash: str, act_hash: str, *args, **kwargs
    ):
        super().__init__(f"{exp_hash} != {act_hash}", *args, **kwargs)
        self.file = file
        self.key = key
        self.exp_hash = exp_hash
        self.act_hash = act_hash


hash_assertion_failures = []


def _parse_hash_file(file: Path) -> Dict[str, str]:
    res = {}
    for line in file.open():
        if line.strip() == "" or line.strip().startswith("#"):
            continue
        key, h = line.strip().split(":", 1)
        res[key.strip()] = h.strip()
    return res


@pytest.fixture(scope="session")
def root_file_exp_hashes():
    path = Path(
        os.environ.get("ROOT_HASH_FILE", Path(__file__).parent / "root_file_hashes.txt")
    )
    return _parse_hash_file(path)


@pytest.fixture(name="assert_root_hash")
def assert_root_hash(request, root_file_exp_hashes):
    if not helpers.doHashChecks:

        def fn(*args, **kwargs):
            pass

        return fn

    def fn(key: str, file: Path):
        """
        Assertion helper function to check the hashes of root files.
        Do NOT use this function directly by importing, rather use it as a pytest fixture

        Arguments you need to provide:
        key: Explicit lookup key for the expected hash, should be unique per test function
        file: Root file to check the expected hash against
        """
        __tracebackhide__ = True
        gkey = f"{request.node.name}__{key}"
        act_hash = helpers.hash_root.hash_root_file(file)
        if not gkey in root_file_exp_hashes:
            warnings.warn(
                f'Hash lookup key "{key}" not found for test "{request.node.name}"'
            )
            check.equal(act_hash, "[MISSING]")
            exc = RootHashAssertionError(file, gkey, "[MISSING]", act_hash)
            hash_assertion_failures.append(exc)

        else:
            refhash = root_file_exp_hashes[gkey]
            check.equal(act_hash, refhash)
            if act_hash != refhash:
                exc = RootHashAssertionError(file, gkey, refhash, act_hash)
                hash_assertion_failures.append(exc)

    return fn


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    docs_url = "https://acts.readthedocs.io/en/latest/examples/python_bindings.html#root-file-hash-regression-checks"
    if len(hash_assertion_failures) > 0:
        terminalreporter.ensure_newline()
        terminalreporter.section(
            "RootHashAssertionErrors", sep="-", red=True, bold=True
        )
        terminalreporter.line(
            "The ROOT files produced by tests have changed since the last recorded reference."
        )
        terminalreporter.line(
            "This can be be expected if e.g. the underlying algorithm changed, or it can be a test failure symptom."
        )
        terminalreporter.line(
            "Please manually check the output files listed below and make sure that their content is correct."
        )
        terminalreporter.line(
            "If it is, you can update the test reference file Examples/Python/tests/root_file_hashes.txt with the new hashes below."
        )
        terminalreporter.line(f"See {docs_url} for more details")
        terminalreporter.line("")

        for e in hash_assertion_failures:
            terminalreporter.line(f"{e.key}: {e.act_hash}")

    if not helpers.doHashChecks:
        terminalreporter.section("Root file has checks", sep="-", blue=True, bold=True)
        terminalreporter.line(
            "NOTE: Root file hash checks were skipped, enable with ROOT_HASH_CHECKS=on"
        )
        terminalreporter.line(f"See {docs_url} for more details")


def kwargsConstructor(cls, *args, **kwargs):
    return cls(*args, **kwargs)


def configKwConstructor(cls, *args, **kwargs):
    assert hasattr(cls, "Config")
    _kwargs = {}
    if "level" in kwargs:
        _kwargs["level"] = kwargs.pop("level")
    config = cls.Config()
    for k, v in kwargs.items():
        setattr(config, k, v)
    return cls(*args, config=config, **_kwargs)


def configPosConstructor(cls, *args, **kwargs):
    assert hasattr(cls, "Config")
    _kwargs = {}
    if "level" in kwargs:
        _kwargs["level"] = kwargs.pop("level")
    config = cls.Config()
    for k, v in kwargs.items():
        setattr(config, k, v)

    return cls(config, *args, **_kwargs)


@pytest.fixture(params=[configPosConstructor, configKwConstructor, kwargsConstructor])
def conf_const(request):
    return request.param


@pytest.fixture
def rng():
    return acts.examples.RandomNumbers(seed=42)


@pytest.fixture
def basic_prop_seq(rng):
    def _basic_prop_seq_factory(geo, s=None):
        if s is None:
            s = acts.examples.Sequencer(events=10, numThreads=1)

        addParticleGun(
            s,
            ParticleConfig(num=10, pdg=acts.PdgParticle.eMuon, randomizeCharge=True),
            EtaConfig(-4.0, 4.0),
            rnd=rng,
        )

        trkParamExtractor = acts.examples.ParticleTrackParamExtractor(
            level=acts.logging.WARNING,
            inputParticles="particles_generated",
            outputTrackParameters="params_particles_generated",
        )
        s.addAlgorithm(trkParamExtractor)

        nav = acts.Navigator(trackingGeometry=geo)
        stepper = acts.StraightLineStepper()

        prop = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))

        alg = acts.examples.PropagationAlgorithm(
            level=acts.logging.WARNING,
            propagatorImpl=prop,
            sterileLogger=False,
            inputTrackParameters="params_particles_generated",
            outputSummaryCollection="propagation_summary",
        )
        s.addAlgorithm(alg)

        return s, alg

    return _basic_prop_seq_factory


@pytest.fixture
def trk_geo():
    detector = acts.examples.GenericDetector()
    trackingGeometry = detector.trackingGeometry()
    yield trackingGeometry


DetectorConfig = namedtuple(
    "DetectorConfig",
    [
        "detector",
        "trackingGeometry",
        "decorators",
        "geometrySelection",
        "digiConfigFile",
        "name",
    ],
)


@pytest.fixture(params=["generic", pytest.param("odd", marks=pytest.mark.odd)])
def detector_config(request):
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    if request.param == "generic":
        detector = acts.examples.GenericDetector()
        trackingGeometry = detector.trackingGeometry()
        decorators = detector.contextDecorators()
        return DetectorConfig(
            detector,
            trackingGeometry,
            decorators,
            geometrySelection=(
                srcdir
                / "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
            ),
            digiConfigFile=(
                srcdir
                / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
            ),
            name=request.param,
        )
    elif request.param == "odd":
        if not helpers.dd4hepEnabled:
            pytest.skip("DD4hep not set up")

        matDeco = acts.IMaterialDecorator.fromFile(
            srcdir / "thirdparty/OpenDataDetector/data/odd-material-maps.root",
            level=acts.logging.INFO,
        )
        detector = getOpenDataDetector(matDeco)
        trackingGeometry = detector.trackingGeometry()
        decorators = detector.contextDecorators()
        return DetectorConfig(
            detector,
            trackingGeometry,
            decorators,
            digiConfigFile=(
                srcdir
                / "thirdparty/OpenDataDetector/config/odd-digi-smearing-config.json"
            ),
            geometrySelection=(
                srcdir / "thirdparty/OpenDataDetector/config/odd-seeding-config.json"
            ),
            name=request.param,
        )
    else:
        raise ValueError(f"Invalid detector {detector}")


@pytest.fixture
def ptcl_gun(rng):
    def _factory(s):
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
                        phi=(0, 360 * u.degree),
                        randomizeCharge=True,
                        numParticles=2,
                    ),
                )
            ],
            outputParticles="particles_generated",
            outputVertices="vertices_input",
            randomNumbers=rng,
        )

        s.addReader(evGen)

        return evGen

    return _factory


@pytest.fixture
def fatras(ptcl_gun, trk_geo, rng):
    def _factory(s):
        evGen = ptcl_gun(s)

        field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))
        simAlg = acts.examples.FatrasSimulation(
            level=acts.logging.INFO,
            inputParticles=evGen.config.outputParticles,
            outputParticles="particles_simulated",
            outputSimHits="simhits",
            randomNumbers=rng,
            trackingGeometry=trk_geo,
            magneticField=field,
            generateHitsOnSensitive=True,
            emScattering=False,
            emEnergyLossIonisation=False,
            emEnergyLossRadiation=False,
            emPhotonConversion=False,
        )

        s.addAlgorithm(simAlg)

        # Digitization
        digiCfg = acts.examples.DigitizationAlgorithm.Config(
            digitizationConfigs=acts.examples.readDigiConfigFromJson(
                str(
                    Path(__file__).parent.parent.parent.parent
                    / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
                )
            ),
            surfaceByIdentifier=trk_geo.geoIdSurfaceMap(),
            randomNumbers=rng,
            inputSimHits=simAlg.config.outputSimHits,
        )
        digiAlg = acts.examples.DigitizationAlgorithm(digiCfg, acts.logging.INFO)

        s.addAlgorithm(digiAlg)

        return evGen, simAlg, digiAlg

    return _factory


def _do_material_recording(d: Path):
    from material_recording import runMaterialRecording

    s = acts.examples.Sequencer(events=2, numThreads=1)

    with getOpenDataDetector() as detector:
        runMaterialRecording(detector, str(d), tracksPerEvent=100, s=s)

        s.run()


@pytest.fixture(scope="session")
def material_recording_session():
    if not helpers.geant4Enabled:
        pytest.skip("Geantino recording requested, but Geant4 is not set up")

    if not helpers.dd4hepEnabled:
        pytest.skip("DD4hep recording requested, but DD4hep is not set up")

    with tempfile.TemporaryDirectory() as d:
        # explicitly ask for "spawn" as CI failures were observed with "fork"
        spawn_context = multiprocessing.get_context("spawn")
        p = spawn_context.Process(target=_do_material_recording, args=(d,))
        p.start()
        p.join()
        if p.exitcode != 0:
            raise RuntimeError("Failure to exeecute material recording")

        yield Path(d)


@pytest.fixture
def material_recording(material_recording_session: Path, tmp_path: Path):
    target = tmp_path / material_recording_session.name
    shutil.copytree(material_recording_session, target)
    yield target
