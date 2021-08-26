from pathlib import Path
import sys
import os
import tempfile


sys.path += [
    str(Path(__file__).parent.parent.parent.parent / "Examples/Scripts/Python/"),
    str(Path(__file__).parent),
]

import helpers

import pytest

import acts
import acts.examples

u = acts.UnitConstants


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

        nav = acts.Navigator(trackingGeometry=geo)
        # field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))
        # stepper = acts.EigenStepper(field)
        stepper = acts.StraightLineStepper()

        prop = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))
        # rnd = acts.examples.RandomNumbers(seed=42)
        alg = acts.examples.PropagationAlgorithm(
            propagatorImpl=prop,
            level=acts.logging.INFO,
            randomNumberSvc=rng,
            ntests=10,
            sterileLogger=False,
            propagationStepCollection="propagation-steps",
        )
        s.addAlgorithm(alg)
        return s, alg

    return _basic_prop_seq_factory


# def pytest_generate_tests(metafunc):
#     if "trk_geo" in metafunc.fixturenames:
#         geometries = ["generic"]

#         if helpers.dd4hepEnabled:
#             geometries.append("dd4hep")

#         metafunc.parametrize("trk_geo", geometries, indirect=True)


@pytest.fixture
def trk_geo(request):
    detector, geo, contextDecorators = acts.examples.GenericDetector.create()
    yield geo

    # if request.param == "generic":
    #     detector, geo, contextDecorators = acts.examples.GenericDetector.create()
    #     request._detector = detector
    #     return geo
    # elif request.param == "dd4hep":
    #     from acts.examples.dd4hep import DD4hepGeometryService, DD4hepDetector

    #     dd4hepConfig = DD4hepGeometryService.Config(
    #         xmlFileNames=["thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
    #     )
    #     detector = DD4hepDetector()

    #     config = acts.MaterialMapJsonConverter.Config()
    #     matDeco = acts.JsonMaterialDecorator(
    #         rConfig=config,
    #         jFileName="thirdparty/OpenDataDetector/config/odd-material-mapping.config",
    #         level=acts.logging.ERROR,
    #     )

    #     geo, _ = detector.finalize(dd4hepConfig, matDeco)
    #     request._detector = detector
    #     return geo
    # else:
    #     raise ValueError("Invalid detector parameter")


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
            outputParticles="particles_input",
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
        simAlg = acts.examples.FatrasAlgorithm(
            level=acts.logging.INFO,
            inputParticles=evGen.config.outputParticles,
            outputParticlesInitial="particles_initial",
            outputParticlesFinal="particles_final",
            outputSimHits="simhits",
            randomNumbers=rng,
            trackingGeometry=trk_geo,
            magneticField=field,
            generateHitsOnSensitive=True,
        )

        s.addAlgorithm(simAlg)

        # Digitization
        digiCfg = acts.examples.DigitizationConfig(
            acts.examples.readDigiConfigFromJson(
                "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
            ),
            trackingGeometry=trk_geo,
            randomNumbers=rng,
            inputSimHits=simAlg.config.outputSimHits,
        )
        # digiCfg.trackingGeometry = trackingGeometry
        # digiCfg.randomNumbers = rnd
        # digiCfg.inputSimHits = simAlg.config.outputSimHits
        digiAlg = acts.examples.DigitizationAlgorithm(digiCfg, acts.logging.INFO)

        s.addAlgorithm(digiAlg)

        return evGen, simAlg, digiAlg

    return _factory


@pytest.fixture(scope="session")
def geantino_recording():
    if not helpers.geant4Enabled:
        pytest.skip("Geantino recording requested, but Geant4 is not set up")

    if not helpers.dd4hepEnabled:
        pytest.skip("DD4hep recording requested, but Geant4 is not set up")

    from geantino_recording import runGeantinoRecording

    dd4hepSvc = acts.examples.dd4hep.DD4hepGeometryService(
        xmlFileNames=["thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
    )
    dd4hepG4ConstructionFactory = (
        acts.examples.geant4.dd4hep.DD4hepDetectorConstructionFactory(dd4hepSvc)
    )

    with tempfile.TemporaryDirectory() as d:

        s = acts.examples.Sequencer(events=2, numThreads=1)

        runGeantinoRecording(
            dd4hepG4ConstructionFactory, str(d), tracksPerEvent=100, s=s
        )
        s.run()

        del s
        del dd4hepSvc
        del dd4hepG4ConstructionFactory

        yield Path(d)
