import pytest

import acts
import acts.examples
from acts.examples.simulation import addParticleGun, EtaConfig, ParticleConfig


class AssertCollectionExistsAlg(acts.examples.IAlgorithm):
    events_seen = 0

    def __init__(self, collection, *args, **kwargs):
        self.collection = collection
        acts.examples.IAlgorithm.__init__(self, *args, **kwargs)

    def execute(self, ctx):
        assert ctx.eventStore.exists(self.collection)
        self.events_seen += 1
        return acts.examples.ProcessCode.SUCCESS


def test_navigator(conf_const):
    nav = conf_const(acts.Navigator)
    nav = conf_const(acts.Navigator, trackingGeometry=None)


def test_steppers(conf_const, trk_geo):
    with pytest.raises(TypeError):
        acts.examples.PropagationAlgorithm()
    with pytest.raises(ValueError):
        acts.examples.PropagationAlgorithm(level=acts.logging.INFO)

    with pytest.raises(TypeError):
        acts.Propagator()

    nav = acts.Navigator(trackingGeometry=trk_geo)

    with pytest.raises(TypeError):
        acts.Propagator(navigator=nav)

    for stepper in (acts.EigenStepper, acts.AtlasStepper):
        with pytest.raises(TypeError):
            stepper()
        s = stepper(acts.NullBField())
        assert s

        seq = acts.examples.Sequencer(
            events=10, numThreads=1, logLevel=acts.logging.WARNING
        )

        rnd = acts.examples.RandomNumbers(seed=42)

        addParticleGun(
            seq,
            ParticleConfig(num=10, pdg=acts.PdgParticle.eMuon, randomizeCharge=True),
            EtaConfig(-4.0, 4.0),
            rnd=rnd,
        )

        prop = acts.examples.ConcretePropagator(
            acts.Propagator(stepper=s, navigator=nav)
        )

        trkParamExtractor = acts.examples.ParticleTrackParamExtractor(
            level=acts.logging.WARNING,
            inputParticles="particles_generated",
            outputTrackParameters="params_particles_generated",
        )
        seq.addAlgorithm(trkParamExtractor)

        alg = conf_const(
            acts.examples.PropagationAlgorithm,
            level=acts.logging.WARNING,
            propagatorImpl=prop,
            inputTrackParameters="params_particles_generated",
            outputSummaryCollection="propagation_summary",
            sterileLogger=False,
        )
        seq.addAlgorithm(alg)

        chkAlg = AssertCollectionExistsAlg(
            "propagation_summary", "chk_alg", level=acts.logging.WARNING
        )
        seq.addAlgorithm(chkAlg)

        seq.run()

    assert acts.StraightLineStepper()
