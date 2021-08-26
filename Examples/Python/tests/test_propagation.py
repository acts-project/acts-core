import pytest

import acts
import acts.examples


class AssertCollectionExistsAlg(acts.examples.BareAlgorithm):
    events_seen = 0

    def __init__(self, collection, *args, **kwargs):
        self.collection = collection
        acts.examples.BareAlgorithm.__init__(self, *args, **kwargs)

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

        prop = acts.examples.ConcretePropagator(
            acts.Propagator(stepper=s, navigator=nav)
        )

        # cfg = acts.examples.PropagationAlgorithm.Config()
        # cfg.propagatorImpl = prop
        # acts.examples.PropagationAlgorithm(level=acts.logging.INFO, config=cfg)
        # acts.examples.PropagationAlgorithm(cfg, acts.logging.INFO)

        alg = conf_const(
            acts.examples.PropagationAlgorithm,
            level=acts.logging.ERROR,
            propagatorImpl=prop,
            randomNumberSvc=acts.examples.RandomNumbers(),
            propagationStepCollection="propagation_steps",
            sterileLogger=False,
            ntests=10,
        )

        seq = acts.examples.Sequencer(
            events=10, numThreads=1, logLevel=acts.logging.ERROR
        )
        seq.addAlgorithm(alg)
        chkAlg = AssertCollectionExistsAlg(
            "propagation_steps", "chk_alg", level=acts.logging.ERROR
        )
        seq.addAlgorithm(chkAlg)
        seq.run()

    assert acts.StraightLineStepper()
