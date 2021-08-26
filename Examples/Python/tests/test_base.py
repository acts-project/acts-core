import pytest

import acts

import acts.examples


def test_logging():
    for l in ("VERBOSE", "DEBUG", "INFO", "WARNING", "ERROR", "FATAL"):
        assert hasattr(acts.logging, l)
        assert hasattr(acts.logging.Level, l)


def test_pgd_particle():
    assert len(acts.PdgParticle.__members__) == 16


def test_algebra():
    v3 = acts.Vector3(1, 2, 3)
    with pytest.raises(TypeError):
        acts.Vector3(1, 2, 3, 4)
    with pytest.raises(TypeError):
        acts.Vector3(1, 2)

    v3 = acts.Vector3([1, 2, 3])
    with pytest.raises(TypeError):
        acts.Vector3([1, 2, 3, 4])
    with pytest.raises(TypeError):
        acts.Vector3([1, 2])
    with pytest.raises(TypeError):
        acts.Vector3()

    v4 = acts.Vector4(1, 2, 3, 4)
    with pytest.raises(TypeError):
        v4 = acts.Vector4(1, 2, 3)
    v4 = acts.Vector4([1, 2, 3, 4])
    with pytest.raises(TypeError):
        acts.Vector4([1, 2, 3])
    with pytest.raises(TypeError):
        acts.Vector4()


def test_empty_sequencer(conf_const):
    s = acts.examples.Sequencer()
    with pytest.raises(RuntimeError):
        s.run()

    s = conf_const(acts.examples.Sequencer, events=1)
    s.run()


def test_random_number():
    rnd = acts.examples.RandomNumbers(seed=42)


def test_constructors():
    s1 = acts.examples.Sequencer()
    print(s1)
    s2 = acts.examples.Sequencer()
    print(s2)
