from typing import Type
import pytest

import acts
import acts.examples

u = acts.UnitConstants


def test_null_bfield():
    assert acts.NullBField()


def test_constant_bfield():
    with pytest.raises(TypeError):
        acts.ConstantBField()
    v = acts.Vector3(1, 2, 3)
    cb = acts.ConstantBField(v)
    assert cb


def test_solenoid(conf_const):
    solenoid = conf_const(
        acts.SolenoidBField,
        radius=1200 * u.mm,
        length=6000 * u.mm,
        bMagCenter=2 * u.T,
        nCoils=1194,
    )

    field = acts.solenoidFieldMap(
        rlim=(0, 1200 * u.mm),
        zlim=(-5000 * u.mm, 5000 * u.mm),
        nbins=(10, 10),
        field=solenoid,
    )

    assert isinstance(field, acts.examples.InterpolatedMagneticField2)
