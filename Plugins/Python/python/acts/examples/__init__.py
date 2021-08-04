from acts.ActsPythonBindings._examples import *
from acts import ActsPythonBindings
import acts
from acts._adapter import _patch_config, _patch_detectors

_propagators = []
_concrete_propagators = []
for prefix in ("Eigen", "Atlas", "StraightLine"):
    _propagators.append(getattr(ActsPythonBindings._propagator, f"{prefix}Propagator"))
    _concrete_propagators.append(
        getattr(ActsPythonBindings._propagator, f"{prefix}ConcretePropagator")
    )


def ConcretePropagator(propagator):
    for prop, prop_if in zip(_propagators, _concrete_propagators):
        if isinstance(propagator, prop):
            return prop_if(propagator)

    raise TypeError(f"Unknown propagator {type(propagator).__name__}")


_patch_config(ActsPythonBindings._examples)

_patch_detectors(ActsPythonBindings._examples)
