from pathlib import Path
from typing import Union

from .ActsPythonBindings import *
from .ActsPythonBindings import __version__
from . import ActsPythonBindings
from ._adapter import _patch_config


def Propagator(stepper, navigator):
    for prefix in ("Eigen", "Atlas", "StraightLine"):
        _stepper = getattr(ActsPythonBindings, f"{prefix}Stepper")
        if isinstance(stepper, _stepper):
            return getattr(ActsPythonBindings._propagator, f"{prefix}Propagator")(
                stepper, navigator
            )
    raise TypeError(f"Unknown stepper {type(stepper).__name__}")


_patch_config(ActsPythonBindings)


@staticmethod
def _decoratorFromFile(file: Union[str, Path], **kwargs):
    if isinstance(file, str):
        file = Path(file)

    kwargs.setdefault("level", ActsPythonBindings.logging.INFO)

    if file.suffix in (".json", ".cbor"):
        c = ActsPythonBindings.MaterialMapJsonConverter.Config()
        for k in kwargs.keys():
            if hasattr(c, k):
                setattr(c, k, kwargs.pop(k))

        return ActsPythonBindings.JsonMaterialDecorator(
            jFileName=str(file), rConfig=c, **kwargs
        )
    elif file.suffix == ".root":
        return ActsPythonBindings._examples.RootMaterialDecorator(
            fileName=str(file), **kwargs
        )
    else:
        raise ValueError(f"Unknown file type {file.suffix}")


ActsPythonBindings.IMaterialDecorator.fromFile = _decoratorFromFile
