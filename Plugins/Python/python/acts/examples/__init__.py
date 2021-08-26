import inspect

from acts.ActsPythonBindings._examples import *
from acts import ActsPythonBindings
import acts
from acts._adapter import _patch_config, _patch_detectors, _patchKwargsConstructor

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


def _makeLayerTriplet(*args, **kwargs):

    if len(args) == 1:
        _type = type(args[0])
        negative = central = positive = args[0]
    else:
        negative = kwargs.get("negative")
        central = kwargs.get("central")
        positive = kwargs.get("positive")
        types = []
        if negative is not None:
            types.append(type(negative))
        if central is not None:
            types.append(type(central))
        if positive is not None:
            types.append(type(positive))

        assert all(types[0] == t for t in types), "Inconsistent types"

        _type = types[0]

    # print(_type, negative, central, positive)

    # type_map = {
    #     bool: LayerTripletBool,
    #     str: LayerTripletString,
    #     list: LayerTripletVectorString,
    #     tuple
    # }

    def fill(obj):
        if negative is not None:
            obj.negative = negative
        if central is not None:
            obj.central = central
        if positive is not None:
            obj.positive = positive
        return obj

    if _type == bool:
        return fill(TGeoDetector.Config.LayerTripletBool())
    elif _type == list:
        if all(
            all(isinstance(v, str) for v in vv)
            for vv in (negative, central, positive)
            if vv is not None
        ):
            return fill(TGeoDetector.Config.LayerTripletVectorString())
        else:
            raise TypeError("Invalid types for list input")
    elif _type == tuple:
        if all(
            all(isinstance(v, float) for v in vv)
            for vv in (negative, central, positive)
            if vv is not None
        ):
            negative = Interval(*negative) if negative is not None else None
            central = Interval(*central) if central is not None else None
            positive = Interval(*positive) if positive is not None else None
            return fill(TGeoDetector.Config.LayerTripletInterval())
        else:
            raise TypeError("Invalid types for tuple input")
    elif _type == Interval:
        return fill(TGeoDetector.Config.LayerTripletInterval())
    elif _type == str:
        return fill(TGeoDetector.Config.LayerTripletString())
    elif _type == float:
        return fill(TGeoDetector.Config.LayerTripletDouble())
    else:
        raise TypeError("Unknown type given")


TGeoDetector.Config.LayerTriplet = _makeLayerTriplet


def _process_volume_intervals(kwargs):
    if len(kwargs) == 0:
        return kwargs  # prevent infinite recursion
    _kwargs = kwargs.copy()

    v = TGeoDetector.Config.Volume()
    for name, value in inspect.getmembers(TGeoDetector.Config.Volume):
        if not isinstance(getattr(v, name), Interval):
            continue
        if not name in _kwargs:
            continue
        if not isinstance(_kwargs[name], tuple):
            continue
        _kwargs[name] = Interval(*_kwargs[name])

    return _kwargs


_patchKwargsConstructor(TGeoDetector.Config.Volume, proc=_process_volume_intervals)
