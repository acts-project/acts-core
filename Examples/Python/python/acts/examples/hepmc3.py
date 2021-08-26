from acts._adapter import _patch_config
from acts import ActsPythonBindings

_patch_config(ActsPythonBindings._examples._hepmc3)

from acts.ActsPythonBindings._examples._hepmc3 import *
