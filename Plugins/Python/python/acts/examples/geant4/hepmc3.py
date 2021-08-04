from acts._adapter import _patch_config
from acts.ActsPythonBindingsGeant4 import hepmc3 as _hepmc3

_patch_config(_hepmc3)

from acts.ActsPythonBindingsGeant4.hepmc3 import *
