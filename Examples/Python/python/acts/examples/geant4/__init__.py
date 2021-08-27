import multiprocessing


# Cannot conveniently catch linker errors, so we launch a suprocess to
# try importing and see if it works in order to provide a useful error message
def _import_test():
    from acts import ActsPythonBindingsGeant4


p = multiprocessing.Process(target=_import_test)
p.start()
p.join()
if p.exitcode != 0:
    raise RuntimeError(
        "Error encountered importing Geant4. Likely you need to source $G4DIR/bin/geant4.sh."
    )

from acts._adapter import _patch_config
from acts import ActsPythonBindingsGeant4

_patch_config(ActsPythonBindingsGeant4)

from acts.ActsPythonBindingsGeant4 import *
