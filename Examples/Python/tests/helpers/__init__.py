import os
from typing import List, Union

from acts.examples import BareAlgorithm

geant4Enabled = any(v.startswith("G4") for v in os.environ.keys())

try:
    import ROOT

    rootEnabled = True
except ImportError:
    rootEnabled = False

    if "ROOTSYS" in os.environ:  # ROOT seems to be set up, but no PyROOT
        import warnings

        warnings.warn(
            "ROOT likely built without/with incompatible PyROOT. Skipping tests that need ROOT"
        )

dd4hepEnabled = "DD4hep_DIR" in os.environ

if dd4hepEnabled:
    try:
        import acts.examples.dd4hep
    except ImportError:
        dd4hepEnabled = False

try:
    import acts.examples.hepmc3

    hepmc3Enabled = True
except ImportError:
    hepmc3Enabled = False


class AssertCollectionExistsAlg(BareAlgorithm):
    events_seen = 0
    collections: List[str]

    def __init__(
        self,
        collections: Union[List[str], str],
        name="check_alg",
        level=acts.logging.INFO,
        *args,
        **kwargs,
    ):
        if isinstance(collections, str):
            self.collections = [collections]
        else:
            self.collections = collections
        BareAlgorithm.__init__(self, name=name, level=level, *args, **kwargs)

    def execute(self, ctx):
        for collection in self.collections:
            assert ctx.eventStore.exists(collection), f"{collection} does not exist"
        self.events_seen += 1
        return acts.examples.ProcessCode.SUCCESS
