#!/usr/bin/env python3

import tempfile
from pathlib import Path
import shutil

import acts
from truth_tracking_kalman_refitting import runRefittingKf

from physmon_common import makeSetup

setup = makeSetup()

with tempfile.TemporaryDirectory() as temp:
    s = acts.examples.Sequencer(
        events=10000,
        numThreads=-1,
        logLevel=acts.logging.INFO,
    )

    tp = Path(temp)
    runRefittingKf(
        trackingGeometry=setup.trackingGeometry,
        field=setup.field,
        digiConfigFile=setup.digiConfig,
        outputDir=tp,
        reverseFilteringMomThreshold=float("inf"),
        reverseFilteringCovarianceScaling=100.0,
        s=s,
    )

    s.run()

    perf_file = tp / "performance_kf_refit.root"
    assert perf_file.exists(), "Performance file not found"
    shutil.copy(perf_file, setup.outdir / "performance_trackrefitting.root")
