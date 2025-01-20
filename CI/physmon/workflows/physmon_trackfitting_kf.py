#!/usr/bin/env python3

import tempfile
from pathlib import Path
import shutil

import acts
from truth_tracking_kalman import runTruthTrackingKalman

from physmon_common import makeSetup

setup = makeSetup()

with tempfile.TemporaryDirectory() as temp:
    s = acts.examples.Sequencer(
        events=10000,
        numThreads=-1,
        logLevel=acts.logging.INFO,
    )

    tp = Path(temp)
    runTruthTrackingKalman(
        trackingGeometry=setup.trackingGeometry,
        field=setup.field,
        digiConfigFile=setup.digiConfig,
        outputDir=tp,
        reverseFilteringMomThreshold=float("inf"),
        reverseFilteringCovarianceScaling=100.0,
        s=s,
    )

    s.run()

    perf_file = tp / "performance_kf.root"
    assert perf_file.exists(), "Performance file not found"
    shutil.copy(perf_file, setup.outdir / "performance_trackfitting.root")
