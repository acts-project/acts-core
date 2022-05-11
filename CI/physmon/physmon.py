#!/usr/bin/env python3
from pathlib import Path
import argparse
import tempfile
import shutil
import os
import sys
import subprocess

sys.path += [
    str(Path(__file__).parent.parent.parent / "Examples/Scripts/Python/"),
]

# this has to happen before we import the ACTS module
os.environ["ACTS_LOG_FAILURE_THRESHOLD"] = "FATAL"
import acts.examples

from particle_gun import addParticleGun, MomentumConfig, EtaConfig, ParticleConfig
from fatras import addFatras
from truth_tracking import runTruthTracking
from digitization import addDigitization
from ckf_tracks import runCKFTracks
from seeding import addSeeding, SeedingAlgorithm, TruthSeedRanges
from common import getOpenDataDetector

parser = argparse.ArgumentParser()
parser.add_argument("outdir")
parser.add_argument("--events", type=int, default=10000)
parser.add_argument("--skip", type=int, default=0)

args = parser.parse_args()

outdir = Path(args.outdir)
outdir.mkdir(exist_ok=True)


srcdir = Path(__file__).resolve().parent.parent.parent


u = acts.UnitConstants

matDeco = acts.IMaterialDecorator.fromFile(
    srcdir / "thirdparty/OpenDataDetector/data/odd-material-maps.root",
    level=acts.logging.INFO,
)
detector, trackingGeometry, decorators = getOpenDataDetector(matDeco)
digiConfig = srcdir / "thirdparty/OpenDataDetector/config/odd-digi-smearing-config.json"
geoSel = srcdir / "thirdparty/OpenDataDetector/config/odd-seeding-config.json"


field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

s = acts.examples.Sequencer(
    events=args.events, numThreads=-1, logLevel=acts.logging.INFO, skip=0
)

with tempfile.TemporaryDirectory() as temp:
    tp = Path(temp)
    runTruthTracking(
        trackingGeometry,
        field,
        digiConfigFile=digiConfig,
        outputDir=tp,
        s=s,
    )

    s.run()
    del s

    perf_file = tp / "performance_track_fitter.root"
    assert perf_file.exists(), "Performance file not found"
    shutil.copy(perf_file, outdir / "performance_truth_tracking.root")


for truthSmeared, truthEstimated, label in [
    (True, False, "truth_smeared"),  # if first is true, second is ignored
    (False, True, "truth_estimated"),
    (False, False, "seeded"),
]:
    s = acts.examples.Sequencer(
        events=args.events, numThreads=1, logLevel=acts.logging.INFO, skip=args.skip
    )

    with tempfile.TemporaryDirectory() as temp:
        tp = Path(temp)
        runCKFTracks(
            trackingGeometry,
            decorators=decorators,
            field=field,
            digiConfigFile=digiConfig,
            geometrySelection=geoSel,
            outputDir=tp,
            outputCsv=False,
            truthSmearedSeeded=truthSmeared,
            truthEstimatedSeeded=truthEstimated,
            s=s,
        )

        s.run()
        del s

        perf_file = tp / "performance_ckf.root"
        assert perf_file.exists(), "Performance file not found"
        shutil.copy(perf_file, outdir / f"performance_ckf_tracks_{label}.root")

        if not truthSmeared and not truthEstimated:
            residual_app = srcdir / "build/bin/ActsAnalysisResidualsAndPulls"
            # @TODO: Add try/except
            subprocess.check_call(
                [
                    str(residual_app),
                    "--predicted",
                    "--filtered",
                    "--smoothed",
                    "--silent",
                    "-i",
                    str(tp / "trackstates_ckf.root"),
                    "-o",
                    str(outdir / "acts_analysis_residuals_and_pulls.root"),
                ]
            )
