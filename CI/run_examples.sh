#!/bin/bash
# Run a reasonably large subset of the Acts examples

# Force bash to catch more errors
set -euo pipefail

# Common setup
NUM_EVENTS=100
SRC_DIR=`pwd`
BUILD_DIR=`pwd`/build
DD4HEP_INPUT="--dd4hep-input file:${SRC_DIR}/thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"
timed_run() {
    echo ""
    echo "=== Running $* ==="
    echo ""
    time ${BUILD_DIR}/bin/$*
}
run_example() {
    timed_run $* -n ${NUM_EVENTS}
}

# Data for Geant4-based examples
G4_DATA_URL=https://acts.web.cern.ch/ACTS/ci/G4/G4.data.tar.xz
G4_DATA_DIR=/tmp/G4
#
curl $G4_DATA_URL \
    -O \
    -C - \
    --retry 5
#
mkdir -p $G4_DATA_DIR
tar xf G4.data.tar.xz --strip-components 1 -C $G4_DATA_DIR
rm G4.data.tar.xz
#
export G4ENSDFSTATEDATA=$G4_DATA_DIR/G4ENSDFSTATE
export G4LEVELGAMMADATA=$G4_DATA_DIR/G4LEVELGAMMA
export G4LEDATA=$G4_DATA_DIR/G4LE
export G4PARTICLEXSDATA=$G4_DATA_DIR/G4PARTICLESXS

# Run hello world example
run_example ActsExampleHelloWorld

# Run geometry examples for all geometries, in the configuration suggested by
# the material mapping tutorial
#
# We must try to avoid running examples for all geometries because
# that results in a combinatorial explosion of CI running time. But
# these examples are fast enough.
#
run_geometry_example() {
    timed_run ActsExampleGeometry$* \
                  -n1 \
                  -j1 \
                  --mat-output-file geometry-map-$1 \
                  --output-json \
                  --mat-output-allmaterial true \
                  --mat-output-sensitives false
}
run_geometry_example Aligned
run_geometry_example DD4hep ${DD4HEP_INPUT}
run_geometry_example Empty
run_geometry_example Generic
run_geometry_example Payload
run_geometry_example Telescope
# TODO: Add TGeo geometry example (needs an input file + knowhow)

# Run propagation examples for all geometries
#
# In addition to the geometry examples, we also need one slightly more
# complex example to run with all geometries because geometries with
# conditions have code paths that are only exercised when the
# Sequencer actually runs.
#
run_example ActsExamplePropagationAligned
run_example ActsExamplePropagationDD4hep ${DD4HEP_INPUT}
# FIXME: Disabled because of issue #710
# run_example ActsExamplePropagationEmpty
run_example ActsExamplePropagationGeneric
run_example ActsExamplePropagationPayload
# TODO: Add TGeo propagation example (needs an input file + knowhow)

# Run event generation examples as suggested by the Fatras tutorial
run_example ActsExampleParticleGun \
                --output-dir=data/gen/four_muons \
                --output-csv \
                --gen-phi-degree=0:90 \
                --gen-eta=-2:2 \
                --gen-mom-gev=1:5 \
                --gen-pdg=13 \
                --gen-randomize-charge \
                --gen-nparticles=4
run_example ActsExamplePythia8 \
                --output-dir=data/gen/ttbar_mu140 \
                --output-csv \
                --rnd-seed=42 \
                --gen-cms-energy-gev=14000 \
                --gen-hard-process=Top:qqbar2ttbar=on \
                --gen-npileup=140

# Run Geantino recording example, as in the material mapping tutorial (but with
# 10x less events to keep CI running times reasonable)
#
# FIXME: Currently only works in single-threaded mode, even though it
#        should theoretically be thread-safe as Geant4 usage is
#        protected by a mutex. What most likely happens is that some
#        thread-unsafe Geant4 code is accidentally run outside of the
#        mutex-protected region of the code. See issue #207 .
#
run_geantino_example() {
    timed_run ActsExampleMaterialRecording$* \
                  -n1000 \
                  -j1 \
                  --g4-material-tracks=geant4-material-tracks-$1 \
                  --output-root
}
run_geantino_example DD4hep ${DD4HEP_INPUT}
# TODO: Add GDML version (needs an input file + knowhow)

# Run material validation example (generic-only, see above)
run_example ActsExampleMaterialValidationGeneric

# TODO: Run CKF examples (needs some setup, documented in Acts howto?)
# TODO: Run EventRecording example (requires CSV input, source to be found. Fatras?)
# TODO: Run Fatras examples (needs some setup, documented in Acts howto?)
# TODO: Run HepMC3 examples (needs input files + knowhow)
# TODO: Run MagneticField examples (needs input files + knowhow)
# TODO: Run MaterialMapping examples (needs input files + knowhow)
# TODO: Run Seeding examples (needs input files + knowhow)
# TODO: Run Show(FatrasGeneric|Particles) examples (needs input files + knowhow)
# TODO: Run TruthTracks examples (needs input files + knowhow)
# TODO: Run (AdaptiveMulti|Iterative)?Vertex(Finder|Fitter) examples (currently broken, not ready for CI)

# TODO: Bring back multi-threaded output reproducibility tests

# Run vertex finder tutorial
run_example ActsTutorialVertexFinder
