#!/bin/bash

set -e

# helper function to selectively print and run commands without a subshell
function run() {
    set -x
    "$@"
    # save exit code
    { rec=$?; } 2> /dev/null
    { set +x;   } 2> /dev/null
    # restore exit code
    (exit $rec)
}

export run

run which python3
shopt -s extglob


mode=${1:-all}
if ! [[ $mode = @(all|kf|gsf|gx2f|refit_kf|refit_gsf|fullchains|simulation|gx2f_vs_kf) ]]; then
    echo "Usage: $0 <all|kf|gsf|gx2f|refit_kf|refit_gsf|fullchains|simulation|gx2f_vs_kf> (outdir)"
    exit 1
fi

outdir=${2:-physmon}
mkdir -p $outdir
mkdir -p $outdir/data
mkdir -p $outdir/html
mkdir -p $outdir/logs

refdir=CI/physmon/reference
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}"  )" &> /dev/null && pwd  )

# File to accumulate the histcmp results
histcmp_results=$outdir/histcmp_results.csv
echo -n "" > $histcmp_results

memory_dir=${outdir}/memory
mkdir -p "$memory_dir"

if [ "$(uname)" == "Darwin" ]; then
    function measure {
        label=$1
        slug=$2
        shift
        shift
        echo "Measure Darwin $label ($slug)" >&2
        tmp=$(mktemp)
        echo "+ $@" >&2
        /usr/bin/time -l -o "$tmp" "$@"
        # save exit code
        rec=$?

        of="${memory_dir}/mem_${slug}.csv"
        {
            echo "# spyral-label: $label"
            echo "# spyral-cmd: $*"
            echo "time,rss,vms"
            # in bytes
            grep -E "real" "$tmp" | awk '{printf $1}'
            printf ","
            grep -E "maximum resident set size" "$tmp" | awk '{printf $1}'
            printf ",0\n"
        } > "$of"
        # restore exit code
        (exit $rec)
    }
    export measure
elif [ "$(uname)" == "Linux" ]; then
    function measure {
        label=$1
        slug=$2
        shift
        shift
        echo "Measure Linux $label ($slug)" >&2
        tmp=$(mktemp)
        echo "+ $@" >&2
        /usr/bin/time -v -o  "$tmp" "$@"
        # save exit code
        rec=$?
        # in kbytes
        max_rss=$(grep "Maximum resident set size (kbytes):" "$tmp" | awk '{printf $(NF)}')
        max_rss=$(( 1024 * max_rss ))
        echo "Maximum resident set size: $(printf "%'d" $max_rss) bytes"

        wall_time=$(grep "Elapsed (wall clock)" "$tmp" | awk '{printf $(NF)}')
        wall_time=$(python3 -c "i='${wall_time}';p=i.split(':');p = p if len(p) == 3 else ['0', *p];t=float(p[0])*60*60 + float(p[1])*60 + float(p[2]);print(t)")
        echo "Elapsed (wall clock) time: ${wall_time} seconds"

        of="${memory_dir}/mem_${slug}.csv"
        {
            echo "# spyral-label: $label"
            echo "# spyral-cmd: $*"
            echo "time,rss,vms"
            echo "${wall_time},${max_rss},0"
        } > "$of"
        # restore exit code
        (exit $rec)
    }
    export measure
else
    function measure {
        echo "Not measuring because unknown environment"
        shift
        shift
        "$@"
    }
    export measure
fi

if [ -n "$CI" ]; then
    echo "CI mode, do not abort immediately on failure"
    set +e
fi
ec=0

source $SCRIPT_DIR/setup.sh

function run_physmon_gen() {
    title=$1
    slug=$2

    script=CI/physmon/workflows/physmon_${slug}.py

    mkdir -p $outdir/data/$slug
    mkdir -p $outdir/logs
    measure "$title" "$slug" ${script} $outdir/data/$slug 2>&1 > $outdir/logs/${slug}.log

    this_ec=$?
    ec=$(($ec | $this_ec))

    if [ $this_ec -ne 0 ]; then
        echo "::error::🟥 Dataset generation failed: ${script} -> ec=$this_ec"
    else
        echo "::notice::✅ Dataset generation succeeded: ${script}"
    fi
}

echo "::group::Generate validation dataset"
if [[ "$mode" == "all" || "$mode" == "simulation" ]]; then
    run_physmon_gen "Simulation" "simulation"
fi
if [[ "$mode" == "all" || "$mode" == "kf" ]]; then
    run_physmon_gen "Truth Tracking KF" "trackfitting_kf"
fi
if [[ "$mode" == "all" || "$mode" == "gsf" ]]; then
    run_physmon_gen "Truth Tracking GSF" "trackfitting_gsf"
fi
if [[ "$mode" == "all" || "$mode" == "gx2f" ]]; then
    run_physmon_gen "Truth Tracking GX2F" "trackfitting_gx2f"
fi
if [[ "$mode" == "all" || "$mode" == "refit_kf" ]]; then
    run_physmon_gen "Truth Tracking KF refit" "trackrefitting_kf"
fi
if [[ "$mode" == "all" || "$mode" == "refit_gsf" ]]; then
    run_physmon_gen "Truth Tracking GSF refit" "trackrefitting_gsf"
fi
if [[ "$mode" == "all" || "$mode" == "fullchains" ]]; then
    run_physmon_gen "CKF single muon" "trackfinding_1muon"
    run_physmon_gen "CKF muon 50" "trackfinding_4muon_50vertices"
    run_physmon_gen "CKF ttbar 200" "trackfinding_ttbar_pu200"
fi
if [[ "$mode" == "all" || "$mode" == "gx2f_vs_kf" ]]; then
    run_physmon_gen "Comparison - Truth Tracking GX2F vs KF" "trackfitting_gx2f_vs_kf"
fi
echo "::endgroup::"


function run_histcmp() {
    a=$1
    b=$2
    title=$3
    html_path=$4
    plots_path=$5
    shift 5

    echo "::group::Comparing $a vs. $b"

    if [ ! -f "$a" ]; then
        echo "::error::histcmp failed: File $a does not exist"
        ec=1
    fi

    if [ ! -f "$b" ]; then
        echo "::error::histcmp failed: File $b does not exist"
        ec=1
    fi

    run histcmp $a $b \
        --label-reference=reference \
        --label-monitored=monitored \
        --title="$title" \
        -o $outdir/html/$html_path \
        -p $outdir/html/$plots_path \
        "$@"

    this_ec=$?
    ec=$(($ec | $this_ec))

    if [ $this_ec -ne 0 ]; then
        echo "::error::histcmp failed: ec=$this_ec"
    fi

    echo "\"${title}\",html/${html_path},${this_ec}" >> $histcmp_results

    echo "::endgroup::"
}

function trackfinding() {
    name=$1
    path=$2

    default_config="CI/physmon/config/default.yml"

    if [ -f $refdir/$path/performance_seeding.root ]; then
        run_histcmp \
            $outdir/data/$path/performance_seeding.root \
            $refdir/$path/performance_seeding.root \
            "Seeding ${name}" \
            $path/performance_seeding.html \
            $path/performance_seeding_plots \
            --config $default_config
    fi

    run_histcmp \
        $outdir/data/$path/performance_finding_ckf.root \
        $refdir/$path/performance_finding_ckf.root \
        "CKF finding performance | ${name}" \
        $path/performance_finding_ckf.html \
        $path/performance_finding_ckf_plots \
        --config $default_config

    run_histcmp \
        $outdir/data/$path/performance_fitting_ckf.root \
        $refdir/$path/performance_fitting_ckf.root \
        "CKF fitting performance | ${name}" \
        $path/performance_fitting_ckf.html \
        $path/performance_fitting_ckf_plots \
        --config $default_config


    run Examples/Scripts/generic_plotter.py \
        $outdir/data/$path/tracksummary_ckf.root \
        tracksummary \
        $outdir/data/$path/tracksummary_ckf_hist.root \
        --silent \
        --config CI/physmon/config/tracksummary_ckf.yml
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/data/$path/tracksummary_ckf.root

    run_histcmp \
        $outdir/data/$path/tracksummary_ckf_hist.root \
        $refdir/$path/tracksummary_ckf_hist.root \
        "CKF track summary | ${name}" \
        $path/tracksummary_ckf.html \
        $path/tracksummary_ckf_plots

    if [ -f $refdir/$path/performance_finding_ckf_ambi.root ]; then
        run_histcmp \
            $outdir/data/$path/performance_finding_ckf_ambi.root \
            $refdir/$path/performance_finding_ckf_ambi.root \
            "Ambisolver finding performance | ${name}" \
            $path/performance_finding_ckf_ambi.html \
            $path/performance_finding_ckf_ambi
    fi

    if [ -f $refdir/$path/performance_finding_ckf_ml_solver.root ]; then
        run_histcmp \
            $outdir/data/$path/performance_finding_ckf_ml_solver.root \
            $refdir/$path/performance_finding_ckf_ml_solver.root \
            "ML Ambisolver | ${name}" \
            $path/performance_finding_ckf_ml_solver.html \
            $path/performance_finding_ckf_ml_solver
    fi
}

function vertexing() {
    name=$1
    path=$2
    config=$3

    if [ -f $refdir/$path/performance_vertexing_ivf_notime_hist.root ]; then
        run Examples/Scripts/generic_plotter.py \
            $outdir/data/$path/performance_vertexing_ivf_notime.root \
            vertexing \
            $outdir/data/$path/performance_vertexing_ivf_notime_hist.root \
            --silent \
            --config $config
        ec=$(($ec | $?))

        # remove ntuple file because it's large
        rm $outdir/data/$path/performance_vertexing_ivf_notime.root

        run_histcmp \
            $outdir/data/$path/performance_vertexing_ivf_notime_hist.root \
            $refdir/$path/performance_vertexing_ivf_notime_hist.root \
            "IVF notime | ${name}" \
            $path/performance_vertexing_ivf_notime.html \
            $path/performance_vertexing_ivf_notime_plots
    fi

    run Examples/Scripts/generic_plotter.py \
        $outdir/data/$path/performance_vertexing_amvf_gauss_notime.root \
        vertexing \
        $outdir/data/$path/performance_vertexing_amvf_gauss_notime_hist.root \
        --silent \
        --config $config
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/data/$path/performance_vertexing_amvf_gauss_notime.root

    run_histcmp \
        $outdir/data/$path/performance_vertexing_amvf_gauss_notime_hist.root \
        $refdir/$path/performance_vertexing_amvf_gauss_notime_hist.root \
        "AMVF gauss notime | ${name}" \
        $path/performance_vertexing_amvf_gauss_notime.html \
        $path/performance_vertexing_amvf_gauss_notime_plots

    run Examples/Scripts/generic_plotter.py \
        $outdir/data/$path/performance_vertexing_amvf_grid_time.root \
        vertexing \
        $outdir/data/$path/performance_vertexing_amvf_grid_time_hist.root \
        --silent \
        --config $config
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/data/$path/performance_vertexing_amvf_grid_time.root

    run_histcmp \
        $outdir/data/$path/performance_vertexing_amvf_grid_time_hist.root \
        $refdir/$path/performance_vertexing_amvf_grid_time_hist.root \
        "AMVF grid time | ${name}" \
        $path/performance_vertexing_amvf_grid_time.html \
        $path/performance_vertexing_amvf_grid_time_plots
}

function simulation() {
    suffix=$1

    config="CI/physmon/config/simulation.yml"

    run Examples/Scripts/generic_plotter.py \
        $outdir/data/simulation/particles_${suffix}.root \
        particles \
        $outdir/data/simulation/particles_${suffix}_hist.root \
        --silent \
        --config $config
    ec=$(($ec | $?))

    # remove ntuple file because it's large
    rm $outdir/data/simulation/particles_${suffix}.root

    run_histcmp \
        $outdir/data/simulation/particles_${suffix}_hist.root \
        $refdir/simulation/particles_${suffix}_hist.root \
        "Particles ${suffix}" \
        simulation/particles_${suffix}.html \
        simulation/particles_${suffix}_plots
}

function generation() {
    run Examples/Scripts/generic_plotter.py \
        $outdir/data/simulation/particles_ttbar.root \
        particles \
        $outdir/data/simulation/particles_ttbar_hist.root \
        --silent \
        --config CI/physmon/config/pythia8_ttbar.yml

    # remove ntuple file because it's large
    rm $outdir/data/simulation/particles_ttbar.root

    run_histcmp \
        $outdir/data/simulation/particles_ttbar_hist.root \
        $refdir/simulation/particles_ttbar_hist.root \
        "Particles ttbar" \
        simulation/particles_ttbar.html \
        simulation/particles_ttbar_plots

    run Examples/Scripts/generic_plotter.py \
        $outdir/data/simulation/vertices_ttbar.root \
        vertices \
        $outdir/data/simulation/vertices_ttbar_hist.root \
        --silent \
        --config CI/physmon/config/pythia8_ttbar.yml

    # remove ntuple file because it's large
    rm $outdir/data/simulation/vertices_ttbar.root

    run_histcmp \
        $outdir/data/simulation/vertices_ttbar_hist.root \
        $refdir/simulation/vertices_ttbar_hist.root \
        "Vertices ttbar" \
        simulation/vertices_ttbar.html \
        simulation/vertices_ttbar_plots
}

if [[ "$mode" == "all" || "$mode" == "simulation" ]]; then
    simulation fatras
    simulation geant4

    generation
fi

if [[ "$mode" == "all" || "$mode" == "kf" ]]; then
    run_histcmp \
        $outdir/data/trackfitting_kf/performance_trackfitting.root \
        $refdir/trackfitting_kf/performance_trackfitting.root \
        "Truth tracking (KF)" \
        trackfitting_kf/performance_trackfitting.html \
        trackfitting_kf/performance_trackfitting_plots \
        --config CI/physmon/config/trackfitting_kf.yml
fi

if [[ "$mode" == "all" || "$mode" == "gsf" ]]; then
    run_histcmp \
        $outdir/data/trackfitting_gsf/performance_trackfitting.root \
        $refdir/trackfitting_gsf/performance_trackfitting.root \
        "Truth tracking (GSF)" \
        trackfitting_gsf/performance_trackfitting.html \
        trackfitting_gsf/performance_trackfitting_plots \
        --config CI/physmon/config/trackfitting_gsf.yml
fi

if [[ "$mode" == "all" || "$mode" == "gx2f" ]]; then
    run_histcmp \
        $outdir/data/trackfitting_gx2f/performance_trackfitting.root \
        $refdir/trackfitting_gx2f/performance_trackfitting.root \
        "Truth tracking (GX2F)" \
        trackfitting_gx2f/performance_trackfitting.html \
        trackfitting_gx2f/performance_trackfitting_plots \
        --config CI/physmon/config/trackfitting_gx2f.yml
fi

if [[ "$mode" == "all" || "$mode" == "kf_refit" ]]; then
    run_histcmp \
        $outdir/data/trackrefitting_kf/performance_trackrefitting.root \
        $refdir/trackrefitting_kf/performance_trackrefitting.root \
        "Truth tracking (KF refit)" \
        trackrefitting_kf/performance_trackrefitting.html \
        trackrefitting_kf/performance_trackrefitting_plots \
        --config CI/physmon/config/trackfitting_kf.yml
fi

if [[ "$mode" == "all" || "$mode" == "gsf_refit" ]]; then
    run_histcmp \
        $outdir/data/trackrefitting_gsf/performance_trackrefitting.root \
        $refdir/trackrefitting_gsf/performance_trackrefitting.root \
        "Truth tracking (GSF refit)" \
        trackrefitting_gsf/performance_trackrefitting.html \
        trackrefitting_gsf/performance_trackrefitting_plots \
        --config CI/physmon/config/trackfitting_gsf.yml
fi

if [[ "$mode" == "all" || "$mode" == "fullchains" ]]; then
    trackfinding "trackfinding | single muon | truth smeared seeding" trackfinding_1muon/truth_smeared
    trackfinding "trackfinding | single muon | truth estimated seeding" trackfinding_1muon/truth_estimated
    trackfinding "trackfinding | single muon | default seeding" trackfinding_1muon/seeded
    trackfinding "trackfinding | single muon | orthogonal seeding" trackfinding_1muon/orthogonal

    trackfinding "trackfinding | 4 muon x 50 vertices | default seeding" trackfinding_4muon_50vertices
    vertexing "trackfinding | 4 muon x 50 vertices | default seeding" trackfinding_4muon_50vertices CI/physmon/config/vertexing_4muon_50vertices.yml

    trackfinding "trackfinding | ttbar with 200 pileup | default seeding" trackfinding_ttbar_pu200
    vertexing "trackfinding | ttbar with 200 pileup | default seeding" trackfinding_ttbar_pu200 CI/physmon/config/vertexing_ttbar_pu200.yml
fi

if [[ "$mode" == "all" || "$mode" == "gx2f_vs_kf" ]]; then
    run_histcmp \
        $outdir/data/trackfitting_gx2f_vs_kf/performance_trackfitting_gx2f.root \
        $outdir/data/trackfitting_gx2f_vs_kf/performance_trackfitting_kf.root \
        "Comparison - Truth tracking (GX2F vs KF)" \
        trackfitting_gx2f_vs_kf/performance_trackfitting.html \
        trackfitting_gx2f_vs_kf/performance_trackfitting_plots \
        --config CI/physmon/config/info_only.yml \
        --label-reference=KF \
        --label-monitored=GX2F
fi

run CI/physmon/summary.py $histcmp_results \
  --md $outdir/summary.md \
  --html $outdir/summary.html
ec=$(($ec | $?))

exit $ec
