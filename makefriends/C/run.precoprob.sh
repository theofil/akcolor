#!/bin/bash
set -euo pipefail

# Produces friends for the precoprob (modified probability-of-color-reconnection,
# 0.95 -> 0.25) Herwig samples into friends/summer26/precoprob. Same xs, same boson
# tags, same makefriends binary/preselection as the nominal summer26 production
# (run.summer26.sh) -- only the input files and output directory differ. Herwig only:
# there is no MG5+Pythia precoprob variant.

WORKDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$WORKDIR"

set +u
source "$WORKDIR/init.source"
set -u
make

# ── Output directory ──────────────────────────────────────────────────────────
OUTDIR="friends/summer26/precoprob"
mkdir -p "$OUTDIR"

# ── Temp directory for all transient files (local disk to avoid XRootD issues) ──
JOBDIR="/tmp/${USER}_akcolor_precoprob_$$"
mkdir -p "$JOBDIR"

# ── Parse arguments ───────────────────────────────────────────────────────────
GO_FAST=false
GO_FAST_N=100
NJOBS=1           # samples processed in parallel; each sample runs all its jobs at once; override with --jobs N
for (( i=1; i<=$#; i++ )); do
    arg="${!i}"
    if [[ "$arg" == "--goFast" ]]; then
        GO_FAST=true
        next=$((i+1))
        if [[ $next -le $# && "${!next}" =~ ^[0-9]+$ ]]; then
            GO_FAST_N="${!next}"
            i=$next
        fi
    elif [[ "$arg" == "--jobs" ]]; then
        next=$((i+1))
        if [[ $next -le $# && "${!next}" =~ ^[0-9]+$ ]]; then
            NJOBS="${!next}"
            i=$next
        fi
    fi
done

# ── Campaign source directories ───────────────────────────────────────────────
CAMPAIGN41="/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00041/merged"
CAMPAIGN42="/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00042/merged"

# ── Symlinks (inside tmp) ─────────────────────────────────────────────────────
declare -A ALL_SYMLINKS
# campaign-00041 — precoprob Herwig
ALL_SYMLINKS["VBFH_herwig_preco.root"]="$CAMPAIGN41/VBFH_herwig_preco.root"
ALL_SYMLINKS["VBFW_herwig_preco.root"]="$CAMPAIGN41/VBFW_herwig_preco.root"
ALL_SYMLINKS["VBFZ_herwig_preco.root"]="$CAMPAIGN41/VBFZ_herwig_preco.root"
ALL_SYMLINKS["QCDWjj_herwig_preco.root"]="$CAMPAIGN41/QCDWjj_herwig_preco.root"
ALL_SYMLINKS["QCDZjj_herwig_preco.root"]="$CAMPAIGN41/QCDZjj_herwig_preco.root"
# campaign-00042 — precoprob Herwig QCDHjj
ALL_SYMLINKS["QCDHjj_herwig_preco.root"]="$CAMPAIGN42/QCDHjj_herwig_preco.root"

for link in "${!ALL_SYMLINKS[@]}"; do
    if [ ! -L "$JOBDIR/$link" ]; then
        ln -s "${ALL_SYMLINKS[$link]}" "$JOBDIR/$link"
        echo "Created symlink: $JOBDIR/$link -> ${ALL_SYMLINKS[$link]}"
    fi
done

# Cross sections (pb) -- identical to nominal summer26 (same process/boson):
# VBFH  herwig: 2.97189
# VBFW  herwig: 7.494
# VBFZ  herwig: 1.099
# QCDHjj herwig: 5.06051
# QCDWjj herwig: 1733.3
# QCDZjj herwig: 340.2

if $GO_FAST; then
    # ── goFast: single foreground run, GO_FAST_N events per sample ───────────
    FAST_OUTDIR="$WORKDIR/fast_tmp_precoprob"
    mkdir -p "$FAST_OUTDIR"
    echo "[$(date +%H:%M:%S)] --goFast: 1 x ${GO_FAST_N} events per sample → $FAST_OUTDIR"
    ./makefriends "$JOBDIR/VBFH_herwig_preco.root"   --genWeight evweight --xs 2.97189  --boson H --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/VBFH_herwig.friend.root"
    ./makefriends "$JOBDIR/VBFW_herwig_preco.root"   --genWeight evweight --xs 7.494    --boson W --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/VBFW_herwig.friend.root"
    ./makefriends "$JOBDIR/VBFZ_herwig_preco.root"   --genWeight evweight --xs 1.099    --boson Z --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/VBFZ_herwig.friend.root"
    ./makefriends "$JOBDIR/QCDHjj_herwig_preco.root" --genWeight evweight --xs 5.06051  --boson H --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/QCDHjj_herwig.friend.root"
    ./makefriends "$JOBDIR/QCDWjj_herwig_preco.root" --genWeight evweight --xs 1733.3   --boson W --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/QCDWjj_herwig.friend.root"
    ./makefriends "$JOBDIR/QCDZjj_herwig_preco.root" --genWeight evweight --xs 340.2    --boson Z --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/QCDZjj_herwig.friend.root"
    echo "[$(date +%H:%M:%S)] Done. Output in $FAST_OUTDIR"
else
    # ── Full production: NJOBS samples in parallel, each sample fires all jobs at once ──
    run_sample() {
        local file=$1 njobs=$2 tot=$3 xs=$4 outfile=$5 boson=$6
        local stem="${file%.root}"
        echo "[$(date +%H:%M:%S)] START $file ($njobs jobs, xs=${xs} pb, boson=${boson})"
        local job_pids=()
        for (( B=0; B<njobs; B++ )); do
            (cd "$JOBDIR" && "$WORKDIR/makefriends" "$file" --genWeight evweight --xs "$xs" \
                --boson "$boson" --skip $((B * tot)) --totEve "$tot" --jobId "$B") &
            job_pids+=($!)
        done
        for pid in "${job_pids[@]}"; do wait "$pid"; done
        hadd -f "$OUTDIR/$outfile" "$JOBDIR"/${stem}*_*.root
        echo "[$(date +%H:%M:%S)] DONE  $outfile"
    }

    # file:njobs:tot:xs:outfile:boson  (outfile drops the _preco suffix so spva_plots.py
    # can pair it 1:1 by sample name against friends/summer26/<name>.friend.root)
    samples=(
        "VBFH_herwig_preco.root:20:50000:2.97189:VBFH_herwig.friend.root:H"
        "VBFW_herwig_preco.root:20:50000:7.494:VBFW_herwig.friend.root:W"
        "VBFZ_herwig_preco.root:20:50000:1.099:VBFZ_herwig.friend.root:Z"
        "QCDHjj_herwig_preco.root:20:50000:5.06051:QCDHjj_herwig.friend.root:H"
        "QCDWjj_herwig_preco.root:20:50000:1733.3:QCDWjj_herwig.friend.root:W"
        "QCDZjj_herwig_preco.root:20:50000:340.2:QCDZjj_herwig.friend.root:Z"
    )

    fail=0
    i=0
    while (( i < ${#samples[@]} )); do
        batch_pids=()
        for (( j=i; j < i+NJOBS && j < ${#samples[@]}; j++ )); do
            IFS=: read -r file njobs tot xs outfile boson <<< "${samples[$j]}"
            run_sample "$file" "$njobs" "$tot" "$xs" "$outfile" "$boson" & batch_pids+=($!)
        done
        for pid in "${batch_pids[@]}"; do wait "$pid" || fail=1; done
        (( i += NJOBS ))
    done
    (( fail )) && { echo "[ERROR] One or more samples failed; check output above."; exit 1; }

    # ── Remove tmp (symlinks + per-job intermediates) ─────────────────────────
    echo "Removing tmp directory..."
    rm -rf "$JOBDIR"

    echo "[$(date +%H:%M:%S)] Done."
fi
