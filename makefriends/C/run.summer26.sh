#!/bin/bash
set -euo pipefail

WORKDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$WORKDIR"

set +u
source "$WORKDIR/init.source"
set -u
make

# ── Output directory ──────────────────────────────────────────────────────────
OUTDIR="friends/summer26"
mkdir -p "$OUTDIR"

# ── Temp directory for all transient files (local disk to avoid XRootD issues) ──
JOBDIR="/tmp/${USER}_akcolor_$$"
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
CAMPAIGN23="/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00023/merged"
CAMPAIGN30="/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00030/merged"
CAMPAIGN31="/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/campaigns/campaign-00031/merged"

# ── Symlinks (inside tmp) ─────────────────────────────────────────────────────
declare -A ALL_SYMLINKS
# campaign-00023 — Herwig + MG5+Pythia pairs
ALL_SYMLINKS["VBFH_herwig.root"]="$CAMPAIGN23/VBFH_herwig.root"
ALL_SYMLINKS["VBFH_mg5_pythia.root"]="$CAMPAIGN23/VBFH_mg5_pythia.root"
ALL_SYMLINKS["VBFW_herwig.root"]="$CAMPAIGN23/VBFW_herwig.root"
ALL_SYMLINKS["VBFW_mg5_pythia.root"]="$CAMPAIGN23/VBFW_mg5_pythia.root"
ALL_SYMLINKS["VBFZ_herwig.root"]="$CAMPAIGN23/VBFZ_herwig.root"
ALL_SYMLINKS["VBFZ_mg5_pythia.root"]="$CAMPAIGN23/VBFZ_mg5_pythia.root"
ALL_SYMLINKS["QCDWjj_herwig.root"]="$CAMPAIGN23/QCDWjj_herwig.root"
ALL_SYMLINKS["QCDWjj_mg5_pythia.root"]="$CAMPAIGN23/QCDWjj_mg5_pythia.root"
ALL_SYMLINKS["QCDZjj_herwig.root"]="$CAMPAIGN23/QCDZjj_herwig.root"
ALL_SYMLINKS["QCDZjj_mg5_pythia.root"]="$CAMPAIGN23/QCDZjj_mg5_pythia.root"
# campaign-00030 — QCDHjj Herwig only
ALL_SYMLINKS["QCDHjj_herwig.root"]="$CAMPAIGN30/QCDHjj_herwig.root"
# campaign-00031 — QCDHjj MG5+Pythia
ALL_SYMLINKS["QCDHjj_mg5_pythia.root"]="$CAMPAIGN31/QCDHjj_mg5_pythia.root"

for link in "${!ALL_SYMLINKS[@]}"; do
    if [ ! -L "$JOBDIR/$link" ]; then
        ln -s "${ALL_SYMLINKS[$link]}" "$JOBDIR/$link"
        echo "Created symlink: $JOBDIR/$link -> ${ALL_SYMLINKS[$link]}"
    fi
done

# Cross sections (pb):
# VBFH  herwig: 2.97189   mg5_pythia: 2.88990
# VBFW  herwig: 7.494      mg5_pythia: 7.20333
# VBFZ  herwig: 1.099      mg5_pythia: 1.08581
# QCDHjj herwig: 5.06051   mg5_pythia: 4.982 (camp-00031)
# QCDWjj herwig: 1733.3    mg5_pythia: 1646.805
# QCDZjj herwig: 340.2     mg5_pythia: 316.363

if $GO_FAST; then
    # ── goFast: single foreground run, GO_FAST_N events per sample ───────────
    # Output goes to fast_tmp/ to avoid overwriting production files.
    FAST_OUTDIR="$WORKDIR/fast_tmp"
    mkdir -p "$FAST_OUTDIR"
    echo "[$(date +%H:%M:%S)] --goFast: 1 x ${GO_FAST_N} events per sample → $FAST_OUTDIR"
    ./makefriends "$JOBDIR/VBFH_herwig.root"        --genWeight evweight --xs 2.97189  --boson H --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/VBFH_herwig.friend.root"
    ./makefriends "$JOBDIR/VBFH_mg5_pythia.root"    --genWeight evweight --xs 2.88990  --boson H --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/VBFH_mg5_pythia.friend.root"
    ./makefriends "$JOBDIR/VBFW_herwig.root"        --genWeight evweight --xs 7.494    --boson W --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/VBFW_herwig.friend.root"
    ./makefriends "$JOBDIR/VBFW_mg5_pythia.root"    --genWeight evweight --xs 7.20333  --boson W --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/VBFW_mg5_pythia.friend.root"
    ./makefriends "$JOBDIR/VBFZ_herwig.root"        --genWeight evweight --xs 1.099    --boson Z --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/VBFZ_herwig.friend.root"
    ./makefriends "$JOBDIR/VBFZ_mg5_pythia.root"    --genWeight evweight --xs 1.08581  --boson Z --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/VBFZ_mg5_pythia.friend.root"
    ./makefriends "$JOBDIR/QCDHjj_herwig.root"      --genWeight evweight --xs 5.06051  --boson H --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/QCDHjj_herwig.friend.root"
    ./makefriends "$JOBDIR/QCDHjj_mg5_pythia.root"  --genWeight evweight --xs 4.982    --boson H --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/QCDHjj_mg5_pythia.friend.root"
    ./makefriends "$JOBDIR/QCDWjj_herwig.root"      --genWeight evweight --xs 1733.3   --boson W --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/QCDWjj_herwig.friend.root"
    ./makefriends "$JOBDIR/QCDWjj_mg5_pythia.root"  --genWeight evweight --xs 1646.805 --boson W --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/QCDWjj_mg5_pythia.friend.root"
    ./makefriends "$JOBDIR/QCDZjj_herwig.root"      --genWeight evweight --xs 340.2    --boson Z --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/QCDZjj_herwig.friend.root"
    ./makefriends "$JOBDIR/QCDZjj_mg5_pythia.root"  --genWeight evweight --xs 316.363  --boson Z --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/QCDZjj_mg5_pythia.friend.root"
    echo "[$(date +%H:%M:%S)] goFast runs complete. Converting to HDF5..."
    python3 root_to_hdf5.py "$FAST_OUTDIR"
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

    # file:njobs:tot:xs:outfile:boson
    samples=(
        "VBFH_herwig.root:20:50000:2.97189:VBFH_herwig.friend.root:H"
        "VBFH_mg5_pythia.root:20:50000:2.88990:VBFH_mg5_pythia.friend.root:H"
        "VBFW_herwig.root:20:50000:7.494:VBFW_herwig.friend.root:W"
        "VBFW_mg5_pythia.root:20:50000:7.20333:VBFW_mg5_pythia.friend.root:W"
        "VBFZ_herwig.root:20:50000:1.099:VBFZ_herwig.friend.root:Z"
        "VBFZ_mg5_pythia.root:20:50000:1.08581:VBFZ_mg5_pythia.friend.root:Z"
        "QCDHjj_herwig.root:20:50000:5.06051:QCDHjj_herwig.friend.root:H"
        "QCDHjj_mg5_pythia.root:20:50000:4.982:QCDHjj_mg5_pythia.friend.root:H"
        "QCDWjj_herwig.root:20:50000:1733.3:QCDWjj_herwig.friend.root:W"
        "QCDWjj_mg5_pythia.root:20:50000:1646.805:QCDWjj_mg5_pythia.friend.root:W"
        "QCDZjj_herwig.root:20:50000:340.2:QCDZjj_herwig.friend.root:Z"
        "QCDZjj_mg5_pythia.root:20:50000:316.363:QCDZjj_mg5_pythia.friend.root:Z"
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

    # ── Convert final friend.root files to HDF5 ───────────────────────────────
    echo "[$(date +%H:%M:%S)] Converting *.friend.root to HDF5 in $OUTDIR..."
    python3 root_to_hdf5.py "$OUTDIR"

    # ── Remove tmp (symlinks + per-job intermediates) ─────────────────────────
    echo "Removing tmp directory..."
    rm -rf "$JOBDIR"

    echo "[$(date +%H:%M:%S)] Done."
fi
