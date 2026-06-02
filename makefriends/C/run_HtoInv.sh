#!/bin/bash
set -euo pipefail

WORKDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$WORKDIR"

make

# ── Parse arguments ───────────────────────────────────────────────────────────
GO_FAST=false
GO_FAST_N=100
for (( i=1; i<=$#; i++ )); do
    arg="${!i}"
    if [[ "$arg" == "--goFast" ]]; then
        GO_FAST=true
        next=$((i+1))
        if [[ $next -le $# && "${!next}" =~ ^[0-9]+$ ]]; then
            GO_FAST_N="${!next}"
            i=$next
        fi
    fi
done

# ── Symlinks ──────────────────────────────────────────────────────────────────
declare -A ALL_SYMLINKS
ALL_SYMLINKS["VBFHtoInv.root"]="/eos/home-a/apapaefs/Projects/Pull/Summer2025/Herwig_13.6/LHC-Matchbox-VBFH-13.6.root"
ALL_SYMLINKS["VBFWtoInv.root"]="/eos/home-a/apapaefs/Projects/Pull/Summer2025/Herwig_13.6/LHC-Matchbox-VBFW-13.6.root"
ALL_SYMLINKS["VBFZtoInv.root"]="/eos/home-a/apapaefs/Projects/Pull/Summer2025/Herwig_13.6/LHC-Matchbox-VBFZ-13.6.root"
ALL_SYMLINKS["QCDHtoInv.root"]="/eos/home-a/apapaefs/Projects/Pull/Summer2025/Herwig_13.6/LHC-LHEMG5-ggHjj-13.6.root"
ALL_SYMLINKS["QCDWtoInv.root"]="/eos/home-a/apapaefs/Projects/Pull/Summer2025/Herwig_13.6/LHC-LHEMG5-ppWjj-13.6.root"
ALL_SYMLINKS["QCDZtoInv.root"]="/eos/home-a/apapaefs/Projects/Pull/Summer2025/Herwig_13.6/LHC-LHEMG5-ppZjj-13.6.root"

CREATED_SYMLINKS=()
for link in "${!ALL_SYMLINKS[@]}"; do
    if [ ! -L "$link" ]; then
        ln -s "${ALL_SYMLINKS[$link]}" "$link"
        CREATED_SYMLINKS+=("$link")
        echo "Created symlink: $link -> ${ALL_SYMLINKS[$link]}"
    fi
done

if $GO_FAST; then
    # ── goFast: single foreground run, GO_FAST_N events per sample ───────────
    echo "[$(date +%H:%M:%S)] --goFast: 1 x ${GO_FAST_N} events per sample, running sequentially..."
    ./makefriends VBFHtoInv.root --genWeight evweight --xs 3.901   --totEve "$GO_FAST_N" --output VBFHtoInv.friend.root
    ./makefriends QCDHtoInv.root --genWeight evweight --xs 2.26114 --totEve "$GO_FAST_N" --output QCDHtoInv.friend.root
    echo "[$(date +%H:%M:%S)] goFast runs complete."
else
    # ── Launch jobs (inline so bash `wait` can reap them) ────────────────────
    launch_jobs() {
        local file=$1 njobs=$2 tot=$3 xs=$4
        for (( JOB=0; JOB<njobs; JOB++ )); do
            ./makefriends "$file" --genWeight evweight --xs "$xs" \
                --skip $((JOB * tot)) --totEve "$tot" --jobId "$JOB" &
        done
    }

    echo "[$(date +%H:%M:%S)] Launching VBFHtoInv jobs (20 x 50000 events, xs=3.901 pb)..."
    launch_jobs VBFHtoInv.root 20 50000 3.901

    echo "[$(date +%H:%M:%S)] Launching QCDHtoInv jobs (20 x 50000 events, xs=2.26114 pb)..."
    launch_jobs QCDHtoInv.root 20 50000 2.26114

    # ── Wait for all background jobs to finish ────────────────────────────────
    echo "Waiting for all jobs to complete..."
    wait
    echo "[$(date +%H:%M:%S)] All jobs finished."

    # ── Merge with hadd (always overwrite — we just produced the inputs) ──────
    echo "Merging VBFHtoInv output files..."
    hadd -f VBFHtoInv.friend.root VBFHtoInv*_*.root

    echo "Merging QCDHtoInv output files..."
    hadd -f QCDHtoInv.friend.root QCDHtoInv*_*.root
fi

# ── Collect output files into friends/ ───────────────────────────────────────
mkdir -p friends
mv VBFHtoInv.friend.root friends/
mv QCDHtoInv.friend.root friends/
echo "Moved slimfriend files to friends/"

# ── Remove per-job intermediate files ────────────────────────────────────────
echo "Removing intermediate per-job files..."
rm -f *friend_*.root

# ── Remove symlinks ───────────────────────────────────────────────────────────
for link in "${!ALL_SYMLINKS[@]}"; do
    if [ -L "$link" ]; then
        rm -f "$link"
        echo "Removed symlink: $link"
    fi
done

echo "[$(date +%H:%M:%S)] Done."
