#!/bin/bash
set -euo pipefail

# nov26 production (see nov26.md): inclusive friends -- every input event is written, the
# dijet gate is stored as a flag -- with <= 2 jets + constituents (CMS-PF-like particle type
# and charge), built by makefriends_nov26. Z channel only for now: VBFZ + QCDZjj, Herwig
# (campaign-00023) and MG5+Pythia (pythia-dipole-recoil), same xs as summer26.
# summer26 (makefriends.cpp, run.summer26.sh, friends/summer26) is not touched.
#
#   ./run.nov26.sh --goFast 2000   # quick test: first N events per sample -> fast_tmp_nov26/
#   ./run.nov26.sh                 # full production: 20 jobs x 50k events per sample -> friends/nov26/

WORKDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$WORKDIR"

set +u
source "$WORKDIR/init.source"
set -u
make makefriends_nov26

# ── Output directory ──────────────────────────────────────────────────────────
OUTDIR="friends/nov26"
mkdir -p "$OUTDIR"

# ── Temp directory for all transient files (local disk to avoid XRootD issues) ──
JOBDIR="/tmp/${USER}_akcolor_nov26_$$"
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
PYDIPOLE="/eos/user/a/apapaefs/Projects/Pull/pythia-dipole-recoil/merged"

# ── Symlinks (inside tmp) ─────────────────────────────────────────────────────
declare -A ALL_SYMLINKS
# campaign-00023 — Herwig
ALL_SYMLINKS["VBFZ_herwig.root"]="$CAMPAIGN23/VBFZ_herwig.root"
ALL_SYMLINKS["QCDZjj_herwig.root"]="$CAMPAIGN23/QCDZjj_herwig.root"
# pythia-dipole-recoil — MG5+Pythia
for proc in VBFZ QCDZjj; do
    ALL_SYMLINKS["${proc}_mg5_pythia.root"]="$PYDIPOLE/$proc/${proc}_mg5_pythia_dipolerecoil.root"
done

for link in "${!ALL_SYMLINKS[@]}"; do
    if [ ! -L "$JOBDIR/$link" ]; then
        ln -s "${ALL_SYMLINKS[$link]}" "$JOBDIR/$link"
        echo "Created symlink: $JOBDIR/$link -> ${ALL_SYMLINKS[$link]}"
    fi
done

# Cross sections (pb), identical to summer26:
# VBFZ   herwig: 1.099    mg5_pythia: 1.08581
# QCDZjj herwig: 340.2    mg5_pythia: 316.363

# file:njobs:tot:xs:outfile:boson
samples=(
    "VBFZ_herwig.root:20:50000:1.099:VBFZ_herwig.friend.root:Z"
    "VBFZ_mg5_pythia.root:20:50000:1.08581:VBFZ_mg5_pythia.friend.root:Z"
    "QCDZjj_herwig.root:20:50000:340.2:QCDZjj_herwig.friend.root:Z"
    "QCDZjj_mg5_pythia.root:20:50000:316.363:QCDZjj_mg5_pythia.friend.root:Z"
)

# Outputs are always deleted before they are rewritten: overwriting a file in place on EOS
# leaves hidden .sys.v#.<name> version dirs behind.

if $GO_FAST; then
    # ── goFast: single foreground run, first GO_FAST_N events per sample ──────
    # Output goes to fast_tmp_nov26/ to avoid overwriting production files.
    FAST_OUTDIR="$WORKDIR/fast_tmp_nov26"
    mkdir -p "$FAST_OUTDIR"
    echo "[$(date +%H:%M:%S)] --goFast: 1 x ${GO_FAST_N} events per sample → $FAST_OUTDIR"
    fast_files=()
    for s in "${samples[@]}"; do
        IFS=: read -r file njobs tot xs outfile boson <<< "$s"
        rm -f "$FAST_OUTDIR/$outfile" "$FAST_OUTDIR/${outfile%.friend.root}.h5"
        ./makefriends_nov26 "$JOBDIR/$file" --genWeight evweight --xs "$xs" --boson "$boson" \
            --totEve "$GO_FAST_N" --output "$FAST_OUTDIR/$outfile"
        fast_files+=("$FAST_OUTDIR/$outfile")
    done
    echo "[$(date +%H:%M:%S)] goFast runs complete. Converting to HDF5..."
    python3 root_to_hdf5_nov26.py "${fast_files[@]}"
    rm -rf "$JOBDIR"
    echo "[$(date +%H:%M:%S)] Done. Output in $FAST_OUTDIR"
else
    # ── Full production: NJOBS samples in parallel, each sample fires all jobs at once ──
    run_sample() {
        local file=$1 njobs=$2 tot=$3 xs=$4 outfile=$5 boson=$6
        local stem="${file%.root}"
        echo "[$(date +%H:%M:%S)] START $file ($njobs jobs, xs=${xs} pb, boson=${boson})"
        local job_pids=() parts=()
        for (( B=0; B<njobs; B++ )); do
            (cd "$JOBDIR" && "$WORKDIR/makefriends_nov26" "$file" --genWeight evweight --xs "$xs" \
                --boson "$boson" --skip $((B * tot)) --totEve "$tot" --jobId "$B") &
            job_pids+=($!)
            parts+=("$JOBDIR/${stem}slimfriend_${B}.root")
        done
        local rc=0
        for pid in "${job_pids[@]}"; do wait "$pid" || rc=1; done
        if (( rc )); then echo "[ERROR] a makefriends_nov26 job failed for $file"; return 1; fi
        # parts are listed in numeric job order (a glob would put _10 before _2), so friend
        # entry i is input entry i (checked downstream via the evIdx branch)
        rm -f "$OUTDIR/$outfile"
        hadd -f "$OUTDIR/$outfile" "${parts[@]}"
        echo "[$(date +%H:%M:%S)] DONE  $outfile"
    }

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

    # ── Convert this run's friend.root files (explicit list) to HDF5 ──────────
    out_files=()
    for s in "${samples[@]}"; do
        IFS=: read -r file njobs tot xs outfile boson <<< "$s"
        out_files+=("$OUTDIR/$outfile")
    done
    echo "[$(date +%H:%M:%S)] Converting ${#out_files[@]} friend.root files to HDF5..."
    python3 root_to_hdf5_nov26.py "${out_files[@]}"

    # ── Remove tmp (symlinks + per-job intermediates) ─────────────────────────
    echo "Removing tmp directory..."
    rm -rf "$JOBDIR"

    echo "[$(date +%H:%M:%S)] Done."
fi
