#!/bin/bash
set -euo pipefail

WORKDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$WORKDIR"

# ── Parse arguments ───────────────────────────────────────────────────────────
GO_FAST=false
for arg in "$@"; do
    [[ "$arg" == "--goFast" ]] && GO_FAST=true
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
    # ── goFast: single foreground run, 1000 events per sample ─────────────────
    echo "[$(date +%H:%M:%S)] --goFast: 1 x 1000 events per sample, running sequentially..."
    ./makefriends VBFHtoInv.root --genWeight evweight --xs 3.901   --totEve 1000 --output VBFHtoInv.slimfriend.root
    ./makefriends QCDHtoInv.root --genWeight evweight --xs 2.26114 --totEve 1000 --output QCDHtoInv.slimfriend.root
    echo "[$(date +%H:%M:%S)] goFast runs complete."
else
    # ── Create jobs.sh dynamically ────────────────────────────────────────────
    cat > jobs.sh << 'JOBSEOF'
#!/bin/bash

usage() {
    echo "Usage: $0 -f FILE -n nJOBS -t TOT -x XS"
    echo ""
    echo "  -f FILE     Input ROOT file"
    echo "  -n nJOBS    Number of parallel jobs"
    echo "  -t TOT      Events per job"
    echo "  -x XS       Cross-section in pb"
    exit 1
}

FILE=""
nJOBS=0
TOT=0
XS=0.0

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -f|--file)  FILE="$2";  shift ;;
        -n|--njobs) nJOBS="$2"; shift ;;
        -t|--tot)   TOT="$2";   shift ;;
        -x|--xs)    XS="$2";    shift ;;
        --help)     usage ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
    shift
done

if [[ -z "$FILE" || "$nJOBS" -eq 0 || "$TOT" -eq 0 ]]; then
    echo "Error: missing required arguments."
    usage
fi

for (( JOB=0; JOB<nJOBS; JOB++ )); do
    SKIP=$((JOB * TOT))
    ./makefriends "$FILE" --genWeight evweight --xs "$XS" --skip "$SKIP" --totEve "$TOT" --jobId "$JOB" &
done
JOBSEOF
    chmod +x jobs.sh
    echo "Created jobs.sh"

    # ── Launch jobs ───────────────────────────────────────────────────────────
    echo "[$(date +%H:%M:%S)] Launching VBFHtoInv jobs (20 x 50000 events, xs=3.901 pb)..."
    ./jobs.sh -f VBFHtoInv.root -n 20 -t 50000 -x 3.901

    echo "[$(date +%H:%M:%S)] Launching QCDHtoInv jobs (20 x 50000 events, xs=2.26114 pb)..."
    ./jobs.sh -f QCDHtoInv.root -n 20 -t 50000 -x 2.26114

    # ── Wait for all output files to exist AND be closed ──────────────────────
    # ROOT creates the file immediately on open but only finalizes it on close.
    # fuser returns 0 if any process holds the file open, non-zero when closed.
    wait_for_files() {
        local stem=$1   # e.g. VBFHtoInv
        local n=$2
        while true; do
            local files=()
            while IFS= read -r -d '' f; do
                files+=("$f")
            done < <(find . -maxdepth 1 -name "${stem}slimfriend_*.root" -print0 2>/dev/null)

            local count=${#files[@]}
            local closed=0
            for f in "${files[@]}"; do
                fuser "$f" &>/dev/null || closed=$((closed + 1))
            done

            echo "[$(date +%H:%M:%S)] $stem: $count/$n files exist, $closed closed"
            [ "$count" -ge "$n" ] && [ "$closed" -ge "$n" ] && break
            sleep 30
        done
    }

    echo "Waiting for jobs to complete..."
    wait_for_files VBFHtoInv 20
    wait_for_files QCDHtoInv 20
    echo "[$(date +%H:%M:%S)] All jobs finished."

    # ── Merge with hadd ───────────────────────────────────────────────────────
    HADD_FORCE=""
    EXISTING=()
    for f in VBFHtoInv.slimfriend.root QCDHtoInv.slimfriend.root; do
        [ -f "$f" ] && EXISTING+=("$f")
    done

    if [ ${#EXISTING[@]} -gt 0 ]; then
        echo "WARNING: the following output files already exist:"
        for f in "${EXISTING[@]}"; do echo "  $f"; done
        read -r -p "Overwrite? [y/N] " answer
        case "$answer" in
            [yY]|[yY][eE][sS]) HADD_FORCE="-f" ;;
            *) echo "Aborted."; exit 1 ;;
        esac
    fi

    echo "Merging VBFHtoInv output files..."
    hadd $HADD_FORCE VBFHtoInv.slimfriend.root VBFHtoInv*_*.root

    echo "Merging QCDHtoInv output files..."
    hadd $HADD_FORCE QCDHtoInv.slimfriend.root QCDHtoInv*_*.root
fi

# ── Collect output files into friends/ ───────────────────────────────────────
mkdir -p friends
mv VBFHtoInv.slimfriend.root friends/
mv QCDHtoInv.slimfriend.root friends/
echo "Moved slimfriend files to friends/"

# ── Remove per-job intermediate files ────────────────────────────────────────
echo "Removing intermediate per-job files..."
rm -f *slimfriend_*.root

# ── Remove dynamically created jobs.sh ───────────────────────────────────────
rm -f jobs.sh
echo "Removed jobs.sh"

# ── Remove symlinks ───────────────────────────────────────────────────────────
for link in "${!ALL_SYMLINKS[@]}"; do
    if [ -L "$link" ]; then
        rm -f "$link"
        echo "Removed symlink: $link"
    fi
done

echo "[$(date +%H:%M:%S)] Done."
