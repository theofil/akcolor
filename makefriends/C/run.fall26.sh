#!/bin/bash
# Fall2026 friends (see fall26.md): one makefriends_fall26 process per input chunk, merged per
# sample in chunk order, then converted to HDF5.
#
#   ./run.fall26.sh                  # full production -> friends/fall26/<sample>.{friend.root,h5}
#   ./run.fall26.sh --goFast 2000    # first chunk of each sample, first 2000 events -> fast_tmp_fall26/
#
#   PROCESSES="VBFZ QCDZjj" (default)  processes to run, both generators each
#   NPROC=40 (default)                 makefriends_fall26 processes at a time
#
# Prerequisite: `source init.source` (LCG ROOT/FastJet, hadd, python3) and `make makefriends_fall26`.
# The chunk list, offsets and weight factors come from fall26/manifest.py (Fall2026 manifest.json).
# Outputs are always deleted before they are rewritten: overwriting a file in place on EOS
# leaves hidden .sys.v#.<name> version dirs behind.
set -euo pipefail

WORKDIR=$(cd "$(dirname "$0")" && pwd)
PROCESSES=${PROCESSES:-"VBFZ QCDZjj"}
NPROC=${NPROC:-40}
GOFAST_N=""
if [[ "${1:-}" == "--goFast" ]]; then
    GOFAST_N=${2:?"--goFast needs a number of events"}
fi

[[ -x "$WORKDIR/makefriends_fall26" ]] || { echo "build it first: make makefriends_fall26"; exit 1; }
command -v hadd >/dev/null || { echo "hadd not found: source init.source"; exit 1; }

if [[ -n "$GOFAST_N" ]]; then
    OUTDIR="$WORKDIR/fast_tmp_fall26"
else
    OUTDIR="$WORKDIR/friends/fall26"
fi
TMPDIR_F26="$OUTDIR/tmp"
rm -rf "$TMPDIR_F26"
mkdir -p "$TMPDIR_F26"

# job lines: sample path pb_per_raw_weight chunk evOffset boson charge events
# shellcheck disable=SC2086
python3 "$WORKDIR/fall26/manifest.py" jobs $PROCESSES > "$TMPDIR_F26/jobs.all"
if [[ -n "$GOFAST_N" ]]; then
    awk '!seen[$1]++' "$TMPDIR_F26/jobs.all" > "$TMPDIR_F26/jobs.txt"   # first chunk per sample
else
    cp "$TMPDIR_F26/jobs.all" "$TMPDIR_F26/jobs.txt"
fi
samples=$(awk '{print $1}' "$TMPDIR_F26/jobs.txt" | uniq)
for s in $samples; do mkdir -p "$TMPDIR_F26/$s"; done
echo "[$(date +%H:%M:%S)] $(wc -l < "$TMPDIR_F26/jobs.txt") chunk jobs, $NPROC at a time, output in $OUTDIR"

run_one() {   # sample path pb_per_raw_weight chunk evOffset boson charge events
    local sample=$1 path=$2 factor=$3 chunk=$4 off=$5 boson=$6 q=$7
    local base; base=$(basename "$path" .root)
    local out="$TMPDIR_F26/$sample/$base.friend.root" log="$TMPDIR_F26/$sample/$base.log"
    local extra=()
    [[ "$q" != 0 ]] && extra+=(--charge "$q")
    [[ -n "$GOFAST_N" ]] && extra+=(--totEve "$GOFAST_N")
    rm -f "$out"
    if "$WORKDIR/makefriends_fall26" "$path" --pbPerRawWeight "$factor" --boson "$boson" \
            --chunk "$chunk" --evOffset "$off" "${extra[@]}" --output "$out" > "$log" 2>&1; then
        echo "[$(date +%H:%M:%S)] ok     $sample $base"
    else
        echo "[$(date +%H:%M:%S)] FAILED $sample $base (log: $log)"
        return 1
    fi
}
export -f run_one
export WORKDIR TMPDIR_F26 GOFAST_N

if ! xargs -P "$NPROC" -L 1 bash -c 'run_one "$@"' _ < "$TMPDIR_F26/jobs.txt"; then
    echo "[ERROR] at least one chunk failed; the logs are in $TMPDIR_F26/<sample>/"
    exit 1
fi

# merge each sample in job (= offset) order, so friend entry i is dataset row i (evIdx checks it)
out_files=()
for s in $samples; do
    parts=()
    while read -r sample path _; do
        [[ "$sample" == "$s" ]] && parts+=("$TMPDIR_F26/$s/$(basename "$path" .root).friend.root")
    done < "$TMPDIR_F26/jobs.txt"
    rm -f "$OUTDIR/$s.friend.root" "$OUTDIR/$s.h5"
    if (( ${#parts[@]} == 1 )); then
        mv "${parts[0]}" "$OUTDIR/$s.friend.root"
    else
        hadd -f "$OUTDIR/$s.friend.root" "${parts[@]}" > "$TMPDIR_F26/hadd_$s.log"
    fi
    echo "[$(date +%H:%M:%S)] merged ${#parts[@]} chunks -> $OUTDIR/$s.friend.root"
    out_files+=("$OUTDIR/$s.friend.root")
done

echo "[$(date +%H:%M:%S)] Converting ${#out_files[@]} friend.root files to HDF5..."
python3 "$WORKDIR/root_to_hdf5_fall26.py" "${out_files[@]}"

cat "$TMPDIR_F26"/*/*.log > "$OUTDIR/run.fall26.chunks.log"
rm -rf "$TMPDIR_F26"
echo "[$(date +%H:%M:%S)] Done."
