#!/bin/bash

# Function to display help
usage() {
    echo "Usage: $0 -f FILE -n nJOBS -t TOT -x XS"
    echo ""
    echo "Arguments:"
    echo "  -f FILE     Path to the input file"
    echo "  -n nJOBS    Number of jobs"
    echo "  -t TOT      Total events per job"
    echo "  -x XS       Cross-section value"
    echo "  --help      Display this help message"

    echo "A concrete example:"
    echo "./jobs.sh -f ~/files/akcolor/Spring2024/LHC-Matchbox-VBFH-3.root -n 3 -t 2 -x 3.601"
    exit 1
}

# Default values (in case user does not provide them)
FILE=""
nJOBS=0
TOT=0
XS=0.0

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -f|--file) FILE="$2"; shift ;;
        -n|--njobs) nJOBS="$2"; shift ;;
        -t|--tot) TOT="$2"; shift ;;
        -x|--xs) XS="$2"; shift ;;
        --help) usage ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Ensure that required arguments are provided
if [[ -z "$FILE" || -z "$nJOBS" || -z "$TOT" || -z "$XS" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# Main loop for processing jobs
for (( JOB=0; JOB<nJOBS; JOB++ )); do
    SKIP=$((JOB * TOT))
    python makefriends.py "$FILE" --genWeight evweight --xs $XS --skip $SKIP --totEve $TOT --jobId $JOB &
done


#TOT=10
#JOB=0;SKIP=$((JOB*$TOT)); python makefriends.py $ifile  --genWeight evweight --xs 3.601 --skip $SKIP  --totEve $TOT --jobId $JOB
#JOB=1;SKIP=$((JOB*$TOT)); python makefriends.py $ifile  --genWeight evweight --xs 3.601 --skip $SKIP  --totEve $TOT --jobId $JOB
