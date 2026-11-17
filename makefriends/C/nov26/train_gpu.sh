#!/bin/bash
# condor wrapper: train fold model $2 of net $1 (NNjPF | ParTj), then score its inference sub-folds
source /cvmfs/sft.cern.ch/lcg/views/LCG_106_cuda/x86_64-el9-gcc11-opt/setup.sh
cd /eos/home-t/theofil/work/akcolor/makefriends/C/nov26
python3 train_fold.py --net "$1" --fold "$2"
