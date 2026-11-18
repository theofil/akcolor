#!/bin/bash
# condor wrapper: train_fold.py with the job arguments, e.g. --net NNjPF --fold 0
# [--campaign fall26]: train one fold model, then score its inference sub-folds
source /cvmfs/sft.cern.ch/lcg/views/LCG_106_cuda/x86_64-el9-gcc11-opt/setup.sh
cd /eos/home-t/theofil/work/akcolor/makefriends/C/nov26
python3 train_fold.py "$@"
