#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_106_cuda/x86_64-el9-gcc11-opt/setup.sh
module load lxbatch/eossubmit

cd /eos/home-t/theofil/work/akcolor/makefriends/C/testAUC

echo "=== phi_visualize.py ===" && python phi_visualize.py
echo "=== phi_gradient.py ===" && python phi_gradient.py
echo "=== phi_pysr.py ===" && python phi_pysr.py
echo "=== observable_roc.py ===" && python observable_roc.py
echo "=== done ==="
