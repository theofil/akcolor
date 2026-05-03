#!/bin/zsh
PYTHON=/Users/konstantinostheofilatos/anaconda3/envs/root/bin/python
LAST=$(ls -t *.py | head -1)
echo "Running: $LAST"
$PYTHON "$LAST"
