#!/bin/bash
# Step 0b long runs: 100 seeds per case (N=1e5, 100 inactive + 200 active, LMAX=12),
# then CRN central differences, 100 seeds per parameter, h = 0.1 cm.
cd "$(dirname "$0")"
PY=../.venv/bin/python
for c in 1 2 3; do
  $PY run3d.py sens $c 100 100000 100 200 1 >> runs_sens_log.txt 2>&1
done
$PY run3d.py fd 1 100 100000 100 200 0.1 R 1001 >> runs_fd_log.txt 2>&1
$PY run3d.py fd 2 100 100000 100 200 0.1 a,b 1001 >> runs_fd_log.txt 2>&1
$PY run3d.py fd 3 100 100000 100 200 0.1 a,tx,c,zhi 1001 >> runs_fd_log.txt 2>&1
echo ALLDONE >> runs_fd_log.txt
