#!/bin/bash
# pilot: 4 seeds per case, N=1e5, 100 inactive + 100 active; FD 4 seeds h=0.1 cm
cd "$(dirname "$0")/.."
PY=../.venv/bin/python
for c in 1 2 3; do
  $PY run3d.py sens $c 4 100000 100 100 1
done
$PY run3d.py fd 1 4 100000 100 100 0.1 R 1001
$PY run3d.py fd 2 4 100000 100 100 0.1 a,b 1001
$PY run3d.py fd 3 4 100000 100 100 0.1 a,tx,c,zhi 1001
