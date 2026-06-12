#!/usr/bin/env bash
# Quick visualization: convert 3d_mesh.dat -> mesh.xml -> mesh.plt
set -e
cd "$(dirname "$0")/.."

python scripts/dat2xml.py output/3d_mesh.dat output/left_tip_block.dat output/right_tip_block.dat \
    -o output/mesh.xml --kleft 24 --kright 68 --no-fix-negative \
    --num-modes 2 --num-points 3

FieldConvert output/mesh.xml output/mesh.plt
echo "Done → output/mesh.plt"
