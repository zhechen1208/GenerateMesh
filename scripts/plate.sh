#!/usr/bin/env bash
set -e
cd "$(dirname "$0")/.."

mkdir -p build output
gfortran -g -O0 -fcheck=all -fbacktrace -ffree-line-length-none \
    src/3d_plate.f90 -o build/3d_plate
./build/3d_plate
