#!/usr/bin/env bash
# Build & run mesh
set -e
cd "$(dirname "$0")/.."

make mesh && mpirun -np ${NPROC:-8} ./build/mesh
