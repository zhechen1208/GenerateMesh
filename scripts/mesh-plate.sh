#!/usr/bin/env bash
# Full pipeline: 3d_plate → mesh → visualization
set -e
cd "$(dirname "$0")/.."
mkdir -p build output

# ---- user parameters (edit these) ----
J_IN=61
K_WING=81
K_TRANS=21
AR=4
DOMAIN=10
THICKNESS=0.05
THETA_N=10
DZ_WING_TIP=0.045
DZ_TRANS_TIP=0.040
IMAX=31
IARC=10
NPROC=8
# -----------------------------------

# 1. build & run 3d_plate
gfortran -g -O0 -fcheck=all -fbacktrace -ffree-line-length-none \
    src/3d_plate.f90 -o build/3d_plate

cat > build/plate_input.txt << EOF
${J_IN}
${K_WING}
${K_TRANS}
${AR}
${DOMAIN}
${THICKNESS}
${THETA_N}
${DZ_WING_TIP}
${DZ_TRANS_TIP}
EOF

./build/3d_plate < build/plate_input.txt | tee build/plate_output.log

# 2. parse output
JMAX=$(grep -oP 'jmax=\s*\K[0-9]+' build/plate_output.log | tail -1)
KMAX_P=$(grep -oP 'kmax=\s*\K[0-9]+' build/plate_output.log | tail -1)
KLEFT=$(grep -oP 'kleft=\s*\K[0-9]+' build/plate_output.log | tail -1)
KRIGHT=$(grep -oP 'kright=\s*\K[0-9]+' build/plate_output.log | tail -1)

[[ -z "$JMAX" || -z "$KMAX_P" || -z "$KLEFT" || -z "$KRIGHT" ]] && {
    echo "Error: failed to parse plate output" >&2; exit 1; }
echo "Parsed: JMAX=$JMAX KMAX=$KMAX_P KLEFT=$KLEFT KRIGHT=$KRIGHT"

# 3. generate params.inc & build mesh
bash scripts/gen_params.sh $IMAX $JMAX $KMAX_P $KLEFT $KRIGHT $IARC $AR $DOMAIN
make && mpirun -np $NPROC ./build/mesh | tee build/main_output.log

echo "Done → output/3d_mesh.dat"
