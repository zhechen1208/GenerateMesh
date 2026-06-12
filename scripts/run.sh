#!/usr/bin/env bash
# === GenerateMesh 全流程 ===
#   1. 生成翼型内边界 → params.inc 自动更新 kleft/kright
#   2. 编译运行 mesh → 生成体网格
#   3. 合并为 Nektar++ XML
set -e
cd "$(dirname "$0")/.."
mkdir -p output

NPROC=${NPROC:-8}

# ---- 1. beta ----
echo "=== [1/3] 生成翼型内边界 ==="
if [[ ! -f params.inc ]]; then
    cp params.inc.template params.inc
    echo "已从模板创建 params.inc，请检查参数后重新运行"
    exit 1
fi

make beta
./build/beta_function_wing 2>&1 | tee build/beta_output.log

KLEFT=$(grep -oP 'kleft=\s*\K[0-9]+' build/beta_output.log | tail -1)
KRIGHT=$(grep -oP 'kright=\s*\K[0-9]+' build/beta_output.log | tail -1)
[[ -z "$KLEFT" || -z "$KRIGHT" ]] && { echo "Error: 未能解析 kleft/kright" >&2; exit 1; }

sed -i -E "s/kleft=[0-9]+,kright=[0-9]+/kleft=${KLEFT},kright=${KRIGHT}/" params.inc
echo "params.inc: kleft=${KLEFT} kright=${KRIGHT}"

# ---- 2. mesh ----
echo "=== [2/3] 编译并运行 mesh ==="
make mesh
mpirun -np "$NPROC" ./build/mesh

# ---- 3. xml ----
echo "=== [3/3] 合并为 Nektar++ XML ==="
python scripts/dat2xml.py \
    output/3d_mesh.dat \
    output/left_tip_block.dat \
    output/right_tip_block.dat \
    -o output/mesh.xml \
    --kleft "$KLEFT" --kright "$KRIGHT" \
    --num-modes ${NUM_MODES:-2} --num-points ${NUM_POINTS:-3} \
    --no-fix-negative

echo "=== 完成: output/mesh.xml ==="
