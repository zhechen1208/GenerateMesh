# GenerateMesh — 三维 O-H 网格生成

针对 **beta-function wing** 的三维结构化 O-H 网格生成工具链。

## 目录结构

    src/          Fortran 源代码
    scripts/      Shell / Python 脚本
    build/        编译产物
    output/       输入输出数据

## 一键运行

    cp params.inc.template params.inc   # 首次配置：编辑网格大小、迭代参数
    vim beta_params.txt                 # 首次配置：编辑翼型参数
    bash scripts/run.sh                 # beta → mesh → xml 全流程

## 分步执行

    make beta && ./build/beta_function_wing    # 生成翼型内边界
    make mesh && mpirun -np 8 ./build/mesh     # 编译运行 mesh
    python scripts/dat2xml.py \                # 合并为 Nektar++ XML
        output/3d_mesh.dat output/left_tip_block.dat output/right_tip_block.dat \
        -o output/mesh.xml --kleft 24 --kright 68 --no-fix-negative

## 常用命令

    make              # 编译全部（mesh + beta）
    make mesh         # 只编译 mesh
    make beta         # 只编译 beta
    make debug        # debug 版（-g -O0 -fcheck=all）
    make clean        # 清理编译产物
    bash scripts/plate.sh  # 平板内边界（备选）

## XML 中 COMPOSITE 含义

| C ID | 标签 | 说明 |
|------|------|------|
| C[0] | `wall` | 主 block 机翼物面 (i=0, k 在 kleft..kright) |
| C[1] | `main_left_transition_wall` | 左过渡段物面 (i=0, k ≤ kleft) |
| C[2] | `main_right_transition_wall` | 右过渡段物面 (i=0, k ≥ kright) |
| C[3] | `left_tip_wall` | 左 tip block 物面 |
| C[4] | `right_tip_wall` | 右 tip block 物面 |
| C[5] | `farfield` | 远场 (i=imax) |
| C[6] | `kmin` | 展向左端面 (k=0) |
| C[7] | `kmax` | 展向右端面 (k=kmax) |
| C[8] | `jmin` | 周向前端面 (j=0) |
| C[9] | `jmax` | 周向后端面 (j=jmax) |
| C[10] | `fluid` | 流体域 (所有 hex) |
