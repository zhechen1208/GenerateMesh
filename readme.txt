# GenerateMesh — 三维 O-H 网格生成

针对 **beta-function wing** 的三维结构化 O-H 网格生成工具链。

## 目录结构

    src/          Fortran 源代码
    scripts/      Shell 和 Python 脚本
    build/        编译产物（.o, .mod, 可执行文件）
    output/       输入输出数据文件
    params.inc    网格参数（编译前配置）

## 快速开始

### 1. 配置参数

编辑自动生成的 params.inc，或运行：

    bash scripts/gen_params.sh 41 41 91 24 68 5 4.0 40.0
    #                          imax jmax kmax kleft kright iarc span domain

### 2. 编译

    make          # release（-O2）
    make debug    # debug（-g -O0 -fcheck=all）

### 3. 生成翼型内边界

    bash scripts/beta.sh       # beta-function wing（推荐）
    bash scripts/plate.sh      # 平板（备选）

按提示输入参数，生成 output/innerboundary.dat，屏幕输出 kleft/kright。

### 4. 运行网格生成

    mpirun -np 8 ./build/mesh

输出在 output/ 下：3d_mesh.dat, left_tip_block.dat, right_tip_block.dat 等。

### 5. 可视化

    bash scripts/trans.sh      # dat → xml → plt (19MB)

或直接用 Tecplot 打开 output/3d_mesh.dat。

## 一键全流程（平板）

    bash scripts/mesh-plate.sh

## 常用命令

    make clean       # 清理 .o .mod
    make distclean   # 清理全部（含 build/ 和 params.inc）
    make && mpirun -np 8 ./build/mesh    # 编译+运行
