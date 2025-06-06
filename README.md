# 大气压流光（Streamer）多物理量数值模拟器  


## 功能概览
| 模块 | 说明 |
| ---- | ---- |
| **流体动力学** | 轴对称 2-D，可选 OpenMP 并行（`THREAD_NUM`＝64） |
| **等离子体化学** | 电子、正、负离子及振动态共 100 + 粒子；电子碰撞/离子-中性/中性-中性反应≥180 条；RK4 自适应步长 |
| **电场求解** | CUDA-GPU 多重网格 + Red-Black SOR；三层递归（NR×NZ→½→¼） |
| **对流 & 扩散** | MUSCL-SuperBee 高分辨率格式 + 显式扩散 |
| **光电离** | 三项 Helmholtz 近似 |
| **耦合** | 单时间步统一推进：Poisson → 场 → 速度 → 对流/扩散 → 化学 |


## 源代码与数据文件结构
project-root/
├── streamer_solver.cu # 您粘贴的主程序（含 C/C++/CUDA/OpenMP）
├── include/ # vec(), mat(), spline() 等头文件（需自行提供）
├── inputdata/
│ ├── Initial/initial_N2_80p_300K.dat
│ ├── Bolsig_Data/O2_20p.dat
│ ├── i_reaction_300K_modmod0516.dat
│ ├── e_reaction_mod1803.dat
│ ├── n_reaction_modmod.dat
│ ├── mesh_r2_0624.dat
│ ├── mesh_z2.dat
│ └── V_Ono_single_str.dat
└── outputdata/ # 运行时自动生成> **提示**  


---

## 依赖环境
| 组件 | 建议版本 | 说明 |
| ---- | -------- | ---- |
| **CUDA Toolkit** | ≥ 11.4 | 代码调用 `cudaMalloc`, `dim3`, 需要 GPU (SM 5.0+) |
| **C/C++ 编译器** | GCC 10+ / Clang 12+ / MSVC 2019 | 必须支持 C++17 |
| **Make 或 CMake** | (二选一) | 本 README 同时给出两种脚本 |
| **OpenMP 运行时** | 可选 | 如果想启用 `THREAD_NUM`≥2 |
| **操作系统** | Linux x86-64/WSL2 | Windows+CUDA 亦可（实验性） |

---
