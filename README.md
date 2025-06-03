<!-- ================================================================
     Atmospheric-Pressure Streamer Solver – README
     ================================================================

     Copy / rename this file to  `README.md`  in the root of the repo.
     Markdown is GitHub-flavoured (GFM) and renders correctly on GitHub,
     GitLab, Bitbucket, VS Code and most static-site generators.
     ----------------------------------------------------------------
-->

<h1 align="center">Atmospheric-Pressure Streamer Solver</h1>
<p align="center">
  <em>A GPU-accelerated multi-physics code for nanosecond streamer discharges<br>
  in humid air (N<sub>2</sub>/O<sub>2</sub>/H<sub>2</sub>O mixtures)</em>
</p>

<p align="center">
  <!-- update these badges for your repo/CI provider -->
  <a href="LICENSE"><img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
  <img alt="Top language"  src="https://img.shields.io/github/languages/top/YOUR-USER/streamer-solver.svg">
  <img alt="GitHub CI"      src="https://img.shields.io/github/actions/workflow/status/YOUR-USER/streamer-solver/ci.yml?branch=main">
</p>

---

## ✨ Key capabilities

| Module | Highlights |
|--------|------------|
| **Fluid plasma model** | 2-D axisymmetric, 1-T, ⬆️ 20 charged species & 22 vibrational levels |
| **Chemistry** | 180 + reactions (electron impact, ion–ion, ion–neutral, neutral–neutral)<br>4th-order RK (adaptive) |
| **Electric field** | Multigrid Poisson solver with red–black SOR on **CUDA** (sm ≥ 5.0) |
| **Advection / diffusion** | MUSCL-TVD (super-bee limiter) + explicit diffusion (OpenMP / C++ threads) |
| **Photo-ionisation** | Helmholtz approximation with three-term absorption kernel |
| **User I/O** | Flexible voltage waveform spline, VTK/CSV/HDF5 output, runtime checkpointing |
| **Performance** | 10–40× faster than baseline Fortran build on RTX 3080 (8 M cells / 30 ns) |

<p align="center"><img width="650" src="docs/demo.gif" alt="Streamer evolution demo"></p>

---

## 🗂️ Repository layout

├── cmake/ # CMake helpers
├── include/ # public headers (vec.h, spline.h, ...)
├── src/ # C++ / CUDA source files
│ ├── chemistry/ # reaction integrators
│ ├── field/ # Poisson + multigrid kernels
│ └── fluid/ # gas-dynamic solver
├── inputdata/ # ✓ example meshes, reaction tables, voltage waveform
├── docs/ # build + usage guides, figures, publications
├── tests/ # unit tests & CI regression cases
├── examples/ # ready-to-run example configurations
├── .github/workflows/ # CI definitions (build / clang-tidy / unit tests)
├── CMakeLists.txt # build script (GNU/Clang/MSVC/NVIDIA/PGI)
└── README.md # this file


---

## 🚀 Quick start

<details>
<summary>Prerequisites (tested configurations)</summary>

| Dependency | Recommended | Notes |
|------------|-------------|-------|
| **CUDA Toolkit** | ≥ 11.4 | code uses cooperative groups & `__shfl_sync` |
| **C++ compiler** | GCC ≥10 · Clang ≥12 · MSVC 2019 | must support C++17 |
| **CMake** | ≥ 3.18 | presets available |
| **GPU** | Compute capability ≥ 5.0, ≥ 4 GB VRAM | Pascal, Volta, Turing, Ampere, Ada OK |
| **CPU** | Any modern x86-64 | OpenMP 4.5 runtime optional |
| **Linux / WSL 2** | Ubuntu 20.04 LTS+ | Windows & macOS (CUDA on eGPU) experimental |
</details>

```bash
# Clone
git clone https://github.com/YOUR-USER/streamer-solver.git
cd streamer-solver

# Configure + build (Release)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

build/solver                               \
  --mesh-r       inputdata/mesh_r2_0624.dat \
  --mesh-z       inputdata/mesh_z2.dat      \
  --initial      inputdata/Initial/initial_N2_80p_300K.dat \
  --voltage      inputdata/V_Ono_single_str.dat \
  --t-stop       8e-9                       \
  --output       outputdata/run-0001

outputdata/run-0001/
 ├── φ_0008000.vti          # electric potential (ParaView)
 ├── n_electron_0008000.vti # species densities, power deposition, …
 ├── discharge.log          # human-readable summary
 └── chkpt_0008000.h5       # restart checkpoint (HDF5)

 ┌──────────────────────────────────────────────────────────┐
 │                          main()                         │
 └──────────────────────────────────────────────────────────┘
              │
              ▼
   ┌─────────────────────┐      Init meshes, load reactions,
   │ initial_condition() │◄──── constants, allocate buffers
   └─────────────────────┘
              │
              ▼
   ┌─────────────────────┐
   │   time-loop (ns)    │   nstp = 0 … N
   └─────────────────────┘
        │  │     ▲  ▲
        │  │     │  │ RK4 chemistry (CPU / OpenMP)
        │  │     │  └─────────┐
        │  │     │            │ diffusion   (CPU)
        │  │     └────────────┤ advection   (CPU threads)
        │  │                  │ e/ion vel.  (CPU)
        │  └─►  Poisson-MG-SOR│ field solve (GPU, CUDA)
        │                     └──> Ex,Ey,φ,|E|
        ▼
  write VTK / progress bar … (I/O thread)

