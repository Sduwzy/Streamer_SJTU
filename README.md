<!--
   README.md – Atmospheric-Pressure Streamer Solver
   =================================================
   Copy this file into the root of your repository.
-->

<h1 align="center">Atmospheric-Pressure Streamer Solver</h1>
<p align="center">
  <em>GPU-accelerated multi-physics solver for nanosecond discharge streamers in N<sub>2</sub>/O<sub>2</sub>/H<sub>2</sub>O mixtures</em>
</p>
<p align="center">
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue.svg" alt="License"></a>
  <img src="https://img.shields.io/github/languages/top/your-name/streamer-solver.svg" alt="top-lang">
  <img src="https://img.shields.io/github/workflow/status/your-name/streamer-solver/CI" alt="ci">
</p>

---

## ✨ Features

* **2-D axisymmetric fluid model** for electrons, 12 positive ions, 5 negative ions, and 22 vibrational levels  
* **Multi-rate chemistry ( > 180 reactions )** with 4-order Runge–Kutta integration  
* **Multigrid Poisson-SOR** solver implemented in **CUDA** (SM 5.0+)  
* **MUSCL-TVD** advection and explicit diffusion for gas dynamics (OpenMP/C++ threads)  
* Volumetric power deposition, cathode–anode optical/photo-ionisation, and user-defined voltage waveforms  
* Pinned-memory host–device transfers, tiling, and load balancing for **10–40× speed-up** versus CPU-only builds  

<p align="center"><img src="docs/demo.gif" width="650"></p>

---

## 📁 Directory Layout

