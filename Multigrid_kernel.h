#ifndef MULTIGRID_KERNEL_H
#define MULTIGRID_KERNEL_H

#include <cuda_runtime.h>
// 1) Poisson
void Poisson_GPU_function(
    dim3 dimGrid, dim3 dimBlock,
    double* d_phi, double* d_rho, double* d_rh, double* d_temp,
    int NX, int NY,
    double* d_P1, double* d_P2, double* d_P3, double* d_P4, double* d_P5,
    int* d_flag,
    double OMEGA, int itnum
);

// 2) Error of Poisson
void Error_poisson_GPU(
    dim3 dimGrid, dim3 dimBlock,
    double* d_phi, double* d_rho, double* d_rh,
    int NX, int NY,
    double* d_P1, double* d_P2, double* d_P3, double* d_P4, double* d_P5,
    int* d_flag,
    double* d_err
);

// 3) Restriction
void Ristriction_GPU(
    dim3 dimGrid2, dim3 dimBlock2,
    int NX2, int NY2, int NY,
    double* d_err, double* dd_rho, double* dd_Cphi
);

// 4) Interpolation
void Interporation_GPU(
    dim3 dimGrid, dim3 dimBlock,
    int NX, int NY, int NY2,
    int* d_flag, double* d_phi, double* dd_Cphi
);

// 5) Convergence Check
void Convergence_check_GPU(
    int N,
    double* d_temp, double* d_phi,
    double mf,
    double* temp, double* pphi,
    double* error, double* Maxphi
);

// 6) Error_Helm
void Error_Helm_GPU(
    dim3 dimGrid, dim3 dimBlock,
    double* d_phi, double* d_rho, double* d_rh,
    int NX, int NY,
    double* d_P1, double* d_P2, double* d_P3, double* d_P4, double* d_P5,
    int* d_flag,
    int* d_iflag,
    int* d_jflag,
    int* d_oflag,
    double* d_err,
    int pnum,
    double pO2
);

// 7) Helmholtz0
void Helmholtz0_GPU_function(
    dim3 dimGrid, dim3 dimBlock,
    double* d_phi, double* d_rho, double* d_rh, double* d_temp,
    int NX, int NY,
    double* d_P1, double* d_P2, double* d_P3, double* d_P4, double* d_P5,
    int* d_flag,
    int* d_iflag,
    int* d_jflag,
    int* d_oflag,
    double OMEGA,
    int itnum,
    double pO2
);

// 8) Helmholtz1
void Helmholtz1_GPU_function(dim3 dimGrid, dim3 dimBlock,
    double* d_phi, double* d_rho, double* d_rh,
    int NX, int NY,
    double* d_P1, double* d_P2, double* d_P3, double* d_P4, double* d_P5,
    int* d_flag,
    int* d_iflag,
    int* d_jflag,
    int* d_oflag,
    double* d_err,
    int pnum,
    double pO2);

// 9) Helmholtz2
void Helmholtz2_GPU_function(dim3 dimGrid, dim3 dimBlock,
    double* d_phi, double* d_rho, double* d_rh,
    int NX, int NY,
    double* d_P1, double* d_P2, double* d_P3, double* d_P4, double* d_P5,
    int* d_flag,
    int* d_iflag,
    int* d_jflag,
    int* d_oflag,
    double* d_err,
    int pnum,
    double pO2);

// 10) Multi_Helmholtz
void Multi_Helmholtz_GPU_function(dim3 dimGrid, dim3 dimBlock,
    double* d_phi, double* d_rho, double* d_rh,
    int NX, int NY,
    double* d_P1, double* d_P2, double* d_P3, double* d_P4, double* d_P5,
    int* d_flag,
    int* d_iflag,
    int* d_jflag,
    int* d_oflag,
    double* d_err,
    int pnum,
    double pO2);

#endif // MULTIGRID_KERNEL_H
