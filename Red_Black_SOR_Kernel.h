#ifndef RED_BLACK_SOR_KERNEL_H
#define RED_BLACK_SOR_KERNEL_H

// 如果需要 CUDA 运行时头文件，可在这里包含
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>

// 如果在别处有定义 point(...) 宏，可删除这里的重复定义
#ifndef point
#define point(N, i, j) ( (N)*(i) + (j) )
#endif



//-------------------------------------------------
// 1) Red_Black_SOR_Kernel
//-------------------------------------------------
__global__ void Red_Black_SOR_Kernel(
    double *d_phi,
    double *d_rho,
    double *d_rh,
    double *d_temp,
    int RB_control,
    int NX,
    int NY,
    double *d_P1,
    double *d_P2,
    double *d_P3,
    double *d_P4,
    double *d_P5,
    int *d_flag,
    double omega
);

//-------------------------------------------------
// 2) Error_Kernel
//-------------------------------------------------
__global__ void Error_Kernel(
    double *d_phi,
    double *d_rho,
    double *d_rh,
    int NX,
    int NY,
    double *d_P1,
    double *d_P2,
    double *d_P3,
    double *d_P4,
    double *d_P5,
    int *d_flag,
    double *d_err
);

//-------------------------------------------------
// 3) Restriction_Kernel
//-------------------------------------------------
__global__ void Restriction_Kernel(
    int NX2,
    int NY2,
    int NY,
    double *d_Cres,
    double *d_Cres2,
    double *d_cCphi
);

//-------------------------------------------------
// 4) Interporation_kernel
//-------------------------------------------------
__global__ void Interporation_kernel(
    int NX,
    int NY,
    int NY2,
    int *d_flag,
    double *d_phi,
    double *dd_phi
);

#endif // RED_BLACK_SOR_KERNEL_H
