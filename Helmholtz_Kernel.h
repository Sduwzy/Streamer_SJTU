
#ifndef HELMHOLTZ_KERNEL_H
#define HELMHOLTZ_KERNEL_H

#include <cmath>            // 如果需要用到fabs之类


// 这里声明若干 __global__ 函数 (即 CUDA 内核)

__global__ void Helmholtz_Kernel0(
    double* d_phi,
    double* d_rho,
    double* d_rh,
    double* d_temp,
    int     RB_control,
    int     NX,
    int     NY,
    double* d_P1,
    double* d_P2,
    double* d_P3,
    double* d_P4,
    double* d_P5,
    int*    d_flag,
    int*    d_iflag,
    int*    d_jflag,
    int*    d_oflag,
    double  omega,
    double  ppO2
);

__global__ void Helmholtz_Kernel1(
    double* d_phi,
    double* d_rho,
    double* d_rh,
    double* d_temp,
    int     RB_control,
    int     NX,
    int     NY,
    double* d_P1,
    double* d_P2,
    double* d_P3,
    double* d_P4,
    double* d_P5,
    int*    d_flag,
    int*    d_iflag,
    int*    d_jflag,
    int*    d_oflag,
    double  omega,
    double  ppO2
);

__global__ void Helmholtz_Kernel2(
    double* d_phi,
    double* d_rho,
    double* d_rh,
    double* d_temp,
    int     RB_control,
    int     NX,
    int     NY,
    double* d_P1,
    double* d_P2,
    double* d_P3,
    double* d_P4,
    double* d_P5,
    int*    d_flag,
    int*    d_iflag,
    int*    d_jflag,
    int*    d_oflag,
    double  omega,
    double  ppO2
);

__global__ void Helm_Error_Kernel(
    double* d_phi,
    double* d_rho,
    int     NX,
    int     NY,
    double* d_rh,
    double* d_P1,
    double* d_P2,
    double* d_P3,
    double* d_P4,
    double* d_P5,
    int*    d_flag,
    int*    d_iflag,
    int*    d_jflag,
    int*    d_oflag,
    double* d_err,
    int     pnum,
    double  ppO2
);

__global__ void multi_Helmholtz_Kernel(
    double* d_phi,
    double* d_rho,
    double* d_rh,
    int     RB_control,
    int     NX,
    int     NY,
    double* d_P1,
    double* d_P2,
    double* d_P3,
    double* d_P4,
    double* d_P5,
    int*    d_flag,
    int*    d_iflag,
    int*    d_jflag,
    int*    d_oflag,
    double  omega
);

#endif // HELMHOLTZ_KERNEL_H
