#include <cuda_runtime.h>  // 常见的头文件
#include <device_launch_parameters.h>
#include <cmath>

#include "Helmholtz_Kernel.h"   // 包含我们刚才的声明


//====================== Kernel 0 =======================
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
    double  ppO2)
{
    // ... （这里直接粘贴你原先 Helmholtz_Kernel0 的实现） ...
    // 注意，把 point(NY, i, j) 之类的宏写好 or ensure we have #define point(...
    // 也注意 i, j 的取值别超范围
    // ...
}

//====================== Kernel 1 =======================
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
    double  ppO2)
{
    // 原先 Helmholtz_Kernel1(...) 函数的内容
    // ...
}

//====================== Kernel 2 =======================
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
    double  ppO2)
{
    // 原先 Helmholtz_Kernel2(...) 函数的内容
    // ...
}

//====================== Helm_Error_Kernel =======================
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
    double  ppO2)
{
    // 原先 Helm_Error_Kernel(...) 的内容
    // ...
}

//====================== multi_Helmholtz_Kernel =======================
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
    double  omega)
{
    // 原先 multi_Helmholtz_Kernel(...) 的内容
    // ...
}
