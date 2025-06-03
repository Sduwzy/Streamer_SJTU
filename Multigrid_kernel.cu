#include "Multigrid_kernel.h"
#include "Red_Black_SOR_Kernel.h"
#include "Helmholtz_Kernel.h"
#include <cuda_runtime.h>
#include <cstdio>  // for printf, etc.
#include <cmath>   // if needed
#include <cstdlib> // if needed


// 1) Poisson_GPU_function
void Poisson_GPU_function(
    dim3 dimGrid, dim3 dimBlock,
    double *d_phi, double *d_rho, double *d_rh, double *d_temp,
    int NX, int NY,
    double *d_P1, double *d_P2, double *d_P3, double *d_P4, double *d_P5,
    int *d_flag,
    double OMEGA,
    int itnum)
{
    // 原示例:
    //  Gauss-Siedel-SOR
    for (int it = 0; it < itnum; ++it)
    {

        // odd sites (Red Point)
        Red_Black_SOR_Kernel<<<dimGrid, dimBlock>>>(
            d_phi, d_rho, d_rh, d_temp,
            0, NX, NY,
            d_P1, d_P2, d_P3, d_P4, d_P5,
            d_flag,
            OMEGA);
        cudaThreadSynchronize(); // or cudaDeviceSynchronize()

        // even sites (Black Point)
        Red_Black_SOR_Kernel<<<dimGrid, dimBlock>>>(
            d_phi, d_rho, d_rh, d_temp,
            1, NX, NY,
            d_P1, d_P2, d_P3, d_P4, d_P5,
            d_flag,
            OMEGA);
        cudaThreadSynchronize(); // or cudaDeviceSynchronize()
    }
}

// 2) Error_poisson_GPU
void Error_poisson_GPU(
    dim3 dimGrid, dim3 dimBlock,
    double *d_phi, double *d_rho, double *d_rh,
    int NX, int NY,
    double *d_P1, double *d_P2, double *d_P3, double *d_P4, double *d_P5,
    int *d_flag,
    double *d_err)
{
    Error_Kernel<<<dimGrid, dimBlock>>>(
        d_phi, d_rho, d_rh,
        NX, NY,
        d_P1, d_P2, d_P3, d_P4, d_P5,
        d_flag,
        d_err);
    cudaThreadSynchronize(); // or cudaDeviceSynchronize()
}

// 3) Ristriction_GPU
void Ristriction_GPU(
    dim3 dimGrid2, dim3 dimBlock2,
    int NX2, int NY2, int NY,
    double *d_err, double *dd_rho, double *dd_Cphi)
{
    Restriction_Kernel<<<dimGrid2, dimBlock2>>>(
        NX2, NY2, NY,
        d_err, dd_rho, dd_Cphi);
    cudaThreadSynchronize();
}

// 4) Interporation_GPU
void Interporation_GPU(
    dim3 dimGrid, dim3 dimBlock,
    int NX, int NY, int NY2,
    int *d_flag, double *d_phi, double *dd_Cphi)
{
    Interporation_kernel<<<dimGrid, dimBlock>>>(
        NX, NY, NY2,
        d_flag, d_phi, dd_Cphi);
    cudaThreadSynchronize();
}

// 5) Convergence_check_GPU
void Convergence_check_GPU(
    int N, double *d_temp, double *d_phi,
    double mf,
    double *temp, double *pphi,
    double *error, double *Maxphi)
{
    // 原示例：
    cudaMemcpy(temp, d_temp, static_cast<size_t>(mf), cudaMemcpyDeviceToHost);
    cudaMemcpy(pphi, d_phi, static_cast<size_t>(mf), cudaMemcpyDeviceToHost);

    (*error) = 0.0;
    (*Maxphi) = 0.0;

    for (int i = 0; i < N; i++)
    {
        if (temp[i] > (*error))
        {
            (*error) = temp[i];
        }
        if (std::fabs(pphi[i]) > (*Maxphi))
        {
            (*Maxphi) = std::fabs(pphi[i]);
        }
    }
}

// 6) Error_Helm_GPU
void Error_Helm_GPU(
    dim3 dimGrid, dim3 dimBlock,
    double *d_phi, double *d_rho, double *d_rh,
    int NX, int NY,
    double *d_P1, double *d_P2, double *d_P3, double *d_P4, double *d_P5,
    int *d_flag,
    int *d_iflag,
    int *d_jflag,
    int *d_oflag,
    double *d_err,
    int pnum,
    double pO2)
{
    Helm_Error_Kernel<<<dimGrid, dimBlock>>>(
        d_phi, d_rho,
        NX, NY,
        d_rh,
        d_P1, d_P2, d_P3, d_P4, d_P5,
        d_flag, d_iflag, d_jflag, d_oflag,
        d_err,
        pnum, pO2);
    cudaThreadSynchronize();
}

// 7) Helmholtz0_GPU_function
void Helmholtz0_GPU_function(
    dim3 dimGrid, dim3 dimBlock,
    double *d_phi, double *d_rho, double *d_rh, double *d_temp,
    int NX, int NY,
    double *d_P1, double *d_P2, double *d_P3, double *d_P4, double *d_P5,
    int *d_flag,
    int *d_iflag,
    int *d_jflag,
    int *d_oflag,
    double OMEGA,
    int itnum,
    double pO2)
{
    for (int it = 0; it < itnum; ++it)
    {
        // color=1
        Helmholtz_Kernel0<<<dimGrid, dimBlock>>>(
            d_phi, d_rho, d_rh, d_temp,
            1, NX, NY,
            d_P1, d_P2, d_P3, d_P4, d_P5,
            d_flag, d_iflag, d_jflag, d_oflag,
            OMEGA, pO2);
        cudaThreadSynchronize();

        // color=0
        Helmholtz_Kernel0<<<dimGrid, dimBlock>>>(
            d_phi, d_rho, d_rh, d_temp,
            0, NX, NY,
            d_P1, d_P2, d_P3, d_P4, d_P5,
            d_flag, d_iflag, d_jflag, d_oflag,
            OMEGA, pO2);
        cudaThreadSynchronize();
    }
}

// 8) Helmholtz1_GPU_function
void Helmholtz1_GPU_function(
    /* 同理 */
)
{
    // 类似上面的写法
}

// 9) Helmholtz2_GPU_function
void Helmholtz2_GPU_function(...) { /*同理*/ }

// 10) Multi_Helmholtz_GPU_function
void Multi_Helmholtz_GPU_function(
    dim3 dimGrid, dim3 dimBlock,
    double *d_phi, double *d_rho, double *d_rh, double *d_temp,
    int NX, int NY,
    double *d_P1, double *d_P2, double *d_P3, double *d_P4, double *d_P5,
    int *d_flag,
    int *d_iflag,
    int *d_jflag,
    int *d_oflag,
    double OMEGA,
    int itnum)
{
    for (int it = 0; it < itnum; ++it)
    {
        multi_Helmholtz_Kernel<<<dimGrid, dimBlock>>>(
            d_phi, d_rho, d_rh,
            0, NX, NY,
            d_P1, d_P2, d_P3, d_P4, d_P5,
            d_flag, d_iflag, d_jflag, d_oflag,
            1.8);
        cudaThreadSynchronize();

        multi_Helmholtz_Kernel<<<dimGrid, dimBlock>>>(
            d_phi, d_rho, d_rh,
            1, NX, NY,
            d_P1, d_P2, d_P3, d_P4, d_P5,
            d_flag, d_iflag, d_jflag, d_oflag,
            1.8);
        cudaThreadSynchronize();
    }
}
