#ifndef POISSON_H
#define POISSON_H

#include <cstdio>   // 如果需要 printf/fscanf 等
#include <cmath>    // 如果需要 sqrt/fabs 等
#include "memory.h" // 假设您的 mat(), free_mat() 等在此声明

// ----------------------------------------------------------------------
// 在此声明您需要的函数 (来自原 .h):
// ----------------------------------------------------------------------

// 说明：discretization 函数主要做离散化处理
//  - N, M: 网格数量
//  - rhalf, zhalf: 网格坐标
//  - P1,P2,P3,P4,P5: 输出用的系数矩阵
//  - a,b: 与几何或坐标变换相关的长度参数
//  - iflag, jflag, otuflag, flag: 标记数组
void discretization(int N, int M,
                    double* rhalf, double* zhalf,
                    double** P1, double** P2, double** P3, double** P4, double** P5,
                    double a, double b,
                    int** iflag, int** jflag, int** otuflag, int** flag);


// 说明：poiseq 函数使用 SOR(或类似)方法求解 Poisson
//  - num_x, num_y: 网格大小
//  - pphi: 电势分布
//  - flag, iflag, jflag, otuflag: 标志位数组
//  - rhalf, zhalf: 网格坐标
//  - rho: 电荷密度
//  - P1..P5: 离散后系数矩阵
void poiseq(int num_x, int num_y,
            double** pphi, int** flag,
            int** iflag, int** jflag, int** otuflag,
            double* rhalf, double* zhalf,
            double** rho,
            double** P1, double** P2, double** P3, double** P4, double** P5);

#endif // POISSON_H
