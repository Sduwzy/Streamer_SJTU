#ifndef ENERGY_TYPE_H
#define ENERGY_TYPE_H

#include <cstdio>   // std::fopen, std::fscanf, etc.
#include <cmath>    // std::exp, std::log
#include <cstdlib>  // std::exit
#include "memory.h" // 假设您已有 memory.h (vec, mat, free_vec, free_mat等)

// ============================================================================
// 函数声明
// ============================================================================

// 1) 将 O2、N2、H2O 三个数组合并到 M[] (每种 Q=8，共9项：索引0~8)
void convert_M(double* O2, double* N2, double* H2O, double* M);

// 2) 从 M[] 中拆分出 O2[], N2[], H2O[]
void convert_ONH2O(double* O2, double* N2, double* H2O, double* M);

// 3) 计算网格 (i,j) 的分子振动弛豫过程，将改变量存到 dM[] 中
//    - i, j 并未直接用到，但保留以满足您可能的定位需求
void vib_relaxation(
    int i, int j,
    double* M,         // 输入:  该网格处 O2,N2,H2O... “分子振动分布”
    double* dM,        // 输出:  每个振动态的改变量
    double Oz,         // 额外的 O原子数量/密度
    double** kvt,
    double** kvv0, double** kvv1, double** kvv2, double** kvv3, double** kvv4,
    double** re_kvt,
    double** re_kvv0, double** re_kvv1, double** re_kvv2, double** re_kvv3, double** re_kvv4,
    double T,          // 当前温度
    double** E         // 能量数据表
);

// 4) 读取文件, 填充 dE[][], E[][], 以及 c1[], c2[], c3[]
void Read_constant(double** dE, double** E, double* c1, double* c2, double* c3);

#endif // ENERGY_TYPE_H
