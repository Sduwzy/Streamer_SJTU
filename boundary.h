#ifndef BOUNDARY_H
#define BOUNDARY_H

// 如果需要使用 C++标准头文件，在这里包含
#include <cstddef>   // size_t
#include <cstdio>    // printf
#include <cmath>     // sqrt / fabs

// 函数声明
void negative_boundary_condition(int NR,int NZ, double **ne, int **flag);
void positive_boundary_condition(int num_x,int num_y, double **ion, int **flag);
void call_flag(int num_x,int num_y, int **flag, int **totuflag, int **otuflag, int **iflag_array, int **jflag_array);

void mol_boundary(
    int NR,int NZ,int **flag,
    double **ne,double **N2p,double **O2p, double **H2Op,
    double **O2m,double **Om,double **OHm,double **Hm,
    double **N4p,double **O4p,double **N2O2p,
    double **O2pH2O, double **H3Op, double **H3OpH2O,
    double **H3OpH2O2, double **H3OpH2O3,
    double **Ex, double **Ey, double **absE
);

#endif // BOUNDARY_H

