#ifndef CALC_E_H
#define CALC_E_H

// 如果函数内部需要用到一些标准函数声明，可以在这里 include <cmath>, <cstdlib> 等；
// 不过头文件只放声明的话，可以只写以下声明；具体实现在 calcE.cpp

// 声明一个分配/释放二维数组的函数(如需要), 也可在别处统一声明
//double** allocate2DDouble(int rows, int cols);
//void free2DDouble(double** arr, int rows);

// 声明 calc_E(...) 函数
void calc_E(int NR,int NZ,
            double **phi,   // 电位
            double **absE,  // 电场大小
            double **Ey,    // y方向电场
            double **Ex,    // x方向电场
            double air_kg,
            int **totuflag,
            int **otuflag,
            int **iflag,
            int **jflag,
            int **flag,
            double a,
            double b,
            double *rhalf,
            double *zhalf,
            double **Mol,
            int num,
            int Prinstp);

#endif // CALC_E_H

