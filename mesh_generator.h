#ifndef MESH_GENERATOR_H
#define MESH_GENERATOR_H

#include <cstdio>  // FILE*, fopen, fscanf, fclose
#include <cstdlib> // malloc/free (若需要)
#include <cmath>   // 如果有用到 sqrt pow

// mesh_generator 函数声明：
// 根据给定的网格文件(FILE_MESH_R, FILE_MESH_Z)，读取网格节点坐标到r[], z[]，
// 计算中心坐标rh[], zh[]，以及 dr[], dz[], 体积Vol[][], 表面积Sr[][], Sz[][] 等。
//
// 参数含义：
//   FILE_MESH_R, FILE_MESH_Z: 存储 r,z 网格节点数据的文件路径
//   NR, NZ: 网格在 r, z 方向的尺寸
//   r, z:     一维数组，存储网格节点坐标
//   rh, zh:   一维数组，存储网格单元中心(或半点)坐标
//   dr, dz:   一维数组，存储网格宽度
//   Vol, Sr, Sz: 二维数组, 分别表示网格单元体积/表面积(侧面/截面?)
// 
// 函数无返回值，计算结果填入传入的数组中
//
void mesh_generator(const char *FILE_MESH_R, const char *FILE_MESH_Z,
                    int NR, int NZ,
                    double *r, double *z,
                    double *rh, double *zh,
                    double *dr, double *dz,
                    double **Vol, double **Sr, double **Sz);

#endif // MESH_GENERATOR_H
