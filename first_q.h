#ifndef FIRST_Q_H
#define FIRST_Q_H

// 必要的头文件
#include <cmath>  // 如果需要 std::exp / std::pow 等
#include <cstdio> // 如果需要 printf（可选）
#include <cstdlib>// 如果可能用到 exit（可选）

/**
 * @brief 在网格 (NR x NZ) 上分配一些初始电荷密度 (如 ne, N2p, O2p, O2m, Om 等)。
 *
 * @param NR   网格在径向方向的大小
 * @param NZ   网格在轴向方向的大小
 * @param flag 2D 整型数组 (NR x NZ), 标记哪些网格处是电极(或需要特殊处理)
 * @param rh   径向网格坐标 (长度NR)
 * @param zh   轴向网格坐标 (长度NZ)
 * @param ne   输出的电子密度数组 (NR x NZ)
 * @param N2p  输出的 N2+ 密度数组 (NR x NZ)
 * @param O2p  输出的 O2+ 密度数组 (NR x NZ)
 * @param O2m  输出的 O2- 密度数组 (NR x NZ)
 * @param Om   输出的 O-  密度数组 (NR x NZ)
 */
void first_q(int NR, int NZ,
             int**    flag,
             double*  rh,
             double*  zh,
             double** ne,
             double** N2p,
             double** O2p,
             double** O2m,
             double** Om);

#endif // FIRST_Q_H
