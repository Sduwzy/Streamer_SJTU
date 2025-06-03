#include "first_q.h"

// 若需要进一步的头文件，也可引入
#include <cmath>   // std::exp, std::pow

// ----------------------------------------------------------------------------
// first_q(...)
// ----------------------------------------------------------------------------
void first_q(int NR,
             int NZ,
             int**    flag,
             double*  rh,
             double*  zh,
             double** ne,
             double** N2p,
             double** O2p,
             double** O2m,
             double** Om)
{
    // 设定一个最大密度 n_max
    const double n_max = 1.0e9; // [m^-3], 原注释

    // 这里默认一个高斯分布中心
    const double r0  = 0.0e-3;    // 中心 r=0.0 mm
    const double z0  = 13.0e-3;   // 中心 z=13.0 mm
    const double sr0 = 100.0e-6;  // 在 r 方向的标准差 100 um
    const double sz0 = 100.0e-6;  // 在 z 方向的标准差 100 um

    // -----------------------------
    // 1) 给 ne 赋初值
    // -----------------------------
    for(int i = 0; i < NR; ++i)
    {
        for(int j = 0; j < NZ; ++j)
        {
            // 如果该网格处是电极或其它特殊区域 (flag==1)，则 ne=0
            if(flag[i][j])
            {
                ne[i][j] = 0.0;
            }
            else
            {
                // 高斯分布
                double dr = (rh[i] - r0);
                double dz = (zh[j] - z0);

                double val =
                    n_max *
                    std::exp( - (std::pow(dr, 2) / (2.0 * std::pow(sr0, 2))
                               + std::pow(dz, 2) / (2.0 * std::pow(sz0, 2))) );

                ne[i][j] = val;
            }
        }
    }

    // -----------------------------
    // 2) 给 N2p, O2p, O2m, Om 赋初值
    // -----------------------------
    for(int i = 0; i < NR; ++i)
    {
        for(int j = 0; j < NZ; ++j)
        {
            if(flag[i][j])
            {
                // 如果是电极，则所有离子密度设为 0
                N2p[i][j] = 0.0;
                O2p[i][j] = 0.0;
                O2m[i][j] = 0.0;
                Om [i][j] = 0.0;
            }
            else
            {
                // 这里做简单分配：N2p = ne/3, O2p = ne/2 等...
                // 这是原始代码的一个假设
                double val_ne = ne[i][j];
                N2p[i][j] = val_ne / 3.0;
                O2p[i][j] = val_ne / 2.0;
                O2m[i][j] = 0.0; // 暂时设为 0
                Om [i][j] = 0.0; // 暂时设为 0
            }
        }
    }
}
