#ifndef SUPERBEE_HPP_
#define SUPERBEE_HPP_

#include <cmath>
#include <algorithm> // std::max/min

//---------------------------------------------------------------------------
// 内联工具函数
//---------------------------------------------------------------------------

// 三者取最大
inline double max3(double a, double b, double c)
{
    return std::max(a, std::max(b, c));
}

// 两者取最小
inline double min2(double a, double b)
{
    return std::min(a, b);
}

/**
 * super-bee（或 minmod）限制器
 * @param r1   ∇φ_{i  }   左差分
 * @param r2   ∇φ_{i+1}   右差分
 * @param b    b = 2 ⇒ superbee；b = 1 ⇒ minmod
 */
inline double minmod(double r1, double r2, double b)
{
    const double sgn = (r2 < 0.0) ? -1.0 : 1.0;
    return sgn * max3(0.0,
                      min2(sgn * b * r1, std::fabs(r2)),
                      min2(sgn * r1, b * std::fabs(r2)));
}

//---------------------------------------------------------------------------
// MUSCL-SuperBee 核心 3 组函数
// （声明；定义在 superbee.cpp 中）
//---------------------------------------------------------------------------
void MUSCL_superbee_methoed_for_pion(double **rho, double **u, double **v,
                                     double dt, double kappa, double b,
                                     int nx, int ny,
                                     double **Sr, double **Sz, double **Vol,
                                     int **iflag, int **jflag, int **otuflag);

void MUSCL_superbee_methoed_for_e(double **rho, double **u, double **v,
                                  double dt, double kappa, double b,
                                  int nx, int ny,
                                  double **Sr, double **Sz, double **Vol,
                                  int **iflag, int **jflag, int **otuflag);

void MUSCL_superbee_methoed_for_mion(double **rho, double **u, double **v,
                                     double dt, double kappa, double b,
                                     int nx, int ny,
                                     double **Sr, double **Sz, double **Vol,
                                     int **iflag, int **jflag, int **otuflag);

#endif // SUPERBEE_HPP_