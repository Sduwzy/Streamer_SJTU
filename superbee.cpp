#include "superbee.h"

// 项目里已有的矩阵工具函数
extern double **mat(int nx, int ny);
extern void free_mat(double **m, int nx, int ny);

//------------------------------------------------------------------------------
// 下方三个函数与您 C 版代码逐字等价，仅把内部 min/max 宏替换为
//  inline max3/min2 以避免宏污染；其余保留双指针 + 手写循环结构。
//------------------------------------------------------------------------------

//-------------------------------------------------------------------------
// helper-template 宏，用于复用实现（不改变原有变量名，方便对照）
//-------------------------------------------------------------------------
#define IMPLEMENT_MUSCL_BODY(NameSuffix, FluxSign)                                                                                                                                                    \
    void NameSuffix(double **ro, double **u, double **v,                                                                                                                                              \
                    double dt, double kappa, double b, int numx, int numy,                                                                                                                            \
                    double **Sr, double **Sz, double **Vol,                                                                                                                                           \
                    int **iflag, int **jflag, int **otuflag)                                                                                                                                          \
    {                                                                                                                                                                                                 \
        int i, j;                                                                                                                                                                                     \
        double **UL = mat(numx, numy);                                                                                                                                                                \
        double **UR = mat(numx, numy);                                                                                                                                                                \
        double **Fl = mat(numx, numy);                                                                                                                                                                \
        double **Gl = mat(numx, numy);                                                                                                                                                                \
                                                                                                                                                                                                      \
        /* ---------------- 计算 Gl 半通量 ---------------- */                                                                                                                                        \
        for (i = 0; i < numx; ++i)                                                                                                                                                                    \
        {                                                                                                                                                                                             \
            for (j = 1; j < numy; ++j)                                                                                                                                                                \
            {                                                                                                                                                                                         \
                UL[i][j] = ro[i][j] + 0.25 * ((1.0 - kappa) * minmod(ro[i][j + 1] - ro[i][j],                                                                                                         \
                                                                     ro[i][j] - ro[i][j - 1], b) +                                                                                                    \
                                              (1.0 + kappa) * minmod(ro[i][j] - ro[i][j - 1],                                                                                                         \
                                                                     ro[i][j + 1] - ro[i][j], b));                                                                                                    \
                UR[i][j] = ro[i][j] - 0.25 * ((1.0 - kappa) * minmod(ro[i][j] - ro[i][j - 1],                                                                                                         \
                                                                     ro[i][j + 1] - ro[i][j], b) +                                                                                                    \
                                              (1.0 + kappa) * minmod(ro[i][j + 1] - ro[i][j],                                                                                                         \
                                                                     ro[i][j] - ro[i][j - 1], b));                                                                                                    \
            }                                                                                                                                                                                         \
        }                                                                                                                                                                                             \
        j = 0;                                                                                                                                                                                        \
        for (i = 0; i < numx; ++i)                                                                                                                                                                    \
        {                                                                                                                                                                                             \
            UL[i][j] = ro[i][j] + 0.25 * ((1 - kappa) * minmod(ro[i][j + 1] - ro[i][j], ro[i][j] - 0.0, b) + (1 + kappa) * minmod(ro[i][j] - 0.0, ro[i][j + 1] - ro[i][j], b));                       \
            UR[i][j] = ro[i][j] - 0.25 * ((1 - kappa) * minmod(ro[i][j] - 0.0, ro[i][j + 1] - ro[i][j], b) + (1 + kappa) * minmod(ro[i][j + 1] - ro[i][j], ro[i][j] - 0.0, b));                       \
        }                                                                                                                                                                                             \
        for (i = 0; i < numx; ++i)                                                                                                                                                                    \
        {                                                                                                                                                                                             \
            for (j = 1; j < numy; ++j)                                                                                                                                                                \
            {                                                                                                                                                                                         \
                if (v[i][j] FluxSign 0.0)                                                                                                                                                             \
                    Gl[i][j] = dt * Sz[i][j] * (v[i][j] * UL[i][j - 1]);                                                                                                                              \
                else                                                                                                                                                                                  \
                    Gl[i][j] = dt * Sz[i][j] * (v[i][j] * UR[i][j]);                                                                                                                                  \
            }                                                                                                                                                                                         \
        }                                                                                                                                                                                             \
        for (i = 0; i < numx; ++i)                                                                                                                                                                    \
            Gl[i][0] = dt * Sz[i][0] * (v[i][0] * UR[i][0]);                                                                                                                                          \
                                                                                                                                                                                                      \
        /* ---------------- 计算 Fl 半通量 ---------------- */                                                                                                                                        \
        for (i = 1; i < numx; ++i)                                                                                                                                                                    \
        {                                                                                                                                                                                             \
            for (j = 0; j < numy; ++j)                                                                                                                                                                \
            {                                                                                                                                                                                         \
                UL[i][j] = ro[i][j] + 0.25 * ((1 - kappa) * minmod(ro[i + 1][j] - ro[i][j], ro[i][j] - ro[i - 1][j], b) + (1 + kappa) * minmod(ro[i][j] - ro[i - 1][j], ro[i + 1][j] - ro[i][j], b)); \
                UR[i][j] = ro[i][j] - 0.25 * ((1 - kappa) * minmod(ro[i][j] - ro[i - 1][j], ro[i + 1][j] - ro[i][j], b) + (1 + kappa) * minmod(ro[i + 1][j] - ro[i][j], ro[i][j] - ro[i - 1][j], b)); \
            }                                                                                                                                                                                         \
        }                                                                                                                                                                                             \
        i = 0;                                                                                                                                                                                        \
        for (j = 0; j < numy; ++j)                                                                                                                                                                    \
        {                                                                                                                                                                                             \
            UL[i][j] = ro[i][j] + 0.25 * ((1 - kappa) * minmod(ro[i + 1][j] - ro[i][j], ro[i][j] - 0.0, b) + (1 + kappa) * minmod(ro[i][j] - 0.0, ro[i + 1][j] - ro[i][j], b));                       \
            UR[i][j] = ro[i][j] - 0.25 * ((1 - kappa) * minmod(ro[i][j] - 0.0, ro[i + 1][j] - ro[i][j], b) + (1 + kappa) * minmod(ro[i + 1][j] - ro[i][j], ro[i][j] - 0.0, b));                       \
        }                                                                                                                                                                                             \
        for (i = 1; i < numx; ++i)                                                                                                                                                                    \
        {                                                                                                                                                                                             \
            for (j = 0; j < numy; ++j)                                                                                                                                                                \
            {                                                                                                                                                                                         \
                if (u[i][j] > 0.0)                                                                                                                                                                    \
                    Fl[i][j] = dt * Sr[i][j] * (u[i][j] * UL[i - 1][j]);                                                                                                                              \
                else                                                                                                                                                                                  \
                    Fl[i][j] = dt * Sr[i][j] * (u[i][j] * UR[i][j]);                                                                                                                                  \
            }                                                                                                                                                                                         \
        }                                                                                                                                                                                             \
        for (j = 0; j < numy; ++j)                                                                                                                                                                    \
            Fl[0][j] = 0.0;                                                                                                                                                                           \
                                                                                                                                                                                                      \
        /* ---------------- 边界/障碍物修正并更新密度 ---------------- */                                                                                                                             \
        for (i = 0; i < numx - 1; ++i)                                                                                                                                                                \
        {                                                                                                                                                                                             \
            for (j = 0; j < numy - 1; ++j)                                                                                                                                                            \
            {                                                                                                                                                                                         \
                if (iflag[i][j])                                                                                                                                                                      \
                {                                                                                                                                                                                     \
                    Fl[i][j] = dt * Sr[i][j] * u[i][j] * ro[i][j];                                                                                                                                    \
                }                                                                                                                                                                                     \
                if (jflag[i][j])                                                                                                                                                                      \
                {                                                                                                                                                                                     \
                    Gl[i][j + 1] = dt * Sz[i][j + 1] * v[i][j + 1] * ro[i][j];                                                                                                                        \
                }                                                                                                                                                                                     \
                if (otuflag[i][j])                                                                                                                                                                    \
                {                                                                                                                                                                                     \
                    Fl[i][j] = dt * Sr[i][j] * u[i][j] * ro[i][j];                                                                                                                                    \
                    Gl[i][j + 1] = dt * Sz[i][j + 1] * v[i][j + 1] * ro[i][j];                                                                                                                        \
                }                                                                                                                                                                                     \
                ro[i][j] -= (Fl[i + 1][j] - Fl[i][j] + Gl[i][j + 1] - Gl[i][j]) / Vol[i][j];                                                                                                          \
            }                                                                                                                                                                                         \
        }                                                                                                                                                                                             \
                                                                                                                                                                                                      \
        free_mat(UL, numx, numy);                                                                                                                                                                     \
        free_mat(UR, numx, numy);                                                                                                                                                                     \
        free_mat(Fl, numx, numy);                                                                                                                                                                     \
        free_mat(Gl, numx, numy);                                                                                                                                                                     \
    }                                                                                                                                                                                                 \
//--- 宏 END ------------------------------------------------------------------

// 生成三个实例
IMPLEMENT_MUSCL_BODY(MUSCL_superbee_methoed_for_pion, >)
IMPLEMENT_MUSCL_BODY(MUSCL_superbee_methoed_for_e, >=)
IMPLEMENT_MUSCL_BODY(MUSCL_superbee_methoed_for_mion, >)

#undef IMPLEMENT_MUSCL_BODY
