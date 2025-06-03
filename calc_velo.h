
#ifndef CALC_VELOCITY_H
#define CALC_VELOCITY_H

// 如果后续要用到 std::size_t、std::vector 等，可以包含其他头文件
#include <cstddef>

/// \brief 计算电子速度场
///
/// \param num_x   网格在 x 方向上的大小
/// \param num_y   网格在 y 方向上的大小
/// \param absE    电场绝对值数组，大小 [num_x][num_y]
/// \param Ex      x 方向电场分量，大小 [num_x][num_y]
/// \param Ey      y 方向电场分量，大小 [num_x][num_y]
/// \param velox   输出：电子在 x 方向上的速度分量
/// \param veloy   输出：电子在 y 方向上的速度分量
/// \param ne      电子密度场(本函数中未实际用到，但保留了形参)
/// \param rh      网格点 x 坐标(或半径)
/// \param zh      网格点 y 坐标(或高度)
/// \param TdE     BOLSIG/Boltzmann 数据中的 E 轴
/// \param v_elec1 BOLSIG/Boltzmann 数据, 与 TdE 对应
/// \param v_elec2 BOLSIG/Boltzmann 二次微分数据(若使用 spline，需要二阶导)
/// \param boltNum BOLSIG 数据数组大小
void calc_e_velo(
    int num_x,
    int num_y,
    double** absE,
    double** Ex,
    double** Ey,
    double** velox,
    double** veloy,
    double** ne,
    double* rh,
    double* zh,
    double* TdE,
    double* v_elec1,
    double* v_elec2,
    int boltNum
);


/// \brief 计算正/负离子的速度场
///
/// \param num_x 网格在 x 方向上的大小
/// \param num_y 网格在 y 方向上的大小
/// \param vx    输出：离子在 x 方向上的速度分量
/// \param vy    输出：离子在 y 方向上的速度分量
/// \param Ex    x 方向电场分量
/// \param Ey    y 方向电场分量
/// \param mp    0 表示 +离子，1 表示 -离子
void calc_ie_velo(
    int num_x,
    int num_y,
    double** vx,
    double** vy,
    double** Ex,
    double** Ey,
    int mp
);

#endif // CALC_VELOCITY_H
