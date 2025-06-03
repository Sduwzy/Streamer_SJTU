#include "calc_velo.h"

// 若 spline.h 包含 spline(...) splint(...) 等函数声明，需包含：
#include "spline.h"  

// 如果使用到某些数学函数，比如 std::fabs, std::log:
#include <cmath>  
#include <cstdio> // 若需要 printf, ...
#include <cstdlib> // 若需要 exit, ...
// ... 视需要包含更多

void calc_e_velo(
    int num_x,
    int num_y,
    double** absE,
    double** Ex,
    double** Ey,
    double** velox,
    double** veloy,
    double** /*ne*/,  // 保留形参但无需使用时可注释形参名
    double*  /*rh*/,
    double*  /*zh*/,
    double* TdE,
    double* v_elec1,
    double* v_elec2,
    int boltNum
)
{
    // 分子数计算(原始代码放了一些注释/计算，这里简单保留):
    // const double Nm   = ((6.02e+23 * 1000.0) / (0.0820578 * 300)) * 1e-6;
    // const double NN   = ((6.02e+23 * 1000.0) / (0.0820578 * 300));
    // const double coeff= 3.74e22 / Nm;

    // 遍历网格
    for(int i = 0; i < num_x; ++i) {
        for(int j = 0; j < num_y; ++j) {
            // 1) 计算局部电场强度
            double localE = std::fabs(absE[i][j]);

            // 2) 在 v_elec1,v_elec2,Boltzmann数据表上插值,得到 v_elec
            double v_elec = 0.0;
            splint(TdE, v_elec1, v_elec2, boltNum, localE, &v_elec);

            // 3) 根据电场方向设置 x 方向速度
            if (Ex[i][j] < 0.0) {
                velox[i][j] =  v_elec * std::fabs(Ex[i][j]) * 1e-21;
            } else {
                velox[i][j] = -v_elec * std::fabs(Ex[i][j]) * 1e-21;
            }

            // 再次插值(或者可以复用上次 v_elec，如原始代码中那样再写一次)
            // 取决于你对 v_elec 是否需再次获取
            splint(TdE, v_elec1, v_elec2, boltNum, localE, &v_elec);

            // 4) y方向
            if (Ey[i][j] < 0.0) {
                veloy[i][j] =  v_elec * std::fabs(Ey[i][j]) * 1e-21;
            } else {
                veloy[i][j] = -v_elec * std::fabs(Ey[i][j]) * 1e-21;
            }
        }
    }
}


void calc_ie_velo(
    int num_x,
    int num_y,
    double** vx,
    double** vy,
    double** Ex,
    double** Ey,
    int mp
)
{
    // 这里是原逻辑: +离子速度可能用 keisu * E, -离子速度也可能类似
    // 但原代码中却给出了 vx=vy=0.0, 所以我们先保留注释:
    // 
    // 你可以根据需要恢复或修改:
    // const double MOL   = (6.02e+23*1000.0)/(0.0820578*300.0);
    // const double keisu = 2.2e-4*MOL*1e-21;

    if(mp == 0) {
        // +离子
        for(int i = 0; i < num_x; ++i) {
            for(int j = 0; j < num_y; ++j) {
                // vx[i][j] =  keisu * Ex[i][j];
                // vy[i][j] =  keisu * Ey[i][j];
                // 暂时按原代码保持 0
                vx[i][j] = 0.0;
                vy[i][j] = 0.0;
            }
        }
    }
    else if(mp == 1) {
        // -离子
        for(int i = 0; i < num_x; ++i) {
            for(int j = 0; j < num_y; ++j) {
                // vx[i][j] = -keisu * Ex[i][j];
                // vy[i][j] = -keisu * Ey[i][j];
                // 暂时按原代码保持 0
                vx[i][j] = 0.0;
                vy[i][j] = 0.0;
            }
        }
    }
    // 若 mp 取值不止 0 或 1，这里未做检查
}
