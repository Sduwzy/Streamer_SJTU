#ifndef HYPERBOLIC_CURVE_H
#define HYPERBOLIC_CURVE_H

// 如果函数中需要数学函数，如 sqrt, pow, atan2, log, tan
#include <cmath>

// 函数声明：计算双曲型电极在 (r, z) 点的电压
// 参数说明：
//   Vele : 电极电压
//   cr   : 先端曲率半径
//   b    : ギャップ長
//   r, z : 要求电压的坐标
// 返回值 : (r, z)处的电压(若在电极内或接地处，返回电极电压/0)
double hyperbolic_curve(double Vele, double cr, double b, double r, double z);

#endif // HYPERBOLIC_CURVE_H
