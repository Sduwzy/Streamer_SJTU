#ifndef SPLINE_H
#define SPLINE_H

// 如果要用到 `vec()` / `free_vec()` 函数，需要 include "memory.h"
#include "memory.h"  // 若您在 spline 中使用了 vec(), free_vec()

#ifdef __cplusplus
extern "C" {
#endif

// 全局常量 sp6(=1/6)
extern double sp6;

/**
 * @brief 分配 2 次微分 (y2[]) 用于三次样条插值
 *
 *  输入:
 *    x[]: 横坐标数组 (长度 n)
 *    y[]: 对应纵坐标 (长度 n)
 *    n  : 点的个数
 *  输出:
 *    y2[]: 计算后的二阶导数数组(长度 n)，供 splint() 使用
 */
void spline(double x[], double y[], int n, double y2[]);

/**
 * @brief 使用 spline() 得到的 y2[]，对指定 x 做三次样条插值
 *
 *  输入:
 *    xa[], ya[]: 节点 x,y 数组 (长度 n)
 *    y2a[]: spline() 算出的二阶导数 (长度 n)
 *    n : 节点个数
 *    x : 待插值点
 *  输出:
 *    *y: 在 x 处的插值结果
 */
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);

/**
 * @brief 与 splint 类似，但带 klo,khi 可重复插值，减少二分搜索开销
 *
 * 说明:
 *   - *klo, *khi 用来记录二分搜索区间，可在后续插值时复用
 *   - 其他逻辑与 splint 相同
 */
void splint_mod(double xa[], double ya[], double y2a[], int n, double x, double *y, int *klo, int *khi);

#ifdef __cplusplus
}
#endif

#endif // SPLINE_H
