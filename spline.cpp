#include "spline.h"   // 包含头文件
#include <cstdio>     // 如果要用 printf 等
#include <cstdlib>    // 如果要用 exit
#include <cmath>      // 如果要用 sqrt 等

// 原代码的全局常量
double sp6 = 0.1666667; // 1.0 / 6.0

// spline 函数：分配 2 次微分 y2[] (三次样条)
void spline(double x[], double y[], int n, double y2[])
{
    int i;
    double p, sig, *u;

    // 分配一个临时数组 u
    u = vec(n-1);  // 来自 memory.h

    y2[0] = 0.0;
    u[0]  = 0.0;

    // 构建三次样条
    for (i = 1; i < n-1; i++) {
        sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
        p = sig * y2[i-1] + 2.0;
        y2[i] = (sig - 1.0) / p;

        u[i] = ( (y[i+1]-y[i])/(x[i+1]-x[i]) 
               - (y[i]   -y[i-1])/(x[i]   -x[i-1]) );
        u[i] = (6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1]) / p;
    }

    // 边界二阶导为 0
    y2[n-1] = 0.0;

    // 回代
    for (i = n-2; i >= 0; i--) {
        y2[i] = y2[i] * y2[i+1] + u[i];
    }

    free_vec(u, n-1);
}

// splint 函数：给定 (xa, ya, y2a), 长度 n，和目标点 x, 返回 *y
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
    int k, klo, khi;
    double inv_h, h, b, a;

    // 1) 二分搜索：找到 klo,khi 使得 xa[klo] <= x < xa[khi]
    klo = 0;
    khi = n - 1;
    while (khi - klo > 1) {
        k = (khi + klo) >> 1;
        if (xa[k] > x) {
            khi = k;
        } else {
            klo = k;
        }
    }

    // 2) 计算三次样条插值
    h = xa[khi] - xa[klo];
    if (h == 0.0) {
        std::fprintf(stderr, "Bad xa input to routine splint\n");
        return;
    }

    inv_h = 1.0 / h;

    a = (xa[khi] - x)*inv_h;
    b = (x - xa[klo])*inv_h;
    *y = a*ya[klo] + b*ya[khi] 
       + ((a*a*a - a)*y2a[klo] + (b*b*b - b)*y2a[khi]) * (h*h) * sp6;

    // 数据外插：若 x 超出 xa[n-1], 直接返回 ya[n-1]; 若 x < xa[0], 返回 ya[0]
    if( x > xa[n-1] ) *y = ya[n-1];
    if( x < xa[0] )   *y = ya[0];
}

// splint_mod: 与 splint 类似，但额外带 klo,khi 参数，减少重复二分搜索
void splint_mod(double xa[], double ya[], double y2a[], int n, double x, double *y, 
                int *klo, int *khi)
{
    int k;
    double inv_h, h, b, a;

    // 利用 *klo, *khi 做二分搜索
    while ((*khi) - (*klo) > 1) {
        k = ((*khi) + (*klo)) >> 1;
        if (xa[k] > x) {
            *khi = k;
        } else {
            *klo = k;
        }
    }

    // 三次样条插值
    h = xa[*khi] - xa[*klo];
    if (h == 0.0) {
        std::fprintf(stderr, "Bad xa input to routine splint_mod\n");
        return;
    }

    inv_h = 1.0 / h;
    a = (xa[*khi] - x)*inv_h;
    b = (x - xa[*klo])*inv_h;

    *y = a*ya[*klo] + b*ya[*khi] 
       + ((a*a*a - a)*y2a[*klo] + (b*b*b - b)*y2a[*khi])*(h*h)*sp6;

    // 数据外插
    if( x > xa[n-1] ) *y = ya[n-1];
    if( x < xa[0]   ) *y = ya[0];
}
