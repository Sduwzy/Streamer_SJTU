#ifndef REACTION_H
#define REACTION_H

#include <cstdio>   // FILE*, fscanf, etc.
#include <cstdlib>  // exit
#include <cstring>  // strcmp, strcpy
#include <cctype>   // isspace

// ============================================================================
// 函数声明
// ============================================================================

/**
 * @brief   根据粒子名称字符串，在 particle[] 中查找对应的索引。
 * 
 * @param p             要查找的粒子名称 (C字符串)
 * @param particle      字符串数组，保存了所有粒子名称
 * @param n_particles   粒子总数
 * @return int          找到则返回对应索引，否则报错并退出
 */
int particle_number(char* p, char** particle, int n_particles);


/**
 * @brief   从文件中读取化学反应信息，并填充 reactl, reactr, A, B, E, ER 等。
 * 
 * @param fp            已打开的文件指针
 * @param NUM_L         每条反应的左边粒子数量
 * @param NUM_R         每条反应的右边粒子数量
 * @param n_particles   粒子总数
 * @param Rnum          用于存放读取到的“反应编号”数组
 * @param reactl        reactl[i][...] 存每条反应左侧粒子编号
 * @param reactr        reactr[i][...] 存每条反应右侧粒子编号
 * @param A, B, E, ER   存放反应相关的参数数组 (1维)
 * @param particle      粒子名称数组
 * @return int          读取到的反应条目数 n
 */
int Read_reaction(
    FILE* fp, 
    int  NUM_L, 
    int  NUM_R,
    int  n_particles,
    int* Rnum, 
    int** reactl,
    int** reactr,
    double* A,
    double* B,
    double* E,
    double* ER,
    char** particle
);


/**
 * @brief   从文件中读取粒子种类及其初始浓度，存入 particle[] 和 y[]。
 * 
 * @param fp        已打开的初始条件文件指针
 * @param particle  存放粒子名称的字符串数组
 * @param y         存放粒子初始浓度的双精度数组
 * @return int      读取到的粒子总数 (含“虚粒子” - 和 M)
 */
int Initial_condition(FILE* fp, char** particle, double* y);

#endif // REACTION_H
