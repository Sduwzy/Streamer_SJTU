#include "poisson.h"   // 包含对应的头文件
#include <cstdio>      // for std::printf
#include <cmath>       // for std::sqrt, std::fabs
#include <iostream>    // 若需要 std::cout

// 如果您还需要其他头文件，也可引入
// #include "..."


void discretization(int N,int M,
                    double *rhalf,double *zhalf,
                    double **P1,double **P2,double **P3,double **P4,double **P5,
                    double a,double b,
                    int **iflag,int **jflag,int **otuflag,int **flag)
{
    // 1) 分配局部矩阵 h1,h2,h3,h4
    double** h1 = mat(N, M);
    double** h2 = mat(N, M);
    double** h3 = mat(N, M);
    double** h4 = mat(N, M);

    // 2) 先给 h1..h4 赋值
    for(int i=0; i<N; i++){
        for(int j=0; j<M; j++){
            if(i==0) {
                if(j==M-1) {
                    h1[i][j] = h1[i][j-1];
                } else {
                    h1[i][j] = zhalf[j+1] - zhalf[j];
                }
                h2[i][j] = rhalf[i];
                h3[i][j] = rhalf[i+1] - rhalf[i];

                if(j==0) {
                    h4[i][j] = zhalf[j];
                } else {
                    h4[i][j] = zhalf[j] - zhalf[j-1];
                }
            }
            else if(j==0) {
                h1[i][j] = zhalf[j+1] - zhalf[j];
                h2[i][j] = rhalf[i] - rhalf[i-1];

                if(i==N-1) {
                    h3[i][j] = h3[i-1][j];
                } else {
                    h3[i][j] = rhalf[i+1] - rhalf[i];
                }
                h4[i][j] = zhalf[j];
            }
            else if(j==M-1) {
                h1[i][j] = h1[i][j-1];
                h2[i][j] = rhalf[i] - rhalf[i-1];

                if(i==N-1) {
                    h3[i][j] = h3[i-1][j];
                } else {
                    h3[i][j] = rhalf[i+1] - rhalf[i];
                }
                h4[i][j] = zhalf[j] - zhalf[j-1];
            }
            else if(i==N-1) {
                if(j==M-1) {
                    h1[i][j] = h1[i][j-1];
                } else {
                    h1[i][j] = zhalf[j+1] - zhalf[j];
                }
                h2[i][j] = rhalf[i] - rhalf[i-1];
                h3[i][j] = h3[i-1][j];
                h4[i][j] = zhalf[j] - zhalf[j-1];
            }
            else {
                h1[i][j] = zhalf[j+1] - zhalf[j];
                h2[i][j] = rhalf[i] - rhalf[i-1];
                h3[i][j] = rhalf[i+1] - rhalf[i];
                h4[i][j] = zhalf[j] - zhalf[j-1];
            }
        } // end for j
    } // end for i

    // 3) 根据 h1..h4 构造 P1..P5
    double hh1, hh2, hh3, hh4;
    double r0_double;

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < M; j++) {
            if(flag[i][j]) {
                // 若此网格在电极内(或其他标记), 则令系数=0
                P1[i][j] = 0.0; P2[i][j] = 0.0;
                P3[i][j] = 0.0; P4[i][j] = 0.0;
            } else {
                // 判断 otuflag, jflag => hh1
                if(otuflag[i][j] == 1 || jflag[i][j] == 1) {
                    hh1 = b*std::sqrt(1.0 + std::pow(rhalf[i]/a,2)) - zhalf[j];
                } else {
                    hh1 = h1[i][j];
                }

                // 判断 iflag, otuflag => hh2
                if(iflag[i][j] == 1 || otuflag[i][j] == 1) {
                    hh2 = rhalf[i] - a*std::sqrt(std::pow(zhalf[j]/b,2) - 1.0);
                } else {
                    hh2 = h2[i][j];
                }

                hh3 = h3[i][j];

                if(j==0) {
                    hh4 = zhalf[j];
                } else {
                    hh4 = h4[i][j];
                }

                r0_double = 2.0 * rhalf[i];

                // 根据 i 是否为 0 做不同处理
                if(i==0) {
                    P1[i][j] = (hh1 * hh4 * hh3 * hh3) / ( hh3*hh3 + 2.0*hh1*hh4 );
                    P2[i][j] = -1.0 / (hh1*(hh1+hh4));
                    P3[i][j] = 0.0;
                    P4[i][j] = -2.0 / (hh3*hh3);
                    P5[i][j] = -1.0 / (hh4*(hh1+hh4));
                } else {
                    P1[i][j] = hh1*hh2*hh3*hh4 / (r0_double*hh2*hh3 + (r0_double + hh2 - hh3)*hh1*hh4);
                    P2[i][j] = -r0_double / (hh1*(hh1+hh4));
                    P3[i][j] = (-r0_double + hh3) / (hh2*(hh2+hh3));
                    P4[i][j] = -(r0_double + hh2)/(hh3*(hh2+hh3));
                    P5[i][j] = -r0_double/(hh4*(hh1+hh4));
                }
            }
        }
    }

    // 4) 释放局部矩阵
    free_mat(h1, N, M);
    free_mat(h2, N, M);
    free_mat(h3, N, M);
    free_mat(h4, N, M);
}



void poiseq(int num_x,int num_y,
            double **pphi, int **flag,
            int **iflag, int **jflag, int **otuflag,
            double *rhalf, double *zhalf,
            double **rho,
            double **P1,double **P2,double **P3,double **P4,double **P5)
{
    // 收敛阈值(可改小/大)
    double Conv  = 1.0e-5;    
    double OMEGA = 1.9;  // SOR松弛系数

    int loop = 0;
    double MaxPhi = 1.0e-20;
    double MaxErr = 0.0;
    double Prev_phi = 0.0;
    double CurErr = 0.0;

    // 临时数组
    double** gsphi = mat(num_x, num_y);

    do {
        MaxErr = 0.0;
        // 每个循环都要更新全场
        for(int i=0; i<num_x; i++){
            for(int j=0; j<num_y; j++){
                Prev_phi = pphi[i][j]; // 保存上一迭代值

                if(flag[i][j]) {
                    gsphi[i][j] = 0.0;
                } else {
                    // 这里分情况
                    if(i == 0) {
                        gsphi[i][j] = P1[i][j]*(
                            rho[i][j]
                            - P2[i][j]*pphi[i][j+1]
                            - P4[i][j]*pphi[i+1][j]
                            - P5[i][j]*pphi[i][j-1]
                        );
                        // 若 j==0, 需要减 pphi[i][j-1], 但 j-1=-1 => 说明边界为 0?
                        if(j==0){
                            gsphi[i][j] = P1[i][j]*(
                                rho[i][j]
                                - P2[i][j]*pphi[i][j+1]
                                - P4[i][j]*pphi[i+1][j]
                                // pphi[i][j-1] => 0
                            );
                        }
                    }
                    else if(j==0) {
                        gsphi[i][j] = P1[i][j]*(
                            rho[i][j]
                            - P2[i][j]*pphi[i][j+1]
                            - P3[i][j]*pphi[i-1][j]
                            - P4[i][j]*pphi[i+1][j]
                            // pphi[i][j-1] => 0
                        );
                    }
                    else if(i == num_x-1) {
                        // 右边界(或 i==num_x-1)做什么？
                        gsphi[i][j] = gsphi[i-1][j];
                    }
                    else if(j == num_y - 1) {
                        // 顶端?
                        gsphi[i][j] = gsphi[i][j-1];
                    }
                    else {
                        // 常规内点
                        gsphi[i][j] = P1[i][j]*(
                            rho[i][j]
                            - P2[i][j]*pphi[i][j+1]
                            - P3[i][j]*pphi[i-1][j]
                            - P4[i][j]*pphi[i+1][j]
                            - P5[i][j]*pphi[i][j-1]
                        );
                    }
                }

                // SOR 更新
                pphi[i][j] = pphi[i][j] + OMEGA*( gsphi[i][j] - pphi[i][j] );

                // 更新 MaxPhi
                if(std::fabs(pphi[i][j]) > MaxPhi) {
                    MaxPhi = std::fabs(pphi[i][j]);
                }

                // 计算更新差
                CurErr = std::fabs(pphi[i][j] - Prev_phi);
                if(CurErr > MaxErr) {
                    MaxErr = CurErr;
                }
            }
        }

        loop++;
        if(loop % 100 == 0) {
            std::printf("%d\t%e\t%e\n", loop, MaxErr/MaxPhi, MaxPhi);
        }
    } while( (MaxErr/MaxPhi) > Conv );

    std::printf("---Finish_poisson_equation---");
    std::printf("%05d  %e\n", loop, MaxPhi);

    free_mat(gsphi, num_x, num_y);
}
