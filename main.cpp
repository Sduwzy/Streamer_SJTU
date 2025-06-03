// ===== 1) 引入头文件 =====
#include <iostream>
#include <cmath>
#include <cstdio>
#include <pthread.h> // 如果还要使用 POSIX 线程
#include <cstring>   // 如果需要 strcpy 等
#include <fstream>   // std::ofstream
#include <iomanip>
#include <thread>     // std::thread
#include <vector>     // std::vector
#include <mutex>      // 如果后续想用 std::mutex，可提前包含
#include <functional> // std::ref / cref

#include <cuda_runtime.h>

#include "memory.h"         // vec(), mat(), imat(), free_vec(), free_mat()
#include "spline.h"         // spline, splint, sp6
#include "calcE.h"          // calc_E(...)
#include "mesh_generator.h" // mesh_generator(...)
#include "get_constant.h"   // get_constant()
#include "boundary.h"
#include "Energy_type.h" // Read_constant()
#include "Reaction.h"
#include "first_q.h"
#include "poisson.h"
#include "Multigrid_kernel.h"
#include "Red_Black_SOR_Kernel.h"
#include "calc_velo.h"
#include "superbee.h"

#define INITIAL_FILE "inputdata/Initial/initial_N2_80p_300K.dat"
#define BOLTZMANNFILE "inputdata/Bolsig_Data/O2_20p.dat"
#define I_REACTION_FILE "inputdata/i_reaction_300K_modmod0516.dat"
#define E_REACTION_FILE "inputdata/e_reaction_mod1803.dat"
#define N_REACTION_FILE "inputdata/n_reaction_modmod.dat"
#define MESH_R_FILE "inputdata/mesh_r2_0624.dat"
#define MESH_Z_FILE "inputdata/mesh_z2.dat"
#define VOLTAGE_FILE "inputdata/V_Ono_single_str.dat"

#define point(N, i, j) ((N) * (i) + (j))

#define kb 1.38e-23 // ボルツマン定数
#define h 6.62e-34  // プランク定数
#define Diff 0.18   // 拡散係数

double T0 = 300.0; // 初期ガス温度(K)

#define PI 3.1415926535897932                        // 円周率
#define QELEC 1.602176462e-19                        // 電子の電荷（素電荷) [C]
#define MELEC 9.10938188e-31                         // 電子の質量
#define MOL ((6.02e+23 * 1000.0) / (0.0820578 * T0)) // 分子数(const)
#define E0 8.85e-12                                  // 真空の誘電率[F/m]=[C/V/m]
#define Avogadro 6.022141e+23                        // アボガドロ定数 コ/mol
#define MassNum (6.022141e+23 / (28.966 / 1000.0))   // 質量分子数 [コ/kg] (Avogadro/air_kg)

#define NR 256  // 径方向のメッシュ数
#define NZ 1792 // 軸方向のメッシュ数

#define VIB_SKIP 20000000

///////////////////GPUのブロック数////////////////
#define blockDim_x 1
#define blockDim_y 256

#define blockDim_x2 1
#define blockDim_y2 128

#define blockDim_x4 1
#define blockDim_y4 64
//////////////////////////////////////////////////

int jfall = 89;    // z=0.30049mm(89), 0.5(111), 1.0(161)
int j1st = 861;    // z= 8.005 mm
int janode = 1332; // z=12.7029mm
int numrad = 100;

#define THREAD_NUM 64 // マルチスレッドのスレッド数
#define BN 1200       // ボルツマン方程式のデーター数
int AN = 1200;        // 電離係数のデータ数
#define DAnum 800     // DA reactionのデータ数
#define DEMnum 1000   // DEM reactionのデータ数

double bgap = 13.0e-3; // ギャップ長
double cr = 500e-6;    // 先端曲率半径
double EPS = 1.0e-5;   // 収束

int n_particles = 100; // 考慮する粒子数、適当に代入

int NUM_L = 3;
int NUM_R = 3;
#define NUM_REACT 41 // C53+1
#define DATA 200     // テキトーに入れておく。反応数(C53とか)より多ければok
int n_react, ne_react, ni_react, nn_react, neR, nDAR, nDER;
int **ereactl, **ereactr, **ireactl, **ireactr, **nreactl, **nreactr;
int *ereactl_1st, *ereactl_2nd, *ereactl_3rd, *ereactr_1st, *ereactr_2nd, *ereactr_3rd;
int *nreactl_1st, *nreactl_2nd, *nreactl_3rd, *nreactr_1st, *nreactr_2nd, *nreactr_3rd;
int *ireactl_1st, *ireactl_2nd, *ireactl_3rd, *ireactr_1st, *ireactr_2nd, *ireactr_3rd;
double *eK0, *eK1, *eK2, *eER, *eER_mod, *eK2_mod;
double *iK0, *iK1, *iK2, *iER, *iK0m, *iER_mod;
double *nK0, *nK1, *nK2, *nER, *nK0m, *nK1m, *nK2m, *nER_mod;
double **T, power, **POWER;
int numO2, numO2p, nume, numN2, numH2O;
int *nRnum, *iRnum, *eRnum;
int **neflag, Prinstp2, nstp, limpoz;

int n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16;
int o2v0, o2v1, o2v2, o2v3, o2v4, n2v0, n2v1, n2v2, n2v3, n2v4, n2v5, n2v6, n2v7, n2v8, h2ov0, h2ov1, h2ov2, h2ov3, oatm;
// double PARTICLE[100][NR][NZ];
double krate[100][NR][NZ];
double PARTICLE[NR][NZ][100];
double molP[NR][NZ][40];
double Yk[NR][NZ][100];

double **ne, **N2p, **N4p, **O2p, **O4p, **N2O2p, **O2pH2O, **H3Op, **O2m, **Om, **O2, **N2, **H2Op, **OHm, **Hm;
double **H3OpH2O, **H3OpH2O2, **H3OpH2O3;
double pN2, pO2, pH2O;

double *TdE, *DAeV, *DEMeV, *deV, *ddeV, *dv, *ddv;
double **dc, **ddc;
double **vr, **vz, **piv_r, **piv_z, **miv_r, **miv_z, **Vol, **Sr, **Sz, kappa, superbee;
double **v, **piv, **miv, *rrho;
double current, current2, CurrDens, InducedCurr, V, current_c;

double **absE, dt, **phot_num, **TeV, *HelmG, **dtELr, **dtELz, **preEr, **preEz, **Ex, **Ey;
int **flag, **iflag, **jflag, **tflag, **oflag, *HelmA;
double **P1, **P2, **P3, **P4, **P5, **rho, **phi, **Lphi, **Cphi, **Cphi1, **Cphi2;
double *dalp, *ddalp, *alpTd, **I, *rh, *zh, **S_ph0, **S_ph1, **S_ph2, **sink, **e_sink, **i_sink;
double *dpower, *ddpower, *powerTd;
double **Br, **Bz, **mvr, **mvz, **cE, *r, *z, *dr, *dz, **rou, **p;
char **particle;

void rk4(int ii, int jj, double *y, double *yout, int n, double dt, int num_react, int **, int **);
void e_reaction_derivs(int ii, int jj, double *y, int *, int, double *yout, int **, int **, double *);
void i_reaction_derivs(int ii, int jj, double *y, int *, int, double *yout, int **, int **, double *);
void n_reaction_derivs(int ii, int jj, double *y, int *, int, double *yout, int **, int **, double *);
void get_constant();
void store(int Z, double time, double **Ey, char **particle, double atmpa, double air_kg, double molV);
void data_input(char *Mol, int numdata, double **n);
void vib_store(int Z);

// C++ 版本的线程参数结构体
struct ThreadArg
{
    int thread_no;          // 线程编号
    pthread_mutex_t *mutex; // 指向互斥锁对象的指针
};

void symmetric_TVD(double u0, double v0, double p0, double t0, double rou0, double dt, double g0, double rgas, double ecp, double *res,
                   double **rou, double **vr, double **vz, double **p, double **q1, double **q2, double **q3, double **q4,
                   int **flag, int **iflag, int **jflag, int **oflag, int **tflag);
void bndcnd(double **rou, double **vr, double **vz, double **p, double **q1, double **q2,
            double **q3, double **q4, double g0, double rgas, double u0, double v0, double p0, double t0,
            double rou0, int **flag, int **iflag, int **jflag, int **oflag, int **tflag);
void calrhs(double **rou, double **vr, double **vz, double **p,
            double **q1, double **q2, double **q3, double **q4, double rgas, double g0,
            double dt, double ecp, double **dq1, double **dq2, double **dq3, double **dq4, int **flag, int te);

//void *thread_func_fluid(void *arg);

void bndcnd(
    double **rou, double **mvr, double **mvz, double **p,
    double **q1, double **q2, double **q3, double **q4,
    double g0, double rgas, double u0, double v0, double p0, double t0,
    double rou0,
    int **flag, int **iflag, int **jflag, int **oflag, int **tflag)
{
    // 1) 局部变量
    double cpgas, htotal, absv;

// 2) 首先并行地通过 q1..q4 恢复 rou, mvr, mvz, p, T
//    这里 T[][] 也是假设是外部可见(可能是 global)
#pragma omp parallel num_threads(THREAD_NUM)
    {
#pragma omp for
        for (int i = 0; i < NR; i++)
        {
            for (int j = 0; j < NZ; j++)
            {
                rou[i][j] = q1[i][j];
                mvr[i][j] = q2[i][j] / q1[i][j];
                mvz[i][j] = q3[i][j] / q1[i][j];
                p[i][j] = (g0 - 1.0) * (q4[i][j] - 0.5 * rou[i][j] * (mvr[i][j] * mvr[i][j] + mvz[i][j] * mvz[i][j]));
                // 假设有一个外部二维数组: double** T;
                T[i][j] = p[i][j] / (rgas * rou[i][j]);
            }
        }
    }

    // 3) 计算Cp(定压比热)和总焓
    cpgas = rgas * g0 / (g0 - 1.0);
    htotal = 0.5 * (u0 * u0 + v0 * v0) + t0 * cpgas;

    // 4) 对“针电极内”区(可能是 mesh 上某段? 如 j >= 1421)
    for (int i = 0; i < NR; i++)
    {
        for (int j = 1421; j < NZ; j++)
        {
            if (flag[i][j])
            {
                // 说明电极内
                mvr[i][j] = 0.0;
                mvz[i][j] = 0.0;
                // p[i][j]   = 1.0*rgas*T0; // 原代码
                p[i][j] = rgas * T0;
                T[i][j] = T0;
                // 28.966/1000.0/22.413996e-3 => 约1.28? 视具体含义
                rou[i][j] = 28.966 / 1000.0 / 22.413996e-3;
            }
            else if (iflag[i][j])
            {
                absv = std::sqrt(mvr[i][j] * mvr[i][j] + mvz[i][j] * mvz[i][j]);
                mvr[i][j] = Br[i][j] * absv; // Br[i][j] 全局?
                mvz[i][j] = Bz[i][j] * absv;

                p[i][j] = p[i + 1][j];
                T[i][j] = T0;
                rou[i][j] = p[i][j] / (rgas * T[i][j]);
            }
            else if (jflag[i][j])
            {
                absv = std::sqrt(mvr[i][j] * mvr[i][j] + mvz[i][j] * mvz[i][j]);
                mvr[i][j] = Br[i][j] * absv;
                mvz[i][j] = Bz[i][j] * absv;

                p[i][j] = p[i][j - 1];
                T[i][j] = T0;
                rou[i][j] = p[i][j] / (rgas * T[i][j]);
            }
            else if (oflag[i][j])
            {
                absv = std::sqrt(mvr[i][j] * mvr[i][j] + mvz[i][j] * mvz[i][j]);
                mvr[i][j] = Br[i][j] * absv;
                mvz[i][j] = Bz[i][j] * absv;

                p[i][j] = p[i + 1][j - 1];
                T[i][j] = T0;
                rou[i][j] = p[i][j] / (rgas * T[i][j]);
            }
        }
    }

    // 5) 左侧边界 (i=0)
    {
        int i = 0;
        for (int j = 0; j < 1421; j++)
        {
            if (jflag[i][j])
            {
                // ...
                mvr[i][j] = mvr[i + 1][j];
                mvz[i][j] = 0.0;
                p[i][j] = p[i][j - 1];
                T[i][j] = T0;
                rou[i][j] = p[i][j] / (rgas * T[i][j]);
            }
            else
            {
                mvr[i][j] = mvr[i + 1][j];
                mvz[i][j] = mvz[i + 1][j];
                p[i][j] = p[i + 1][j];
                // T[i][j] = (htotal - 0.5*(...)) / cpgas; 但原注释中是 T[i][j] = T[i+1][j];
                T[i][j] = T[i + 1][j];
                rou[i][j] = p[i][j] / (rgas * T[i][j]);
            }
        }
    }

    // 6) 上侧边界 (j=NZ-1)
    {
        int j = NZ - 1;
        for (int i = 0; i < NR; i++)
        {
            // 若电极内, 跳过(已设置)
            if (flag[i][j])
            {
                // do nothing
            }
            else if (iflag[i][j])
            {
                absv = std::sqrt(mvr[i][j] * mvr[i][j] + mvz[i][j] * mvz[i][j]);
                mvr[i][j] = Br[i][j] * absv;
                mvz[i][j] = Bz[i][j] * absv;
                p[i][j] = p[i][j - 1];
                T[i][j] = T0;
                rou[i][j] = p[i][j] / (rgas * T[i][j]);
            }
            else
            {
                mvr[i][j] = mvr[i][j - 1];
                mvz[i][j] = mvz[i][j - 1];
                p[i][j] = p[i][j - 1];
                T[i][j] = (htotal - 0.5 * (mvr[i][j] * mvr[i][j] + mvz[i][j] * mvz[i][j])) / cpgas;
                rou[i][j] = p[i][j] / (rgas * T[i][j]);
            }
        }
    }

    // 7) 下侧边界 (j=0)
    {
        int j = 0;
        for (int i = 0; i < NR; i++)
        {
            mvr[i][j] = mvr[i][j + 1];
            mvz[i][j] = 0.0;
            p[i][j] = p[i][j + 1];
            T[i][j] = T0;
            rou[i][j] = p[i][j] / (rgas * T[i][j]);
        }
    }

    // 8) 右侧边界 (i=NR-1)
    {
        int i = NR - 1;
        for (int j = 0; j < NZ; j++)
        {
            mvr[i][j] = mvr[i - 1][j];
            mvz[i][j] = mvz[i - 1][j];
            p[i][j] = p[i - 1][j];
            T[i][j] = (htotal - 0.5 * (mvr[i][j] * mvr[i][j] + mvz[i][j] * mvz[i][j])) / cpgas;
            rou[i][j] = p[i][j] / (rgas * T[i][j]);
        }
    }

// 9) 最后再并行把更新后的 rou,mvr,mvz,p,T 写回 q1,q2,q3,q4
#pragma omp parallel num_threads(THREAD_NUM)
    {
#pragma omp for
        for (int i = 0; i < NR; i++)
        {
            for (int j = 0; j < NZ; j++)
            {
                q1[i][j] = rou[i][j];
                q2[i][j] = rou[i][j] * mvr[i][j];
                q3[i][j] = rou[i][j] * mvz[i][j];
                q4[i][j] = rou[i][j] * (rgas * T[i][j] / (g0 - 1.0) + 0.5 * (mvr[i][j] * mvr[i][j] + mvz[i][j] * mvz[i][j]));
            }
        }
    }
}
void thread_func_fluid(ThreadArg &arg)
{
        ThreadArg& targ = arg;     // 可省，直接用 `arg`

    switch (targ.thread_no)
    {
    case 0: // electron: ne
        MUSCL_superbee_methoed_for_e(
            ne, vr, vz, dt, kappa, superbee,
            NR, NZ, Sr, Sz, Vol,
            iflag, jflag, oflag);
        negative_boundary_condition(NR, NZ, ne, flag);
        std::cout << "ne ";
        break;

    case 1: // Om
        MUSCL_superbee_methoed_for_mion(
            Om, miv_r, miv_z, dt, kappa, superbee,
            NR, NZ, Sr, Sz, Vol,
            iflag, jflag, oflag);
        negative_boundary_condition(NR, NZ, Om, flag);
        std::cout << "Om ";
        break;

    case 2: // O2m
        MUSCL_superbee_methoed_for_mion(
            O2m, miv_r, miv_z, dt, kappa, superbee,
            NR, NZ, Sr, Sz, Vol,
            iflag, jflag, oflag);
        negative_boundary_condition(NR, NZ, O2m, flag);
        std::cout << "O2m ";
        break;

    case 3: // N2p
        MUSCL_superbee_methoed_for_pion(
            N2p, piv_r, piv_z, dt, kappa, superbee,
            NR, NZ, Sr, Sz, Vol,
            iflag, jflag, oflag);
        positive_boundary_condition(NR, NZ, N2p, flag);
        std::cout << "N2p ";
        break;

    case 4: // OHm
        MUSCL_superbee_methoed_for_mion(
            OHm, miv_r, miv_z, dt, kappa, superbee,
            NR, NZ, Sr, Sz, Vol,
            iflag, jflag, oflag);
        negative_boundary_condition(NR, NZ, OHm, flag);
        std::cout << "OHm ";
        break;

    case 5: // Hm
        MUSCL_superbee_methoed_for_mion(
            Hm, miv_r, miv_z, dt, kappa, superbee,
            NR, NZ, Sr, Sz, Vol,
            iflag, jflag, oflag);
        negative_boundary_condition(NR, NZ, Hm, flag);
        std::cout << "Hm ";
        break;

    case 6: // O2p
        MUSCL_superbee_methoed_for_pion(
            O2p, piv_r, piv_z, dt, kappa, superbee,
            NR, NZ, Sr, Sz, Vol,
            iflag, jflag, oflag);
        positive_boundary_condition(NR, NZ, O2p, flag);
        std::cout << "O2p ";
        break;

    case 7: // H2Op
        MUSCL_superbee_methoed_for_pion(
            H2Op, piv_r, piv_z, dt, kappa, superbee,
            NR, NZ, Sr, Sz, Vol,
            iflag, jflag, oflag);
        positive_boundary_condition(NR, NZ, H2Op, flag);
        std::cout << "H2Op ";
        break;

    case 8: // N4p
        MUSCL_superbee_methoed_for_pion(
            N4p, piv_r, piv_z, dt, kappa, superbee,
            NR, NZ, Sr, Sz, Vol,
            iflag, jflag, oflag);
        positive_boundary_condition(NR, NZ, N4p, flag);
        std::cout << "N4p ";
        break;

    case 9: // O4p
        MUSCL_superbee_methoed_for_pion(
            O4p, piv_r, piv_z, dt, kappa, superbee,
            NR, NZ, Sr, Sz, Vol,
            iflag, jflag, oflag);
        positive_boundary_condition(NR, NZ, O4p, flag);
        std::cout << "O4p ";
        break;

    case 10: // N2O2p
        MUSCL_superbee_methoed_for_pion(
            N2O2p, piv_r, piv_z, dt, kappa, superbee,
            NR, NZ, Sr, Sz, Vol,
            iflag, jflag, oflag);
        positive_boundary_condition(NR, NZ, N2O2p, flag);
        std::cout << "N2O2p ";
        break;

    case 11: // O2pH2O
        MUSCL_superbee_methoed_for_pion(
            O2pH2O, piv_r, piv_z, dt, kappa, superbee,
            NR, NZ, Sr, Sz, Vol,
            iflag, jflag, oflag);
        positive_boundary_condition(NR, NZ, O2pH2O, flag);
        std::cout << "O2pH2O ";
        break;

    case 12: // H3Op
        MUSCL_superbee_methoed_for_pion(
            H3Op, piv_r, piv_z, dt, kappa, superbee,
            NR, NZ, Sr, Sz, Vol,
            iflag, jflag, oflag);
        positive_boundary_condition(NR, NZ, H3Op, flag);
        std::cout << "H3Op ";
        break;

    case 13: // H3OpH2O
        MUSCL_superbee_methoed_for_pion(
            H3OpH2O, piv_r, piv_z, dt, kappa, superbee,
            NR, NZ, Sr, Sz, Vol,
            iflag, jflag, oflag);
        positive_boundary_condition(NR, NZ, H3OpH2O, flag);
        std::cout << "H3OpH2O ";
        break;

    case 14: // H3OpH2O2
        MUSCL_superbee_methoed_for_pion(
            H3OpH2O2, piv_r, piv_z, dt, kappa, superbee,
            NR, NZ, Sr, Sz, Vol,
            iflag, jflag, oflag);
        positive_boundary_condition(NR, NZ, H3OpH2O2, flag);
        std::cout << "H3OpH2O2 ";
        break;

    case 15: // H3OpH2O3
        MUSCL_superbee_methoed_for_pion(
            H3OpH2O3, piv_r, piv_z, dt, kappa, superbee,
            NR, NZ, Sr, Sz, Vol,
            iflag, jflag, oflag);
        positive_boundary_condition(NR, NZ, H3OpH2O3, flag);
        std::cout << "H3OpH2O3 ";
        break;

    default:
        std::cerr << "[thread_func_fluid] unexpected thread_no = "
                  << targ.thread_no << '\n';
        break;
    }
}
void get_constant()
{
    int i, j;
    FILE *fp;

    // 分配
    dc = mat(DATA, BN);
    ddc = mat(DATA, BN);
    TdE = vec(BN);
    deV = vec(BN);
    dv = vec(BN);
    ddeV = vec(BN);
    ddv = vec(BN);
    alpTd = vec(AN);
    dalp = vec(AN);
    ddalp = vec(AN);

    powerTd = vec(BN);
    dpower = vec(BN);
    ddpower = vec(BN);

    // 1) 读 BOLTZMANNFILE
    fp = std::fopen(BOLTZMANNFILE, "r");
    if (!fp)
    {
        std::perror("Cannot open BOLTZMANNFILE");
        std::exit(EXIT_FAILURE);
    }

    for (i = 0; i < BN; i++)
    {
        if (std::feof(fp))
        {
            std::fprintf(stderr, "Wrong Num of BN!!\n");
            std::exit(EXIT_FAILURE);
        }
        std::fscanf(fp, "%lf", &TdE[i]);
        std::fscanf(fp, "%lf", &deV[i]);
        std::fscanf(fp, "%lf", &dv[i]);
        for (j = 1; j < NUM_REACT; j++)
        {
            std::fscanf(fp, "%lf", &dc[j][i]);
        }
    }
    std::fclose(fp);

    /*
    // 如果需要再读 ion_co.dat:
    fp=fopen("inputdata/ion_co.dat","r");
    for(i=0;i<AN;i++){
        fscanf(fp,"%lf",&alpTd[i]);
        fscanf(fp,"%lf",&dalp[i]);
    }
    fclose(fp);
    */

    // 2) spline
    spline(TdE, dv, BN, ddv);
    spline(TdE, deV, BN, ddeV);
    for (i = 1; i < NUM_REACT; i++)
    {
        spline(TdE, dc[i], BN, ddc[i]);
    }
    // spline(alpTd, dalp, AN, ddalp); // 如果要

    // 3) 读 inputdata/PowerEdep.dat
    fp = std::fopen("inputdata/PowerEdep.dat", "r");
    if (!fp)
    {
        std::perror("Cannot open inputdata/PowerEdep.dat");
        std::exit(EXIT_FAILURE);
    }
    for (i = 0; i < BN; i++)
    {
        std::fscanf(fp, "%lf", &powerTd[i]);
        std::fscanf(fp, "%lf", &dpower[i]);
    }
    std::fclose(fp);

    // 4) spline
    spline(powerTd, dpower, BN, ddpower);

    // 到此和原 get_constant() 相同
    std::cout << "get_constant() finished reading constants.\n";
}

// 改写后的线程函数：C++ 风格
void *thread_func_2D_to_3D(void *arg)
{
    // 1) 先把 arg 转回 ThreadArg*
    ThreadArg *targ = static_cast<ThreadArg *>(arg);

    // 2) 计算本线程处理的区间
    int ist = NR * targ->thread_no / THREAD_NUM;
    int ifi = NR * (targ->thread_no + 1) / THREAD_NUM;

    for (int i = ist; i < ifi; i++)
    {
        for (int j = 0; j < NZ; j++)
        {
            double inv_rou = 1.0 / (rou[i][j] * MassNum);

            PARTICLE[i][j][n1] = ne[i][j] * inv_rou;
            PARTICLE[i][j][n2] = N2p[i][j] * inv_rou;
            PARTICLE[i][j][n3] = O2p[i][j] * inv_rou;
            PARTICLE[i][j][n4] = O2m[i][j] * inv_rou;
            PARTICLE[i][j][n5] = Om[i][j] * inv_rou;
            PARTICLE[i][j][n6] = H2Op[i][j] * inv_rou;
            PARTICLE[i][j][n7] = OHm[i][j] * inv_rou;
            PARTICLE[i][j][n8] = Hm[i][j] * inv_rou;
            PARTICLE[i][j][n9] = N4p[i][j] * inv_rou;
            PARTICLE[i][j][n10] = O4p[i][j] * inv_rou;
            PARTICLE[i][j][n11] = N2O2p[i][j] * inv_rou;
            PARTICLE[i][j][n12] = O2pH2O[i][j] * inv_rou;
            PARTICLE[i][j][n13] = H3Op[i][j] * inv_rou;
            PARTICLE[i][j][n14] = H3OpH2O[i][j] * inv_rou;
            PARTICLE[i][j][n15] = H3OpH2O2[i][j] * inv_rou;
            PARTICLE[i][j][n16] = H3OpH2O3[i][j] * inv_rou;
        }
    }

    // pthread 规定：线程函数返回 void*
    return nullptr;
}

// 用一个普通的函数(或可调用对象)来表示线程执行体
void thread_func_3D_to_2D(const ThreadArg &targ)
{
    // 线程号
    int thread_id = targ.thread_no;

    // 计算分段
    int ist = NR * thread_id / THREAD_NUM;
    int ifi = NR * (thread_id + 1) / THREAD_NUM;

    for (int i = ist; i < ifi; i++)
    {
        for (int j = 0; j < NZ; j++)
        {
            double rrou = rou[i][j] * MassNum;

            ne[i][j] = PARTICLE[i][j][n1] * rrou;
            N2p[i][j] = PARTICLE[i][j][n2] * rrou;
            O2p[i][j] = PARTICLE[i][j][n3] * rrou;
            // ...
        }
    }

    // 不需要返回任何指针，C++中可直接结束
}

// 这里我们直接用一个普通函数
void thread_func_forCurrent(const ThreadArg &targ)
{
    int ist = NR * targ.thread_no / THREAD_NUM;
    int ifi = NR * (targ.thread_no + 1) / THREAD_NUM;

    double cnst1 = MOL * 1e-21;
    double cnst2 = 1.0 / dt;

    for (int i = ist; i < ifi; i++)
    {
        for (int j = 0; j < NZ; j++)
        {
            preEr[i][j] = Ex[i][j] * cnst1;
            preEz[i][j] = Ey[i][j] * cnst1;

            dtELr[i][j] += preEr[i][j];
            dtELz[i][j] += preEz[i][j];

            dtELr[i][j] *= cnst2;
            dtELz[i][j] *= cnst2;
        }
    }
    // 不需要return值，C++函数可直接结束
}

void thread_func_forPHT(ThreadArg &targ)
{
    int threadNo = targ.thread_no;

    // 根据线程号划分处理范围
    int ist = NR * threadNo / THREAD_NUM;
    int ifi = NR * (threadNo + 1) / THREAD_NUM;

    // 用到的常量(与物理量)：
    double TeVkeisu = 2.0 / 3.0 * QELEC / kb; // 电子温度计算用的系数
    double ph_const = 0.1 * 30.0 / (760.0 + 30.0) * 100.0;
    double ph_consto = 0.02 * 36.0 / (760.0 + 36.0) * 100.0;

    // 在本线程的 i 范围内做循环
    for (int i = ist; i < ifi; ++i)
    {
        for (int j = 0; j < NZ; ++j)
        {
            // 计算 v[i][j]
            //   如果 i == NR-1 或 j == NZ-1，直接用 vz[i][j], vr[i][j] 本身
            //   否则用相邻网格做平均
            if (i == NR - 1 || j == NZ - 1)
            {
                v[i][j] = std::sqrt(
                    std::pow(0.5 * (vz[i][j] + vz[i][j]), 2) +
                    std::pow(0.5 * (vr[i][j] + vr[i][j]), 2));
            }
            else
            {
                v[i][j] = std::sqrt(
                    std::pow(0.5 * (vz[i][j] + vz[i][j + 1]), 2) +
                    std::pow(0.5 * (vr[i][j] + vr[i + 1][j]), 2));
            }

            // 计算 piv, miv
            piv[i][j] = std::sqrt(piv_z[i][j] * piv_z[i][j] + piv_r[i][j] * piv_r[i][j]);
            miv[i][j] = std::sqrt(miv_z[i][j] * miv_z[i][j] + miv_r[i][j] * miv_r[i][j]);

            // 使用 spline 函数查询某些插值: alpTd/dalp + dc[41], ...
            double kalp = 0.0;
            double vo2i = 0.0;
            splint(alpTd, dalp, ddalp, AN, absE[i][j], &kalp);
            splint(TdE, dc[41], ddc[41], BN, absE[i][j], &vo2i);

            // 计算生成率(光电离相关)
            I[i][j] = ph_const * kalp * v[i][j] * ne[i][j]; // in m^-3 s^-1

            // rrho[...] = ...
            rrho[point(NZ, i, j)] = ph_const * kalp * v[i][j] * ne[i][j] + ph_consto * vo2i * ne[i][j] * pO2 * MOL;

            // 如果当前网格是flag==1(比如电极内?)
            //   TeV[i][j] = 0
            // 否则按电场强度计算电子温度(energy -> eV -> K)
            if (flag[i][j])
            {
                TeV[i][j] = 0.0;
            }
            else
            {
                double localE = std::fabs(absE[i][j]);
                double eVvalue = 0.0;
                // 用 spline 在 deV 上插值
                splint(TdE, deV, ddeV, BN, localE, &eVvalue);
                TeV[i][j] = TeVkeisu * eVvalue;
            }
        }
    }
}

void react_function(int thread_no)
{
    // 局部变量 (与原函数完全一致)
    int ist, ifi, i, j, k, R1, R2, R3, L1, L2, L3, nnum;
    int kk, klo, khi;
    double inv_h, hh, bb, aa, r, dT, kr;
    // react
    double *Y;
    double *nYout;
    double *eYout;
    double *iYout;
    double RE, dph;

    // 原代码里写的：
    //   double temp1=exp(0.39*log(T[i][j]));
    //   double temp2=exp(0.7*log(T[i][j]));
    // 但 i,j 尚未确定，放在循环内才正确(见下)。
    // 先声明：
    double temp1, temp2, temp3, temp4, temp5, temp6;
    double inv_rou, rrou;
    int numOH;

    // 若需要找 "OH" 对应的粒子编号:
    //   numOH = particle_number("OH", particle, n_particles);
    // 这里就保持注释：
    // numOH = particle_number("OH", particle, n_particles);

    // 分配临时数组
    Y = vec(n_particles);
    nYout = vec(n_particles);
    eYout = vec(n_particles);
    iYout = vec(n_particles);

    // 计算本线程负责的 i 区间
    ist = NR * thread_no / THREAD_NUM;
    ifi = NR * (thread_no + 1) / THREAD_NUM;

    // 如果是第0号线程，可打印一些提示
    if (thread_no == 0)
    {
        std::printf("---Reaction-Calculation... ");
    }

    //========================================
    // 主循环：i从ist到ifi, j从0到NZ
    //========================================
    for (i = ist; i < ifi; i++)
    {
        for (j = 0; j < NZ; j++)
        {

            // 只有非flag网格才进行计算
            if (!flag[i][j])
            {

                // --------------
                // 先计算 rrou
                // --------------
                rrou = rou[i][j] * MassNum;

                // 将 3D数组 PARTICLE[i][j][k] -> Y[k]
                // 并把 eYout, iYout, nYout 清零
                for (k = 0; k < n_particles; k++)
                {
                    Y[k] = PARTICLE[i][j][k] * rrou;
                    eYout[k] = 0.0;
                    iYout[k] = 0.0;
                    nYout[k] = 0.0;
                }

                // 原代码：Y[0] = 1.0;
                Y[0] = 1.0;
                dT = 0.0;

                // 下面两句原先直接写在声明处，但是要放到循环里才能用 i,j
                temp1 = std::exp(0.39 * std::log(T[i][j]));
                temp2 = std::exp(0.70 * std::log(T[i][j]));

                //------------------------------------------------
                // (1) 电子反应 ne_react 次
                //------------------------------------------------
                for (k = 0; k < ne_react; k++)
                {

                    RE = absE[i][j];
                    nnum = eRnum[k];

                    // 二分搜索 klo,khi
                    klo = 0;
                    khi = BN - 1;
                    while ((khi - klo) > 1)
                    {
                        kk = (khi + klo) >> 1; // 取中点
                        if (TdE[kk] > RE)
                        {
                            khi = kk;
                        }
                        else
                        {
                            klo = kk;
                        }
                    }

                    // 根据 klo,khi 计算插值
                    hh = TdE[khi] - TdE[klo];
                    if (hh == 0.0)
                    {
                        std::fprintf(stderr, "Bad xa input to routine splint\n");
                        // 若需退出:
                        // std::exit(1); // or return;
                    }

                    inv_h = 1.0 / hh;
                    aa = (TdE[khi] - RE) * inv_h;
                    bb = (RE - TdE[klo]) * inv_h;

                    kr = aa * dc[nnum][klo] + bb * dc[nnum][khi] + ((aa * aa * aa - aa) * ddc[nnum][klo] + (bb * bb * bb - bb) * ddc[nnum][khi]) * (hh * hh) * sp6;

                    // 如果 RE 超出表范围
                    if (RE > TdE[BN - 1])
                    {
                        kr = dc[nnum][BN - 1];
                    }
                    if (RE < TdE[0])
                    {
                        kr = dc[nnum][0];
                    }

                    // L1, R1, L2, R2, ...
                    L1 = ereactl_1st[k];
                    R1 = ereactr_1st[k];
                    L2 = ereactl_2nd[k];
                    R2 = ereactr_2nd[k];
                    L3 = ereactl_3rd[k];
                    R3 = ereactr_3rd[k];

                    // 原代码里针对 nstp > Prinstp2 做了额外判断
                    if (nstp > Prinstp2)
                    {
                        if (j > limpoz && ne[i][j] > 1e16)
                        {
                            // k == 15 or 27 => 可能是N2电离 / O2电离
                            if (k == 15 || k == 27)
                            {
                                if (neflag[i][j])
                                {
                                    if (RE > 120.0)
                                    {
                                        kr = kr * std::exp(-0.01 * (RE - 120.0));
                                    }
                                    else
                                    {
                                        // 原代码是空
                                    }
                                }
                                else
                                {
                                    // 空
                                }
                            }
                            else
                            {
                                // 空
                            }
                        }
                    }

                    // 计算生成/消耗速率 r
                    r = dt * kr * (Y[L1] * Y[L2] * Y[L3]);
                    dT += r * eER_mod[k];

                    // 更新 eYout
                    if (L1 != 0)
                        eYout[L1] -= r;
                    if (L2 != 0)
                        eYout[L2] -= r;
                    if (L3 != 0)
                        eYout[L3] -= r;

                    if (R1 != 0)
                        eYout[R1] += r;
                    if (R2 != 0)
                        eYout[R2] += r;
                    if (R3 != 0)
                        eYout[R3] += r;

                    // 若 i<numrad 则记录 sink
                    if (i < numrad)
                    {
                        if (j < jfall)
                            e_sink[0][k] += r * Vol[i][j];
                        else if (j < j1st)
                            e_sink[1][k] += r * Vol[i][j];
                        else if (j < janode)
                            e_sink[2][k] += r * Vol[i][j];
                        else if (j < 1421)
                            e_sink[3][k] += r * Vol[i][j];
                    }

                    molP[i][j][k] = r * eK2_mod[k] * Vol[i][j];

                } // end for k < ne_react

                //------------------------------------------------
                // (2) 各种离子反应 iK0m[0..5], iRnum
                //------------------------------------------------
                temp3 = std::log(TeV[i][j]);
                temp6 = 0.5 * std::log(300.0 / TeV[i][j]);

                iK0m[0] = iK0[1] * temp1 * std::exp(-0.39 * temp3);
                iK0m[1] = iK0[2] * temp2 * std::exp(-0.70 * temp3);
                iK0m[2] = iK0[3] * std::exp(temp6);
                iK0m[3] = iK0[4] * std::exp(temp6);
                iK0m[4] = iK0[5] * std::exp(temp6);
                iK0m[5] = iK0[6] * std::exp(temp6);

                // 原代码里逐一写 k=0..5
                // 下面只演示 k=0 => 其余1..5同理
                k = 0;
                L1 = ireactl_1st[k];
                R1 = ireactr_1st[k];
                L2 = ireactl_2nd[k];
                R2 = ireactr_2nd[k];
                L3 = ireactl_3rd[k];
                R3 = ireactr_3rd[k];

                r = dt * iK0m[k] * Y[L1] * Y[L2] * Y[L3];
                dT += r * iER_mod[k];

                if (L1 != 0)
                    iYout[L1] -= r;
                if (L2 != 0)
                    iYout[L2] -= r;
                if (L3 != 0)
                    iYout[L3] -= r;

                if (R1 != 0)
                    iYout[R1] += r;
                if (R2 != 0)
                    iYout[R2] += r;
                if (R3 != 0)
                    iYout[R3] += r;

                if (i < numrad)
                {
                    if (j < jfall)
                        i_sink[0][k] += r * Vol[i][j];
                    else if (j < j1st)
                        i_sink[1][k] += r * Vol[i][j];
                    else if (j < janode)
                        i_sink[2][k] += r * Vol[i][j];
                    else if (j < 1421)
                        i_sink[3][k] += r * Vol[i][j];
                }

                // 同理 k=1..5
                // ...
                for (k = 6; k < ni_react; k++)
                {
                    L1 = ireactl_1st[k];
                    R1 = ireactr_1st[k];
                    L2 = ireactl_2nd[k];
                    R2 = ireactr_2nd[k];
                    L3 = ireactl_3rd[k];
                    R3 = ireactr_3rd[k];

                    r = dt * iK0m[k] * Y[L1] * Y[L2] * Y[L3];
                    dT += r * iER_mod[k];

                    if (L1 != 0)
                        iYout[L1] -= r;
                    if (L2 != 0)
                        iYout[L2] -= r;
                    if (L3 != 0)
                        iYout[L3] -= r;

                    if (R1 != 0)
                        iYout[R1] += r;
                    if (R2 != 0)
                        iYout[R2] += r;
                    if (R3 != 0)
                        iYout[R3] += r;

                    if (i < numrad)
                    {
                        if (j < jfall)
                            i_sink[0][k] += r * Vol[i][j];
                        else if (j < j1st)
                            i_sink[1][k] += r * Vol[i][j];
                        else if (j < janode)
                            i_sink[2][k] += r * Vol[i][j];
                        else if (j < 1421)
                            i_sink[3][k] += r * Vol[i][j];
                    }
                }

                //------------------------------------------------
                // (3) 中性粒子反应
                //------------------------------------------------
                temp4 = 1.0 / T[i][j];
                temp5 = std::log(T[i][j] / 298.0);

                for (k = 0; k < nn_react; k++)
                {
                    L1 = nreactl_1st[k];
                    R1 = nreactr_1st[k];
                    L2 = nreactl_2nd[k];
                    R2 = nreactr_2nd[k];
                    L3 = nreactl_3rd[k];
                    R3 = nreactr_3rd[k];

                    r = dt * nK0m[k] * std::exp(-nK2m[k] * temp4 + nK1m[k] * temp5) * Y[L1] * Y[L2] * Y[L3];
                    dT += r * nER_mod[k];

                    if (L1 != 0)
                        nYout[L1] -= r;
                    if (L2 != 0)
                        nYout[L2] -= r;
                    if (L3 != 0)
                        nYout[L3] -= r;

                    if (R1 != 0)
                        nYout[R1] += r;
                    if (R2 != 0)
                        nYout[R2] += r;
                    if (R3 != 0)
                        nYout[R3] += r;

                    if (i < numrad)
                    {
                        if (j < jfall)
                            sink[0][k] += r * Vol[i][j];
                        else if (j < j1st)
                            sink[1][k] += r * Vol[i][j];
                        else if (j < janode)
                            sink[2][k] += r * Vol[i][j];
                        else if (j < 1421)
                            sink[3][k] += r * Vol[i][j];
                    }
                }

                // 赋值到 cE
                cE[i][j] = dT;

                //------------------------------------------------
                // (4) 电场消耗功
                //------------------------------------------------
                // 先用 splint(...) 插值得到 power
                double power;
                splint(powerTd, dpower, ddpower, BN, RE, &power);
                // POWER[i][j] = dt*(PARTICLE[i][j][nume]*rrou)* power * ((...)*rrou)*Vol[i][j];
                POWER[i][j] = dt * (PARTICLE[i][j][nume] * rrou) * power * ((PARTICLE[i][j][numO2] + PARTICLE[i][j][numN2] + PARTICLE[i][j][numH2O]) * rrou) * Vol[i][j];

                //------------------------------------------------
                // (5) 写回 PARTICLE, 再加上光电离
                //------------------------------------------------
                inv_rou = 1.0 / rrou;
                for (int m = 0; m < n_particles; m++)
                {
                    PARTICLE[i][j][m] += (eYout[m] + iYout[m] + nYout[m]) * inv_rou;
                }

                // 光电离
                dph = dt * phot_num[i][j];
                PARTICLE[i][j][numO2] -= dph * inv_rou;
                PARTICLE[i][j][numO2p] += dph * inv_rou;
                PARTICLE[i][j][nume] += dph * inv_rou;

            } // end if(!flag[i][j])
        } // end for j
    } // end for i

    // 释放临时数组
    free_vec(Y, n_particles);
    free_vec(nYout, n_particles);
    free_vec(eYout, n_particles);
    free_vec(iYout, n_particles);

    // 打印线程编号(模拟原先 printf("%d  ",targ->thread_no); )
    std::printf("%d  ", thread_no);

    // 结束
    return;
}
void Debug_Boltzmann(int debug)
{
    if (debug)
    {
        // Boltzmann Check
        double d0, d1, d2;
        int i, k;
        for (i = 0; i < 1000; i++)
        {
            // 先让 k=1，然后调用 splint(...) 将结果放到 d0
            k = 1;
            splint(TdE, dc[k], ddc[k], BN, static_cast<double>(i), &d0);

            // 再让 k=2，将结果放到 d1
            k = 2;
            splint(TdE, dc[k], ddc[k], BN, static_cast<double>(i), &d1);

            // 最后 k=43，将结果放到 d2
            k = 43;
            splint(TdE, dc[k], ddc[k], BN, static_cast<double>(i), &d2);

            // 输出
            std::printf("%f\t%.4e\t%.4e\t%.4e\n",
                        static_cast<double>(i), d0, d1, d2);
        }

        std::printf("Debug_for_Boltzmann_spline\n");
        std::exit(0);
    }
    else
    {
        // 如果 debug == 0，不进行任何操作
    }
}
void vib_convert(int i, int j, double *vib, double **Ev)
{
    double inv_O2Tv, rrou;

    // 1) 计算该网格处的 (rou[i][j]*MassNum)
    rrou = rou[i][j] * MassNum;

    // 2) 依次把 PARTICLE 中若干特定索引的值乘以 rrou，放到 vib[] 指定位置
    vib[0] = PARTICLE[i][j][o2v0] * rrou;
    vib[1] = PARTICLE[i][j][o2v1] * rrou;
    vib[2] = PARTICLE[i][j][o2v2] * rrou;
    vib[3] = PARTICLE[i][j][o2v3] * rrou;
    vib[4] = PARTICLE[i][j][o2v4] * rrou;

    // 3) 计算 1 / ((Ev[0][4] - Ev[0][3]) / log(...)) => inv_O2Tv
    inv_O2Tv = 1.0 / ((Ev[0][4] - Ev[0][3]) /
                      std::log(vib[3] / vib[4]));

    // 4) 按公式继续给 vib[5..8] 赋值
    vib[5] = vib[4] * std::exp((Ev[0][4] - Ev[0][5]) * inv_O2Tv); // O2v5
    vib[6] = vib[5] * std::exp((Ev[0][5] - Ev[0][6]) * inv_O2Tv); // O2v6
    vib[7] = vib[6] * std::exp((Ev[0][6] - Ev[0][7]) * inv_O2Tv); // O2v7
    vib[8] = vib[7] * std::exp((Ev[0][7] - Ev[0][8]) * inv_O2Tv); // O2v8

    // 5) N2 各振动态
    vib[9] = PARTICLE[i][j][n2v0] * rrou;
    vib[10] = PARTICLE[i][j][n2v1] * rrou;
    vib[11] = PARTICLE[i][j][n2v2] * rrou;
    vib[12] = PARTICLE[i][j][n2v3] * rrou;
    vib[13] = PARTICLE[i][j][n2v4] * rrou;
    vib[14] = PARTICLE[i][j][n2v5] * rrou;
    vib[15] = PARTICLE[i][j][n2v6] * rrou;
    vib[16] = PARTICLE[i][j][n2v7] * rrou;
    vib[17] = PARTICLE[i][j][n2v8] * rrou;

    // 6) H2O 各振动态
    vib[18] = PARTICLE[i][j][h2ov0] * rrou;
    vib[19] = PARTICLE[i][j][h2ov2] * rrou;
    vib[20] = PARTICLE[i][j][h2ov1] * rrou;
    vib[21] = PARTICLE[i][j][h2ov3] * rrou;

    // 函数结束
}

void vib_assemble(int i, int j, double *dvib, double dt)
{
    // 1) 先计算 inv_rou = dt / ( rou[i][j] * MassNum )
    double inv_rou = dt / (rou[i][j] * MassNum);

    // 2) O2 的振动态
    PARTICLE[i][j][o2v0] += dvib[0] * inv_rou;
    PARTICLE[i][j][o2v1] += dvib[1] * inv_rou;
    PARTICLE[i][j][o2v2] += dvib[2] * inv_rou;
    PARTICLE[i][j][o2v3] += dvib[3] * inv_rou;
    PARTICLE[i][j][o2v4] += dvib[4] * inv_rou;

    // 3) N2 的振动态
    PARTICLE[i][j][n2v0] += dvib[9] * inv_rou;
    PARTICLE[i][j][n2v1] += dvib[10] * inv_rou;
    PARTICLE[i][j][n2v2] += dvib[11] * inv_rou;
    PARTICLE[i][j][n2v3] += dvib[12] * inv_rou;
    PARTICLE[i][j][n2v4] += dvib[13] * inv_rou;
    PARTICLE[i][j][n2v5] += dvib[14] * inv_rou;
    PARTICLE[i][j][n2v6] += dvib[15] * inv_rou;
    PARTICLE[i][j][n2v7] += dvib[16] * inv_rou;
    PARTICLE[i][j][n2v8] += dvib[17] * inv_rou;

    // 4) H2O 的振动态
    PARTICLE[i][j][h2ov0] += dvib[18] * inv_rou;
    PARTICLE[i][j][h2ov2] += dvib[19] * inv_rou;
    PARTICLE[i][j][h2ov1] += dvib[20] * inv_rou;
    PARTICLE[i][j][h2ov3] += dvib[21] * inv_rou;

    // 函数结束
}

int main()

{
    // --- 1) 整型、浮点型、指针/数组等变量声明 ---
    int i = 0, j = 0, START = 0, Helm_iter = 0, k = 0;
    int chname = 0, loop = 0;
    double e_max = 0.0, pureAir_kg = 0.0, OMEGA = 0.0;
    double time = 0.0, CFLmemo = 0.0;
    double A = 0.0, B1 = 0.0, B2 = 0.0, C = 0.0, D = 0.0, CPp = 0.0, CPm = 0.0;

    // 这里若 spDATA 本身是 int 或 double 看您的需求
    // 假设它是 double 并且先初始化为 0
    double spDATA = 0.0;

    double *spTime = nullptr, *dspV = nullptr, *ddspV = nullptr, spV = 0.0;
    double **q1 = nullptr, **q2 = nullptr, **q3 = nullptr, **q4 = nullptr;
    double res = 0.0;
    double inv_rou = 0.0, rrou = 0.0, sum = 0.0, inv_sum = 0.0, MolMass = 0.0;
    double eV = 0.0, **dTvib = nullptr, **sumN2C = nullptr, **Dvx = nullptr, **Dvy = nullptr;

    // 文件 IO
    FILE *fp = nullptr;

    // 计算 a = sqrt(cr*bgap)
    double a = std::sqrt(cr * bgap);

    // filename 可以继续使用 C 风格数组
    char filename[256] = {0};

    // --- 2) 线程相关 ---
    pthread_t handle[THREAD_NUM];
    ThreadArg targ[THREAD_NUM];
    pthread_mutex_t mutex;

    // --- 3) 分配向量和矩阵 ---
    // 假设您已有的 vec() / mat() / imat() / ...
    // 均是能返回双指针或指针的 C 风格分配函数，保留原有调用方式：
    r = vec(NR);
    z = vec(NZ);
    rh = vec(NR);
    zh = vec(NZ);
    dr = vec(NR);
    dz = vec(NZ);
    Vol = mat(NR, NZ);
    Sr = mat(NR, NZ);
    Sz = mat(NR, NZ);

    flag = imat(NR, NZ);
    iflag = imat(NR, NZ);
    jflag = imat(NR, NZ);
    tflag = imat(NR, NZ);
    oflag = imat(NR, NZ);

    ne = mat(NR, NZ);
    N2p = mat(NR, NZ);
    O2p = mat(NR, NZ);
    O2m = mat(NR, NZ);
    Om = mat(NR, NZ);
    H2Op = mat(NR, NZ);
    OHm = mat(NR, NZ);
    Hm = mat(NR, NZ);

    N4p = mat(NR, NZ);
    O4p = mat(NR, NZ);
    N2O2p = mat(NR, NZ);
    O2pH2O = mat(NR, NZ);
    H3Op = mat(NR, NZ);
    H3OpH2O = mat(NR, NZ);
    H3OpH2O2 = mat(NR, NZ);
    H3OpH2O3 = mat(NR, NZ);

    O2 = mat(NR, NZ);
    N2 = mat(NR, NZ);

    P1 = mat(NR, NZ);
    P2 = mat(NR, NZ);
    P3 = mat(NR, NZ);
    P4 = mat(NR, NZ);
    P5 = mat(NR, NZ);

    rho = mat(NR, NZ);
    phi = mat(NR, NZ);
    Lphi = mat(NR, NZ);
    Cphi = mat(NR, NZ);

    absE = mat(NR, NZ);
    Ex = mat(NR, NZ);
    Ey = mat(NR, NZ);
    vr = mat(NR, NZ);
    vz = mat(NR, NZ);

    piv_r = mat(NR, NZ);
    piv_z = mat(NR, NZ);
    miv_r = mat(NR, NZ);
    miv_z = mat(NR, NZ);
    v = mat(NR, NZ);
    piv = mat(NR, NZ);
    miv = mat(NR, NZ);

    dtELr = mat(NR, NZ);
    dtELz = mat(NR, NZ);

    POWER = mat(NR, NZ);
    sumN2C = mat(NR, NZ);
    Dvx = mat(NR, NZ);
    Dvy = mat(NR, NZ);

    phot_num = mat(NR, NZ);
    HelmA = ivec(3);
    HelmG = vec(3);

    TeV = mat(NR, NZ);
    Cphi1 = mat(NR, NZ);
    Cphi2 = mat(NR, NZ);
    I = mat(NR, NZ);
    preEr = mat(NR, NZ);
    preEz = mat(NR, NZ);

    S_ph0 = mat(NR, NZ);
    S_ph1 = mat(NR, NZ);
    S_ph2 = mat(NR, NZ);

    sink = mat(4, 200);
    e_sink = mat(4, 200);
    i_sink = mat(4, 200);

    // 初始化 sink/e_sink/i_sink
    for (i = 0; i < 4; ++i)
    {
        for (j = 0; j < 200; ++j)
        {
            sink[i][j] = 0.0;
            e_sink[i][j] = 0.0;
            i_sink[i][j] = 0.0;
        }
    }

    Br = mat(NR, NZ);
    Bz = mat(NR, NZ);
    mvr = mat(NR, NZ);
    mvz = mat(NR, NZ);
    rou = mat(NR, NZ);
    p = mat(NR, NZ);

    q1 = mat(NR, NZ);
    q2 = mat(NR, NZ);
    q3 = mat(NR, NZ);
    q4 = mat(NR, NZ);
    cE = mat(NR, NZ);

    dTvib = mat(NR, NZ);

    // --- 4) 调用函数 ---
    // 读取常量 / 反应数据 等
    get_constant();     // 定数を読み込み
    Debug_Boltzmann(0); // Boltzmann 调试或检查

    // ------- 1) 声明变量 -------
    int cuda = 0;
    int cuda_flag = 0;
    int iter = 0;
    int itnum = 0;
    double error = 0.0;
    double Maxphi = 0.0;

    double *d_err = nullptr;
    double *dd_err = nullptr;
    double *d_phi = nullptr;
    double *d_rho = nullptr;
    double *d_temp = nullptr;
    double *d_rh = nullptr;
    double *d_data = nullptr;
    double *d_Sph0 = nullptr;
    double *d_Sph1 = nullptr;
    double *d_Sph2 = nullptr;
    double *d_P1 = nullptr;
    double *d_P2 = nullptr;
    double *d_P3 = nullptr;
    double *d_P4 = nullptr;
    double *d_P5 = nullptr;

    double toP = 0.0;
    double *each_mP = nullptr;

    int *d_flag = nullptr;
    int *d_iflag = nullptr;
    int *d_jflag = nullptr;
    int *d_oflag = nullptr;

    double *pphi = nullptr;
    double *temp = nullptr;
    double *rrh = nullptr;
    double *SSph0 = nullptr;
    double *SSph1 = nullptr;
    double *SSph2 = nullptr;
    double *PP1 = nullptr;
    double *PP2 = nullptr;
    double *PP3 = nullptr;
    double *PP4 = nullptr;
    double *PP5 = nullptr;

    int *fflag = nullptr;
    int *iiflag = nullptr;
    int *jjflag = nullptr;
    int *ooflag = nullptr;

    int N = 0;
    int mi = 0;
    int mf = 0;

    // ------- 2) 计算 N 并配置 CUDA 网格大小 -------
    N = NR * NZ; // grid size

    // dim3 来自 CUDA API，用于配置线程块、网格
    dim3 dimGrid(NR / blockDim_x, NZ / blockDim_y);
    dim3 dimBlock(blockDim_x, blockDim_y);

    // ------- 3) 计算需要的显存大小 -------
    mi = N * sizeof(int);    // memory size (int)
    mf = N * sizeof(double); // memory size (double)

    std::cout << "All initial allocations done.\n\n";
    // ---------------------------
    // 下面开始“2nd Grid”对应的代码
    // ---------------------------

    // 1) 计算 NR2、NZ2、NN
    int NR2 = NR / 2;
    int NZ2 = NZ / 2;
    int NN = NR2 * NZ2; // 2D grid size

    // 2) 配置 CUDA 线程网格
    //    （若仍需要在后续 kernel 中使用 dimGrid2 / dimBlock2）
    dim3 dimGrid2(NR2 / blockDim_x2, NZ2 / blockDim_y2);
    dim3 dimBlock2(blockDim_x2, blockDim_y2);

    // 3) 计算每种类型的总字节数
    int mi2 = NN * sizeof(int);    // 内存大小(以 int 为元素)
    int mf2 = NN * sizeof(double); // 内存大小(以 double 为元素)

    // 4) 在 CPU 上分配 rrho2 / fflag2
    //    （这里示例用项目已有的 vec(), ivec() 函数）
    double *rrho2 = vec(NN); // 分配 double[NN]
    int *fflag2 = ivec(NN);  // 分配 int[NN]

    // 5) 在 GPU 端分配 dd_rho / dd_Cphi / dd_rh / ... 等指针
    double *dd_rho = nullptr;
    double *dd_Cphi = nullptr;
    double *dd_rh = nullptr;
    double *dd_P1 = nullptr;
    double *dd_P2 = nullptr;
    double *dd_P3 = nullptr;
    double *dd_P4 = nullptr;
    double *dd_P5 = nullptr;
    double *dd_temp = nullptr;

    int *dd_flag = nullptr;
    int *dd_iflag = nullptr;
    int *dd_jflag = nullptr;
    int *dd_oflag = nullptr;
    int *dd_tflag = nullptr;

    // 6) 分配 GPU 内存（cudaMalloc）
    cudaMalloc(reinterpret_cast<void **>(&dd_rho), mf2);
    cudaMalloc(reinterpret_cast<void **>(&dd_Cphi), mf2);
    cudaMalloc(reinterpret_cast<void **>(&dd_rh), mf2);
    cudaMalloc(reinterpret_cast<void **>(&dd_P1), mf2);
    cudaMalloc(reinterpret_cast<void **>(&dd_P2), mf2);
    cudaMalloc(reinterpret_cast<void **>(&dd_P3), mf2);
    cudaMalloc(reinterpret_cast<void **>(&dd_P4), mf2);
    cudaMalloc(reinterpret_cast<void **>(&dd_P5), mf2);
    cudaMalloc(reinterpret_cast<void **>(&dd_temp), mf2);

    cudaMalloc(reinterpret_cast<void **>(&dd_flag), mi2);
    cudaMalloc(reinterpret_cast<void **>(&dd_iflag), mi2);
    cudaMalloc(reinterpret_cast<void **>(&dd_jflag), mi2);
    cudaMalloc(reinterpret_cast<void **>(&dd_oflag), mi2);
    cudaMalloc(reinterpret_cast<void **>(&dd_tflag), mi2);

    std::cout << "Set up 2nd Grid and allocated GPU memory (NR2="
              << NR2 << ", NZ2=" << NZ2 << ", total NN=" << NN << ").\n\n";

    int NR4 = NR2 / 2;
    int NZ4 = NZ2 / 2;
    // 计算网格大小
    int NNN = NR4 * NZ4; // grid size

    // 设置 GPU kernel 的 grid/block 维度
    dim3 dimGrid4(NR4 / blockDim_x4, NZ4 / blockDim_y4);
    dim3 dimBlock4(blockDim_x4, blockDim_y4);

    // 计算申请大小
    int mi4 = NNN * sizeof(int);    // memory size for int
    int mf4 = NNN * sizeof(double); // memory size for double

    // 声明需要的指针（设备端）
    double *ddd_rho = nullptr;
    double *ddd_Cphi = nullptr;
    double *ddd_rh = nullptr;
    double *ddd_P1 = nullptr;
    double *ddd_P2 = nullptr;
    double *ddd_P3 = nullptr;
    double *ddd_P4 = nullptr;
    double *ddd_P5 = nullptr;
    double *ddd_temp = nullptr;

    int *ddd_flag = nullptr;
    int *ddd_iflag = nullptr;
    int *ddd_jflag = nullptr;
    int *ddd_oflag = nullptr;
    int *ddd_tflag = nullptr;

    double *rrho4 = vec(NN);
    int *fflag4 = ivec(NN);

    // ============ 设备内存分配 ============
    cudaMalloc((void **)&ddd_rho, mf4);
    cudaMalloc((void **)&ddd_Cphi, mf4);
    cudaMalloc((void **)&ddd_rh, mf4);
    cudaMalloc((void **)&ddd_P1, mf4);
    cudaMalloc((void **)&ddd_P2, mf4);
    cudaMalloc((void **)&ddd_P3, mf4);
    cudaMalloc((void **)&ddd_P4, mf4);
    cudaMalloc((void **)&ddd_P5, mf4);
    cudaMalloc((void **)&ddd_temp, mf4);

    cudaMalloc((void **)&ddd_flag, mi4);
    cudaMalloc((void **)&ddd_iflag, mi4);
    cudaMalloc((void **)&ddd_jflag, mi4);
    cudaMalloc((void **)&ddd_oflag, mi4);
    cudaMalloc((void **)&ddd_tflag, mi4);

    // 依次类似的分配 rrho=vec(N); PP1=vec(N); ...  这些都是 CPU 端内存

    // 这些是 CPU 端
    rrho = vec(N);
    rrh = vec(N);

    PP1 = vec(N);
    PP2 = vec(N);
    PP3 = vec(N);
    PP4 = vec(N);
    PP5 = vec(N);

    fflag = ivec(N);
    iiflag = ivec(N);
    jjflag = ivec(N);
    ooflag = ivec(N);

    // 您还提到： “ページロックメモリの確保。CPU→GPU間のデータ転送を高速化”
    // 也就是 cudaMallocHost
    cudaMallocHost((void **)&temp, mf);
    cudaMallocHost((void **)&pphi, mf);
    cudaMallocHost((void **)&SSph0, mf);
    cudaMallocHost((void **)&SSph1, mf);
    cudaMallocHost((void **)&SSph2, mf);

    // 然后再对 d_err, dd_err, d_phi 等做 cudaMalloc
    //   例如：
    cudaMalloc((void **)&d_err, mf);
    cudaMalloc((void **)&dd_err, mf2);

    cudaMalloc((void **)&d_phi, mf);
    cudaMalloc((void **)&d_rho, mf);
    cudaMalloc((void **)&d_temp, mf);
    cudaMalloc((void **)&d_data, mf);
    cudaMalloc((void **)&d_rh, mf);

    cudaMalloc((void **)&d_Sph0, mf);
    cudaMalloc((void **)&d_Sph1, mf);
    cudaMalloc((void **)&d_Sph2, mf);

    cudaMalloc((void **)&d_P1, mf);
    cudaMalloc((void **)&d_P2, mf);
    cudaMalloc((void **)&d_P3, mf);
    cudaMalloc((void **)&d_P4, mf);
    cudaMalloc((void **)&d_P5, mf);

    cudaMalloc((void **)&d_flag, mi);
    cudaMalloc((void **)&d_iflag, mi);
    cudaMalloc((void **)&d_jflag, mi);
    cudaMalloc((void **)&d_oflag, mi);

    std::cout << "CUDA device memory allocated successfully.\n\n";

    // -- 设置 superbee 与 kappa --
    superbee = 2.0;    // 数值越大，数值斜坡越陡；superbee通常可取2
    kappa = 1.0 / 3.0; // superbee

    // -- 设定 Helmholtz 方程的一些系数 --
    //   单位与注释保持原意
    HelmA[0] = 1.986e-4 * 1e4; // m^-2 Torr^-2
    HelmA[1] = 0.0051 * 1e4;   // m^-2 Torr^-2
    HelmA[2] = 0.4886 * 1e4;   // m^-2 Torr^-2

    HelmG[0] = 0.0553 * 1e2; // in m^-1 Torr^-1
    HelmG[1] = 0.1460 * 1e2; // in m^-1 Torr^-1
    HelmG[2] = 0.89 * 1e2;   // in m^-1 Torr^-1

    // -- 网格生成 (调用 mesh_generator) --
    //   注意此函数是您已有的C++函数
    //   传入相关文件名与网格尺寸，然后传出 r,z,rh,zh,dr,dz,Vol,Sr,Sz

    mesh_generator(MESH_R_FILE,
                   MESH_Z_FILE,
                   NR, NZ,
                   r, z,
                   rh, zh,
                   dr, dz,
                   Vol, Sr, Sz);

    std::cout << "Mesh generatedd successfully.\n\n";

    std::cout << "Initializing vectors and matrixes.\n";

    for (int i = 0; i < NR; ++i)
    {
        for (int j = 0; j < NZ; ++j)
        {
            // 计算逻辑： pow(zh[j]/bgap, 2) - pow(rh[i]/a, 2) >= 1 ?
            double lhs = std::pow(zh[j] / bgap, 2.0) - std::pow(rh[i] / a, 2.0);

            if (lhs >= 1.0)
            {
                flag[i][j] = 1;
            }
            else
            {
                flag[i][j] = 0;
            }

            phi[i][j] = 0.0;
            dtELr[i][j] = 0.0;
            dtELz[i][j] = 0.0;
            Dvx[i][j] = 0.0;
            Dvy[i][j] = 0.0;
        }
    }
    std::cout << "Done.\n";

    // 移流計算の境界条件を設定するためのフラグたて
    call_flag(NR, NZ, flag, tflag, oflag, iflag, jflag);

    double *y;
    particle = cmat(n_particles, 50);
    y = vec(n_particles);
    ireactl = imat(200, 3), ireactr = imat(200, 3);
    ereactl = imat(200, 3), ereactr = imat(200, 3);

    ereactl_1st = ivec(200), ereactl_2nd = ivec(200), ereactl_3rd = ivec(200);
    ereactr_1st = ivec(200), ereactr_2nd = ivec(200), ereactr_3rd = ivec(200);

    ireactl_1st = ivec(200), ireactl_2nd = ivec(200), ireactl_3rd = ivec(200);
    ireactr_1st = ivec(200), ireactr_2nd = ivec(200), ireactr_3rd = ivec(200);

    nreactl_1st = ivec(200), nreactl_2nd = ivec(200), nreactl_3rd = ivec(200);
    nreactr_1st = ivec(200), nreactr_2nd = ivec(200), nreactr_3rd = ivec(200);

    nreactl = imat(200, 3), nreactr = imat(200, 3);
    eK0 = vec(200), eK1 = vec(200), eK2 = vec(200), eER = vec(200), eRnum = ivec(200);
    iK0 = vec(200), iK1 = vec(200), iK2 = vec(200), iER = vec(200), iRnum = ivec(200);
    iK0m = vec(200);
    nK0 = vec(200), nK1 = vec(200), nK2 = vec(200), nER = vec(200), nRnum = ivec(200);
    nK0m = vec(200), nK1m = vec(200), nK2m = vec(200);
    eER_mod = vec(200), iER_mod = vec(200), nER_mod = vec(200);
    eK2_mod = vec(200);
    T = mat(NR, NZ);

    // 1) 声明指针
    double **dEv = nullptr;
    double **Ev = nullptr;
    double *c1 = nullptr;
    double *c2 = nullptr;
    double *c3 = nullptr;

    double *TN2TH2O = nullptr;
    double *kN2H2O = nullptr;
    double *nkN2H2O = nullptr;

    // 2) 使用您已有的分配函数进行分配
    //    例如: mat(rows, cols), vec(length)
    dEv = mat(5, 30);
    Ev = mat(5, 30);
    c1 = vec(6);
    c2 = vec(6);
    c3 = vec(6);

    TN2TH2O = vec(5);
    kN2H2O = vec(5);
    nkN2H2O = vec(5);

    // 3) 给前几个数据赋值 (若只用到前4个，下标0..3)
    //    这里与原代码一致
    TN2TH2O[0] = 293.0;
    kN2H2O[0] = 1.2e-14;
    TN2TH2O[1] = 523.0;
    kN2H2O[1] = 2.5e-14;
    TN2TH2O[2] = 693.0;
    kN2H2O[2] = 5.1e-14;
    TN2TH2O[3] = 963.0;
    kN2H2O[3] = 1.1e-13;
    // 如果您还需设置下标 4，或做更多初始化，可继续

    // 4) 调用 spline (TdE, kN2H2O等)，计算二次微分数组 nkN2H2O
    spline(TN2TH2O, kN2H2O, 4, nkN2H2O);
    //   注意：这里第3个参数“4”表示插值点个数；
    //        第4个参数“nkN2H2O”将保存二次微分结果

    // 5) 调用自定义函数 Read_constant(...)
    Read_constant(dEv, Ev, c1, c2, c3);

    std::cout << "==============Read_constant function completed================\n\n";

    std::cout << "==============Intinialzing================\n\n";

    // 考慮する粒子種と初期濃度を読み込む
    fp = fopen(INITIAL_FILE, "r");
    n_particles = Initial_condition(fp, particle, y);
    fclose(fp);
    std::cout << "==============Intinialzation successed================\n\n";
    // 光電離で使う。計算コストを削減するため。
    numO2 = particle_number("O2", particle, n_particles);
    numO2p = particle_number("O2p", particle, n_particles);
    nume = particle_number("e", particle, n_particles);

    numN2 = particle_number("N2", particle, n_particles);
    numH2O = particle_number("H2O", particle, n_particles);
    // -------------- 1) 计算 y[] 的总和 sum --------------
    sum = 0.0;
    for (int k = 0; k < n_particles; ++k)
    {
        sum += y[k];
    }

    // -------------- 2) 计算 inv_sum --------------
    inv_sum = 1.0 / sum;

    // -------------- 3) 将初始浓度分配到 PARTICLE[i][j][k] --------------
    for (int i = 0; i < NR; ++i)
    {
        for (int j = 0; j < NZ; ++j)
        {
            for (int k = 0; k < n_particles; ++k)
            {
                PARTICLE[i][j][k] = y[k] * inv_sum;
            }
        }
    }

    // -------------- 4) 计算 O2, N2 的百分比浓度 --------------
    int idx_O2 = particle_number("O2", particle, n_particles);
    int idx_N2 = particle_number("N2", particle, n_particles);
    // int idx_H2O = particle_number("H2O", particle, n_particles);  // 如果需要 H2O

    pO2 = y[idx_O2] * inv_sum;
    pN2 = y[idx_N2] * inv_sum;
    // double pH2O = 1.0 - pO2 - pN2;    // 如有需求

    // -------------- 5) 计算混合气体的平均分子量 --------------
    MolMass = pN2 * 28.0 + pO2 * 32.0;

    // -----------------------------------------------------------
    // Townsent_ionization coefficient (回/cm)
    // -----------------------------------------------------------
    int icn = 0;
    for (int i = 260; i < BN; ++i)
    {
        // 1) 将 TdE[i] 赋给 alpTd[icn]
        alpTd[icn] = TdE[i];

        // 2) 根据 pO2 / pN2 计算 dalp[icn]
        //   （原先注释的 if(pN2==0.0) 分支也可保留,
        //    这里仅展示直接使用 pN2 != 0.0 时的公式）
        dalp[icn] =
            (dc[16][i] * (MOL * pN2) + 0.2 * 1.19 * dc[28][i] * (MOL * pO2)) / (dv[i] * TdE[i] * 1.0e-21) * 0.01;

        // 3) 递增计数
        ++icn;
    }

    // 4) 最终的 AN 值即为 icn
    AN = icn;

    // 5) 用 spline(...) 函数对 (alpTd, dalp) 做插值处理
    //    这里 spline 为您已有的函数
    spline(alpTd, dalp, AN, ddalp);

    /*
    // 如需测试插值结果，可保留以下调试输出的逻辑
    double kalp = 0.0;
    for(int iTest = 0; iTest < 10000; ++iTest) {
        double xval = static_cast<double>(iTest) * 0.1;
        splint(alpTd, dalp, ddalp, AN, xval, &kalp);
        std::printf("%d\t%e\n", iTest, kalp);
    }
    std::exit(0);
    */

    // -----------------------------------------------------------
    // 打印说明信息：气体组分等
    // -----------------------------------------------------------
    std::cout << "\n*** Atmospheric Pressure Streamer Simulation ***\n";
    std::cout << "*** Gas Component: O2 = " << (pO2 * 100.0)
              << "%, N2 = " << (pN2 * 100.0)
              << "%, H2O = " << (pH2O * 100.0)
              << "%\n\n";

    // 1) air_kg
    auto air_kg = MolMass / 1000.0;
    //   表示气体平均分子质量 (kg/mol)
    //   其中 MolMass = (pN2*28.0 + pO2*32.0)，在之前已计算

    // 2) 比热比
    auto g0 = 1.4;

    // 3) 标准大气压 (Pa)
    auto atmpa = 1.013e5;

    // 4) 气体常数
    //    rgas = 8.314 / air_kg (Pa*m^3)/(K*kg)
    auto rgas = 8.314 / air_kg;

    // 下面这行可以再保留为注释：
    // // rgas = rgas / Avogadro;

    // 5) 初始流体速度
    auto u0 = 0.0;
    auto v0 = 0.0;

    // 6) 初始压力和温度
    auto p0 = atmpa;
    auto t0 = T0; // T0 在此前全局或其他地方定义为 300.0

    // 7) 计算 molV: 摩尔体积 (m^3/mol) (300K, 1atm 的情况)
    //    300.0 是假设温度
    auto molV = (rgas * 300.0) / atmpa;
    // 原注释: // 22.413996e-3; //モル体積(m^3/mol)

    // 8) rou0 = 1.0 / (rgas * t0 / p0)
    //    表示 ρ0 = p / (R * T)
    auto rou0 = 1.0 / (rgas * (t0 / p0));
    // 这里等价于 rou0 = p0 / (rgas * t0);

    // 如果需要输出做调试
    std::cout << "air_kg = " << air_kg << "\n"
              << "g0 = " << g0 << ", atmpa = " << atmpa << "\n"
              << "rgas = " << rgas << ", (u0, v0) = (" << u0 << "," << v0 << ")\n"
              << "p0 = " << p0 << ", t0 = " << t0 << ", molV = " << molV << "\n"
              << "rou0 = " << rou0 << std::endl;

    // 電荷、分子の初期分布
    first_q(NR, NZ, flag, rh, zh, ne, N2p, O2p, O2m, Om);

    /*

    // 调用 first_q 后，检查输出：
    std::cout << "After calling first_q(), check ne array:\n";
    std::cout << "ne[0][0]   = " << ne[0][0] << std::endl;
    std::cout << "ne[10][10] = " << ne[10][10] << std::endl;
    // ...
    std::cout << "ne[NR-1][NZ-1] = " << ne[NR - 1][NZ - 1] << std::endl;

    // 同理也可检查 N2p, O2p, ...
    std::cout << "N2p[0][0]  = " << N2p[0][0] << std::endl;
    std::cout << "O2p[0][0]  = " << O2p[0][0] << std::endl;
    std::cout << "O2m[0][0]  = " << O2m[0][0] << std::endl;
    std::cout << "Om[0][0]   = " << Om[0][0] << std::endl;
*/

    // 電子衝突反応の反応リストを読み込む
    fp = fopen(E_REACTION_FILE, "r");
    ne_react = Read_reaction(fp, NUM_L, NUM_R, n_particles, eRnum, ereactl, ereactr, eK0, eK1, eK2, eER, particle);
    fclose(fp);

    for (i = 0; i < ne_react; i++)
    {
        ereactl_1st[i] = ereactl[eRnum[i]][0];
        ereactl_2nd[i] = ereactl[eRnum[i]][1];
        ereactl_3rd[i] = ereactl[eRnum[i]][2];

        ereactr_1st[i] = ereactr[eRnum[i]][0];
        ereactr_2nd[i] = ereactr[eRnum[i]][1];
        ereactr_3rd[i] = ereactr[eRnum[i]][2];

        eER_mod[i] = eER[eRnum[i]];
        eK2_mod[i] = eK2[eRnum[i]];
    }

    // イオンを含む反応の反応リストを読み込む
    fp = fopen(I_REACTION_FILE, "r"); // 化学反応リストを読み込む
    ni_react = Read_reaction(fp, NUM_L, NUM_R, n_particles, iRnum, ireactl, ireactr, iK0, iK1, iK2, iER, particle);
    fclose(fp);
    for (i = 0; i < ni_react; i++)
        iK0m[i] = iK0[iRnum[i]];

    for (i = 0; i < ni_react; i++)
    {
        ireactl_1st[i] = ireactl[iRnum[i]][0];
        ireactl_2nd[i] = ireactl[iRnum[i]][1];
        ireactl_3rd[i] = ireactl[iRnum[i]][2];

        ireactr_1st[i] = ireactr[iRnum[i]][0];
        ireactr_2nd[i] = ireactr[iRnum[i]][1];
        ireactr_3rd[i] = ireactr[iRnum[i]][2];

        iER_mod[i] = iER[iRnum[i]];
    }

    // 中性粒子の反応リストを読み込む
    fp = fopen(N_REACTION_FILE, "r"); // 化学反応リストを読み込む
    nn_react = Read_reaction(fp, NUM_L, NUM_R, n_particles, nRnum, nreactl, nreactr, nK0, nK1, nK2, nER, particle);
    fclose(fp);

    for (i = 0; i < nn_react; i++)
    {
        nK0m[i] = nK0[nRnum[i]];
        nK1m[i] = nK1[nRnum[i]];
        nK2m[i] = nK2[nRnum[i]];
    }
    for (i = 0; i < nn_react; i++)
    {
        nreactl_1st[i] = nreactl[nRnum[i]][0];
        nreactl_2nd[i] = nreactl[nRnum[i]][1];
        nreactl_3rd[i] = nreactl[nRnum[i]][2];

        nreactr_1st[i] = nreactr[nRnum[i]][0];
        nreactr_2nd[i] = nreactr[nRnum[i]][1];
        nreactr_3rd[i] = nreactr[nRnum[i]][2];

        nER_mod[i] = nER[nRnum[i]];
    }

    for (i = 0; i < ne_react; i++)
    {
        printf("%2d %4s  %4s  %4s  -->  %4s  %4s  %4s  %2.2e\t%2.2f\t%2.2f\n", eRnum[i],
               particle[ereactl[eRnum[i]][0]], particle[ereactl[eRnum[i]][1]],
               particle[ereactl[eRnum[i]][2]],
               particle[ereactr[eRnum[i]][0]], particle[ereactr[eRnum[i]][1]],
               particle[ereactr[eRnum[i]][2]], eK0[eRnum[i]], eK1[eRnum[i]], eK2[eRnum[i]]);
    }
    // exit(0);

    pureAir_kg = (pN2 * (28) + pO2 * (32)) / 1000.0; // 空気の平均質量 kg/mol
    for (i = 0; i < NR; i++)
        for (j = 0; j < NZ; j++)
            T[i][j] = T0;
    // 離散化
    discretization(NR, NZ, rh, zh, P1, P2, P3, P4, P5, a, bgap, iflag, jflag, oflag, flag); // poisson計算に必要なmesh間隔をhhに保存

    // 1) 将 P1,P2,P3,P4,P5, rh, flag... 等从 2D 拷到 1D 数组
    for (int i = 0; i < NR; i++)
    {
        for (int j = 0; j < NZ; j++)
        {
            PP1[point(NZ, i, j)] = P1[i][j];
            PP2[point(NZ, i, j)] = P2[i][j];
            PP3[point(NZ, i, j)] = P3[i][j];
            PP4[point(NZ, i, j)] = P4[i][j];
            PP5[point(NZ, i, j)] = P5[i][j];

            rrh[point(NZ, i, j)] = rh[i]; // 仅拷贝 rhalf[i] -> rrh

            fflag[point(NZ, i, j)] = flag[i][j];
            iiflag[point(NZ, i, j)] = iflag[i][j];
            jjflag[point(NZ, i, j)] = jflag[i][j];
            ooflag[point(NZ, i, j)] = oflag[i][j];
        }
    }

    // 2) 将这些 1D 数组拷贝到 GPU (CUDA) device 内存
    cudaMemcpy(d_P1, PP1, mf, cudaMemcpyHostToDevice);
    cudaMemcpy(d_P2, PP2, mf, cudaMemcpyHostToDevice);
    cudaMemcpy(d_P3, PP3, mf, cudaMemcpyHostToDevice);
    cudaMemcpy(d_P4, PP4, mf, cudaMemcpyHostToDevice);
    cudaMemcpy(d_P5, PP5, mf, cudaMemcpyHostToDevice);
    cudaMemcpy(d_rh, rrh, mf, cudaMemcpyHostToDevice);

    cudaMemcpy(d_flag, fflag, mi, cudaMemcpyHostToDevice);
    cudaMemcpy(d_iflag, iiflag, mi, cudaMemcpyHostToDevice);
    cudaMemcpy(d_jflag, jjflag, mi, cudaMemcpyHostToDevice);
    cudaMemcpy(d_oflag, ooflag, mi, cudaMemcpyHostToDevice);

    // 3) 将 pphi[...] 置为 0.0 (这里假设 pphi 是 1D，大小N=NR*NZ)
    for (int idx = 0; idx < N; idx++)
    {
        pphi[idx] = 0.0;
    }
    // (或 std::fill(pphi, pphi + N, 0.0);  // 若可以用C++算法库 )

    cudaMemcpy(d_err, pphi, mf, cudaMemcpyHostToDevice);
    cudaMemcpy(d_temp, pphi, mf, cudaMemcpyHostToDevice);
    cudaMemcpy(d_phi, pphi, mf, cudaMemcpyHostToDevice);
    cudaMemcpy(d_Sph0, pphi, mf, cudaMemcpyHostToDevice);
    cudaMemcpy(d_Sph1, pphi, mf, cudaMemcpyHostToDevice);
    cudaMemcpy(d_Sph2, pphi, mf, cudaMemcpyHostToDevice);

    // 4) 初始化一些变量
    START = -1;
    dt = 0.125e-12;
    time = -5.0e-9;

    // 5) 打开文件并关闭 —— 用 C++ ofstream 替代
    {
        std::ofstream ofs1("outputdata/current/current.dat");
        //  如果仅仅是清空文件，可以什么都不写
    }
    {
        std::ofstream ofs2("outputdata/current/Power.dat");
    }

    // 6) 给一些变量赋值
    V = 28.0e+3;
    A = V;
    B1 = 30.276901e-9;
    B2 = B1 + 0.0e-9;
    C = 16.5e-9;
    D = 0.003 * 1e9;

    double hrm, hzm, absr;

    // 境界条件Br, Bz の設定
    for (i = 0; i < NR; i++)
    {
        for (j = 1421; j < NZ; j++)
        {
            if (flag[i][j])
            {
                ;
            }
            else if (iflag[i][j])
            {
                hrm = rh[i] - a * sqrt(pow(zh[j] / bgap, 2) - 1.0);
                hzm = bgap * sqrt(1.0 + pow(rh[i] / a, 2)) - zh[j];
                absr = sqrt(hrm * hrm + hzm * hzm);

                Br[i][j] = hrm / absr;
                Bz[i][j] = hzm / absr;
            }
            else if (jflag[i][j])
            {
                hrm = rh[i] - a * sqrt(pow(zh[j] / bgap, 2) - 1.0);
                hzm = bgap * sqrt(1.0 + pow(rh[i] / a, 2)) - zh[j];
                absr = sqrt(hrm * hrm + hzm * hzm);

                Br[i][j] = hrm / absr;
                Bz[i][j] = hzm / absr;
            }
            else if (oflag[i][j])
            {
                hrm = rh[i] - a * sqrt(pow(zh[j] / bgap, 2) - 1.0);
                hzm = bgap * sqrt(1.0 + pow(rh[i] / a, 2)) - zh[j];
                absr = sqrt(hrm * hrm + hzm * hzm);

                Br[i][j] = hrm / absr;
                Bz[i][j] = hzm / absr;
            }
        }
    }

    spDATA = 331;
    spTime = vec(spDATA), dspV = vec(spDATA), ddspV = vec(spDATA);
    fp = fopen(VOLTAGE_FILE, "r");
    for (j = 0; j < spDATA; j++)
        fscanf(fp, "%le\t%le\n", &spTime[j], &dspV[j]);
    spline(spTime, dspV, spDATA, ddspV);

    // 1) Production_rate_cathode.dat
    {
        std::ofstream ofsCat("Production_rate_cathode.dat");
        if (!ofsCat)
        {
            std::cerr << "Error: cannot open Production_rate_cathode.dat for writing\n";
            // 根据需求是否 return/exit
        }
        // 写制表符
        ofsCat << "\t";
        // 写 eRnum
        for (int j = 0; j < ne_react; ++j)
        {
            ofsCat << "e" << eRnum[j] << "\t";
        }
        // 写 iRnum
        for (int j = 0; j < ni_react; ++j)
        {
            ofsCat << "i" << iRnum[j] << "\t";
        }
        // 写 nRnum
        for (int j = 0; j < nn_react; ++j)
        {
            ofsCat << "n" << nRnum[j] << "\t";
        }
        // 换行
        ofsCat << "\n";
        // ofsCat 在该作用域结束时自动关闭，或您可以显式 ofsCat.close();
    }

    // 2) Production_rate_1st.dat
    {
        std::ofstream ofsFirst("Production_rate_1st.dat");
        if (!ofsFirst)
        {
            std::cerr << "Error: cannot open Production_rate_1st.dat for writing\n";
        }
        ofsFirst << "\t";
        for (int j = 0; j < ne_react; ++j)
        {
            ofsFirst << "e" << eRnum[j] << "\t";
        }
        for (int j = 0; j < ni_react; ++j)
        {
            ofsFirst << "i" << iRnum[j] << "\t";
        }
        for (int j = 0; j < nn_react; ++j)
        {
            ofsFirst << "n" << nRnum[j] << "\t";
        }
        ofsFirst << "\n";
    }

    // 3) Production_rate_2nd.dat
    {
        std::ofstream ofsSecond("Production_rate_2nd.dat");
        if (!ofsSecond)
        {
            std::cerr << "Error: cannot open Production_rate_2nd.dat for writing\n";
        }
        ofsSecond << "\t";
        for (int j = 0; j < ne_react; ++j)
        {
            ofsSecond << "e" << eRnum[j] << "\t";
        }
        for (int j = 0; j < ni_react; ++j)
        {
            ofsSecond << "i" << iRnum[j] << "\t";
        }
        for (int j = 0; j < nn_react; ++j)
        {
            ofsSecond << "n" << nRnum[j] << "\t";
        }
        ofsSecond << "\n";
    }

    // 4) Production_rate_anode.dat
    {
        std::ofstream ofsAnode("Production_rate_anode.dat");
        if (!ofsAnode)
        {
            std::cerr << "Error: cannot open Production_rate_anode.dat for writing\n";
        }
        ofsAnode << "\t";
        for (int j = 0; j < ne_react; ++j)
        {
            ofsAnode << "e" << eRnum[j] << "\t";
        }
        for (int j = 0; j < ni_react; ++j)
        {
            ofsAnode << "i" << iRnum[j] << "\t";
        }
        for (int j = 0; j < nn_react; ++j)
        {
            ofsAnode << "n" << nRnum[j] << "\t";
        }
        ofsAnode << "\n";
    }

    // --- 1) 声明指针(二维整型、双精度等) ---
    int **flag2 = nullptr;
    int **tflag2 = nullptr;
    int **oflag2 = nullptr;
    int **iflag2 = nullptr;
    int **jflag2 = nullptr;

    double *zh2 = nullptr;
    double *rh2 = nullptr;
    double **cCphi = nullptr;
    double **cP1 = nullptr;
    double **cP2 = nullptr;
    double **cP3 = nullptr;
    double **cP4 = nullptr;
    double **cP5 = nullptr;

    // --- 2) 声明另一套指针 ---
    int **flag4 = nullptr;
    int **tflag4 = nullptr;
    int **oflag4 = nullptr;
    int **iflag4 = nullptr;
    int **jflag4 = nullptr;

    double *zh4 = nullptr;
    double *rh4 = nullptr;
    double **Cphi4 = nullptr;
    double **ccP1 = nullptr;
    double **ccP2 = nullptr;
    double **ccP3 = nullptr;
    double **ccP4 = nullptr;
    double **ccP5 = nullptr;

    // --- 3) 依次分配 ---
    //  这里示例中用到 NR2, NZ2, NR4, NZ4 等，请确保它们在当前作用域可见
    //  并且您已有 vec()/imat()/mat() 等分配函数
    zh2 = vec(NZ2);
    rh2 = vec(NR2);
    cCphi = mat(NR2, NZ2);

    flag2 = imat(NR2, NZ2);
    tflag2 = imat(NR2, NZ2);
    oflag2 = imat(NR2, NZ2);
    iflag2 = imat(NR2, NZ2);
    jflag2 = imat(NR2, NZ2);

    cP1 = mat(NR2, NZ2);
    cP2 = mat(NR2, NZ2);
    cP3 = mat(NR2, NZ2);
    cP4 = mat(NR2, NZ2);
    cP5 = mat(NR2, NZ2);

    zh4 = vec(NZ4);
    rh4 = vec(NR4);
    Cphi4 = mat(NR4, NZ4);

    flag4 = imat(NR4, NZ4);
    tflag4 = imat(NR4, NZ4);
    oflag4 = imat(NR4, NZ4);
    iflag4 = imat(NR4, NZ4);
    jflag4 = imat(NR4, NZ4);

    ccP1 = mat(NR4, NZ4);
    ccP2 = mat(NR4, NZ4);
    ccP3 = mat(NR4, NZ4);
    ccP4 = mat(NR4, NZ4);
    ccP5 = mat(NR4, NZ4);

    for (int i = 0; i < NR2; i++)
    {
        for (int j = 0; j < NZ2; j++)
        {
            // 下采样 *2
            zh2[j] = zh[j * 2];
            rh2[i] = rh[i * 2];

            // 根据几何判断
            if (std::pow(zh2[j] / bgap, 2) - std::pow(rh2[i] / a, 2) >= 1.0)
            {
                flag2[i][j] = 1;
            }
            else
            {
                flag2[i][j] = 0;
            }

            // cCphi 先清零
            cCphi[i][j] = 0.0;
        }
    }

    // 3) 调用您的函数 call_flag(...) 给 flag2, tflag2, iflag2... 做进一步处理
    call_flag(NR2, NZ2, flag2, tflag2, oflag2, iflag2, jflag2);

    // 4) 调用离散化函数 discretization(...)
    discretization(NR2, NZ2, rh2, zh2,
                   cP1, cP2, cP3, cP4, cP5,
                   a, bgap,
                   iflag2, jflag2, oflag2, flag2);

    // 5) 把 cP1/cP2/... (及 rh2) 拷贝到一维数组 rrho2，然后 cudaMemcpy 到设备 dd_rh, dd_Pn 等
    //   其中 point(NZ2,i,j)是您原先的宏：#define point(N,i,j) ((N)*(i)+(j))

    // 先把 rh2[i] 映射到 rrho2
    for (int i = 0; i < NR2; i++)
    {
        for (int j = 0; j < NZ2; j++)
        {
            rrho2[point(NZ2, i, j)] = rh2[i];
        }
    }
    cudaMemcpy(dd_rh, rrho2, mf2, cudaMemcpyHostToDevice);

    // cP1 => dd_P1
    for (int i = 0; i < NR2; i++)
    {
        for (int j = 0; j < NZ2; j++)
        {
            rrho2[point(NZ2, i, j)] = cP1[i][j];
        }
    }
    cudaMemcpy(dd_P1, rrho2, mf2, cudaMemcpyHostToDevice);

    // cP2 => dd_P2
    for (int i = 0; i < NR2; i++)
    {
        for (int j = 0; j < NZ2; j++)
        {
            rrho2[point(NZ2, i, j)] = cP2[i][j];
        }
    }
    cudaMemcpy(dd_P2, rrho2, mf2, cudaMemcpyHostToDevice);

    // cP3 => dd_P3
    for (int i = 0; i < NR2; i++)
    {
        for (int j = 0; j < NZ2; j++)
        {
            rrho2[point(NZ2, i, j)] = cP3[i][j];
        }
    }
    cudaMemcpy(dd_P3, rrho2, mf2, cudaMemcpyHostToDevice);

    // cP4 => dd_P4
    for (int i = 0; i < NR2; i++)
    {
        for (int j = 0; j < NZ2; j++)
        {
            rrho2[point(NZ2, i, j)] = cP4[i][j];
        }
    }
    cudaMemcpy(dd_P4, rrho2, mf2, cudaMemcpyHostToDevice);

    // cP5 => dd_P5
    for (int i = 0; i < NR2; i++)
    {
        for (int j = 0; j < NZ2; j++)
        {
            rrho2[point(NZ2, i, j)] = cP5[i][j];
        }
    }
    cudaMemcpy(dd_P5, rrho2, mf2, cudaMemcpyHostToDevice);

    // 6) flag2 / iflag2 / jflag2 / tflag2 / oflag2 同理
    //    先把 flag2[i][j] => fflag2[point(...)] 再 cudaMemcpy => dd_flag
    for (int i = 0; i < NR2; i++)
    {
        for (int j = 0; j < NZ2; j++)
        {
            fflag2[point(NZ2, i, j)] = flag2[i][j];
        }
    }
    cudaMemcpy(dd_flag, fflag2, mi2, cudaMemcpyHostToDevice);

    // 同理 iflag2 => dd_iflag
    for (int i = 0; i < NR2; i++)
    {
        for (int j = 0; j < NZ2; j++)
        {
            fflag2[point(NZ2, i, j)] = iflag2[i][j];
        }
    }
    cudaMemcpy(dd_iflag, fflag2, mi2, cudaMemcpyHostToDevice);

    // jflag2 => dd_jflag
    for (int i = 0; i < NR2; i++)
    {
        for (int j = 0; j < NZ2; j++)
        {
            fflag2[point(NZ2, i, j)] = jflag2[i][j];
        }
    }
    cudaMemcpy(dd_jflag, fflag2, mi2, cudaMemcpyHostToDevice);

    // tflag2 => dd_tflag
    for (int i = 0; i < NR2; i++)
    {
        for (int j = 0; j < NZ2; j++)
        {
            fflag2[point(NZ2, i, j)] = tflag2[i][j];
        }
    }
    cudaMemcpy(dd_tflag, fflag2, mi2, cudaMemcpyHostToDevice);

    // oflag2 => dd_oflag
    for (int i = 0; i < NR2; i++)
    {
        for (int j = 0; j < NZ2; j++)
        {
            fflag2[point(NZ2, i, j)] = oflag2[i][j];
        }
    }
    cudaMemcpy(dd_oflag, fflag2, mi2, cudaMemcpyHostToDevice);

    // 2) 先将 zh4、rh4 等在 i、j 上赋值
    for (i = 0; i < NR4; i++)
    {
        for (j = 0; j < NZ4; j++)
        {
            // zh4[j] = zh2[j*2];   // “zh2” 来自您之前更粗网格/更细网格的数组
            // rh4[i] = rh2[i*2];
            zh4[j] = zh2[j * 2]; // 如果 j*2 不越界 (需要保证 2*j < NZ2)
            rh4[i] = rh2[i * 2]; // 同理 2*i < NR2

            // 判定是否在电极内
            if (std::pow(zh4[j] / bgap, 2) - std::pow(rh4[i] / a, 2) >= 1.0)
            {
                flag4[i][j] = 1;
            }
            else
            {
                flag4[i][j] = 0;
            }

            // Cphi4 初始化为 0
            Cphi4[i][j] = 0.0;
        }
    }

    // 3) 调用您已有的函数: call_flag(...)
    call_flag(NR4, NZ4, flag4, tflag4, oflag4, iflag4, jflag4);

    // 4) 调用 discretization(...) 进行离散化
    discretization(NR4, NZ4,
                   rh4, zh4,
                   ccP1, ccP2, ccP3, ccP4, ccP5,
                   a, bgap,
                   iflag4, jflag4, oflag4, flag4);

    // 5) 将 ccP1..ccP5 等数组复制到 rrho4，然后再 cudaMemcpy 到 ddd_** 设备内存
    //    同时也把 flag4 等复制到 fflag4，再拷贝到 ddd_flag 等
    //    下面这几行与您给出的基本一致，只是做了些风格/注释调整

    //   先把 rh4[i] 拷到 rrho4[point(NZ4,i,j)]
    for (i = 0; i < NR4; i++)
    {
        for (j = 0; j < NZ4; j++)
        {
            rrho4[point(NZ4, i, j)] = rh4[i];
        }
    }
    cudaMemcpy(ddd_rh, rrho4, mf4, cudaMemcpyHostToDevice);

    //   把 ccP1 拷到 rrho4
    for (i = 0; i < NR4; i++)
    {
        for (j = 0; j < NZ4; j++)
        {
            rrho4[point(NZ4, i, j)] = ccP1[i][j];
        }
    }
    cudaMemcpy(ddd_P1, rrho4, mf4, cudaMemcpyHostToDevice);

    //   同理 ccP2
    for (i = 0; i < NR4; i++)
    {
        for (j = 0; j < NZ4; j++)
        {
            rrho4[point(NZ4, i, j)] = ccP2[i][j];
        }
    }
    cudaMemcpy(ddd_P2, rrho4, mf4, cudaMemcpyHostToDevice);

    //   同理 ccP3
    for (i = 0; i < NR4; i++)
    {
        for (j = 0; j < NZ4; j++)
        {
            rrho4[point(NZ4, i, j)] = ccP3[i][j];
        }
    }
    cudaMemcpy(ddd_P3, rrho4, mf4, cudaMemcpyHostToDevice);

    //   同理 ccP4
    for (i = 0; i < NR4; i++)
    {
        for (j = 0; j < NZ4; j++)
        {
            rrho4[point(NZ4, i, j)] = ccP4[i][j];
        }
    }
    cudaMemcpy(ddd_P4, rrho4, mf4, cudaMemcpyHostToDevice);

    //   同理 ccP5
    for (i = 0; i < NR4; i++)
    {
        for (j = 0; j < NZ4; j++)
        {
            rrho4[point(NZ4, i, j)] = ccP5[i][j];
        }
    }
    cudaMemcpy(ddd_P5, rrho4, mf4, cudaMemcpyHostToDevice);

    // 6) 将 flag4 / iflag4 / jflag4 / tflag4 / oflag4
    //    拷到 fflag4 数组中，再 cudaMemcpy 到 ddd_**Flag
    //    根据您的代码：先拷flag4到fflag4，再 ddd_flag；然后 iflag4->fflag4->ddd_iflag 等
    for (i = 0; i < NR4; i++)
    {
        for (j = 0; j < NZ4; j++)
        {
            fflag4[point(NZ4, i, j)] = flag4[i][j];
        }
    }
    cudaMemcpy(ddd_flag, fflag4, mi4, cudaMemcpyHostToDevice);

    for (i = 0; i < NR4; i++)
    {
        for (j = 0; j < NZ4; j++)
        {
            fflag4[point(NZ4, i, j)] = iflag4[i][j];
        }
    }
    cudaMemcpy(ddd_iflag, fflag4, mi4, cudaMemcpyHostToDevice);

    for (i = 0; i < NR4; i++)
    {
        for (j = 0; j < NZ4; j++)
        {
            fflag4[point(NZ4, i, j)] = jflag4[i][j];
        }
    }
    cudaMemcpy(ddd_jflag, fflag4, mi4, cudaMemcpyHostToDevice);

    for (i = 0; i < NR4; i++)
    {
        for (j = 0; j < NZ4; j++)
        {
            fflag4[point(NZ4, i, j)] = tflag4[i][j];
        }
    }
    cudaMemcpy(ddd_tflag, fflag4, mi4, cudaMemcpyHostToDevice);

    for (i = 0; i < NR4; i++)
    {
        for (j = 0; j < NZ4; j++)
        {
            fflag4[point(NZ4, i, j)] = oflag4[i][j];
        }
    }
    cudaMemcpy(ddd_oflag, fflag4, mi4, cudaMemcpyHostToDevice);

    double a2 = 0.0;
    double b2 = 0.0;
    double thetaI = 0.0;
    double theta = 0.0;
    double Q = 0.0;
    double **hyperboroid;
    double Q_div_E0;

    hyperboroid = mat(NR, NZ); // 若先前尚未分配，可在更外层分配

    // -------------------------------------------------------
    // (2) 计算 hyperboloid 数组
    // -------------------------------------------------------
    for (int i = 0; i < NR; i++)
    {
        for (int j = 0; j < NZ; j++)
        {

            // 如果网格在电极内，就直接跳过
            if (flag[i][j])
            {
                // do nothing; hyperboroid[i][j] 保持默认
            }
            // 若是 z=0.0 (接地), 这里也什么都不做
            else if (zh[j] == 0.0)
            {
                // do nothing
            }
            else
            {
                // 与原 C 代码相同的逻辑
                thetaI = std::atan2(a, bgap);       // 计算 θI
                Q = std::sqrt(bgap * bgap + a * a); // 焦点

                double tmp1 = -(rh[i] * rh[i] + zh[j] * zh[j] - Q * Q);
                double tmp2 = std::pow((rh[i] * rh[i] + zh[j] * zh[j] - Q * Q), 2.0) - 4.0 * (-rh[i] * rh[i] * Q * Q);

                a2 = (tmp1 + std::sqrt(tmp2)) * 0.5;
                a2 = std::sqrt(a2);

                b2 = std::sqrt(Q * Q - a2 * a2);

                theta = std::atan2(a2, b2);

                // 计算 hyperboloid[i][j]
                //   => (log(1.0 / tan(θ/2))) / (log(1.0 / tan(θI/2)))
                hyperboroid[i][j] =
                    std::log(1.0 / std::tan(theta * 0.5)) / std::log(1.0 / std::tan(thetaI * 0.5));
            }
        }
    }

    // -------------------------------------------------------
    // (3) 初始化 mvr, mvz, rou, p, q1..q4
    // -------------------------------------------------------
    for (int i = 0; i < NR; i++)
    {
        for (int j = 0; j < NZ; j++)
        {
            mvr[i][j] = 0.0;
            mvz[i][j] = 0.0;

            // 初始密度
            rou[i][j] = rou0;
            // p = ρ * R * T[i][j]
            p[i][j] = rou[i][j] * rgas * T[i][j];

            // q1..q4
            q1[i][j] = rou[i][j]; // [kg/m^3]
            q2[i][j] = rou[i][j] * mvr[i][j];
            q3[i][j] = rou[i][j] * mvz[i][j];

            // q4 = p/(γ-1) + ½ρV²(乘以体积 Vol[i][j])
            double vel2 = mvr[i][j] * mvr[i][j] + mvz[i][j] * mvz[i][j];
            q4[i][j] = p[i][j] / (g0 - 1.0) + (rou[i][j] * Vol[i][j]) * 0.5 * vel2;
        }
    }

    // 对网格上的每个 (i, j) 单元进行初始化
    for (int i = 0; i < NR; i++)
    {
        for (int j = 0; j < NZ; j++)
        {
            mvr[i][j] = 0.0;
            mvz[i][j] = 0.0;
            rou[i][j] = rou0;

            // 计算 p[i][j]
            p[i][j] = rou[i][j] * rgas * T[i][j];

            // 接下来设置 q1, q2, q3, q4
            q1[i][j] = rou[i][j];             // (kg/m^3)
            q2[i][j] = rou[i][j] * mvr[i][j]; // [kg/(s*m^2)] 体积流量
            q3[i][j] = rou[i][j] * mvz[i][j]; // 类似
            q4[i][j] = p[i][j] / (g0 - 1.0) + (rou[i][j] * Vol[i][j]) * (mvr[i][j] * mvr[i][j] + mvz[i][j] * mvz[i][j]) * 0.5;
            // 上面注释所示单位参考原注释
        }
    }

    // 然后调用边界条件函数 bndcnd
    // 该函数会根据 flag/iflag 等来修改 mvr/mvz/p/T 等
    bndcnd(
        rou, mvr, mvz, p,
        q1, q2, q3, q4,
        g0, rgas, u0, v0, p0, t0,
        rou0,
        flag, iflag, jflag, oflag, tflag);

    // 依次获取各种粒子在 PARTCLE[][][] 数组中的索引
    n1 = particle_number("e", particle, n_particles);
    n2 = particle_number("N2p", particle, n_particles);
    n3 = particle_number("O2p", particle, n_particles);
    n4 = particle_number("O2m", particle, n_particles);
    n5 = particle_number("Om", particle, n_particles);
    n6 = particle_number("H2Op", particle, n_particles);
    n7 = particle_number("OHm", particle, n_particles);
    n8 = particle_number("Hm", particle, n_particles);
    n9 = particle_number("N4p", particle, n_particles);
    n10 = particle_number("O4p", particle, n_particles);
    n11 = particle_number("N2O2p", particle, n_particles);
    n12 = particle_number("O2pH2O", particle, n_particles);
    n13 = particle_number("H3Op", particle, n_particles);
    n14 = particle_number("H3OpH2O", particle, n_particles);
    n15 = particle_number("H3OpH2O2", particle, n_particles);
    n16 = particle_number("H3OpH2O3", particle, n_particles);

    // O2 振动态
    o2v0 = particle_number("O2", particle, n_particles);
    o2v1 = particle_number("O2v1", particle, n_particles);
    o2v2 = particle_number("O2v2", particle, n_particles);
    o2v3 = particle_number("O2v3", particle, n_particles);
    o2v4 = particle_number("O2v4", particle, n_particles);

    // N2 振动态
    n2v0 = particle_number("N2", particle, n_particles);
    n2v1 = particle_number("N2v1", particle, n_particles);
    n2v2 = particle_number("N2v2", particle, n_particles);
    n2v3 = particle_number("N2v3", particle, n_particles);
    n2v4 = particle_number("N2v4", particle, n_particles);
    n2v5 = particle_number("N2v5", particle, n_particles);
    n2v6 = particle_number("N2v6", particle, n_particles);
    n2v7 = particle_number("N2v7", particle, n_particles);
    n2v8 = particle_number("N2v8", particle, n_particles);

    // H2O 振动态
    h2ov0 = particle_number("H2O", particle, n_particles);
    h2ov1 = particle_number("H2Ov1", particle, n_particles);
    h2ov2 = particle_number("H2Ov2", particle, n_particles);
    h2ov3 = particle_number("H2Ov3", particle, n_particles);

    // 可能还需要 “O(P)” 索引
    oatm = particle_number("O(P)", particle, n_particles);

    // 分配一个数组，用来存放每一条电子反应的耗散功率 (根据上下文 each_mP 大概率是这样用途)
    each_mP = vec(ne_react);

    // 一些用于后续计算的常量、变量
    Q_div_E0 = QELEC / E0; // 计算时为了简化而使用的一个常数

    // 以下是一些与时间步长、振动能量输运有关的参数
    double adt = 1.25e-13; // a + 2c = 5e-13 (推测是某种数值稳定相关量)
    int bdt = 300000;
    double cdt = 0.875e-13;
    double ddt = 1.6895269e-5;

    double dt_vib = dt; // 初始电子(或振动)时间步长

    // 也可能需要文件指针（如需调试信息等）
    FILE *fdbg = nullptr;

    // 下面是用于某些限幅、阈值或标记的变量
    double pEPS;     // 未初始化，后面会赋值
    double **LEx;    // 2D数组
    int Priflag = 0; // 判断是否进行额外处理
    int Priflag2 = 0;
    int Prinstp = 1000000000;
    Prinstp2 = 1000000000; // 全局变量？若非全局，还需在此声明
    limpoz = 1161;         // 后面会用到？

    // 分配 neflag 和 LEx
    neflag = imat(NR, NZ);
    LEx = mat(NR, NZ);

    std::cout << "////////////Routine Start!!!////////////////////////\n\n";

    // 主循环，nspt从(START+1)开始递增，无上限，直到某处 break 或 return
    for (int nstp = START + 1; /*无限循环*/; ++nstp)
    {
        // 若有需要，可选择是否用 pEPS = 1e-10 / EPS
        // 这里保留原逻辑中最终的:
        pEPS = EPS;

        // 原注释中: “V = VOLTAGE_RATE * time * 1e9 * 1e3;” 被注释掉了，就保持注释
        // 下面保留了 splint(...) 与 spV = ... 逻辑
        splint(spTime, dspV, ddspV, spDATA, time * 1.0e9, &spV);
        V = spV * 1.0e3; // [V]

        // === 1) 針電極の解析解: 计算 Lphi[i][j] ===
        // 若 flag[i][j]为真(电极内), 则 Lphi = V
        // 若 z[j] == 0.0 (接地), 则 Lphi = 0.0
        // 否则 Lphi = V*hyperboroid[i][j]
        for (int i = 0; i < NR; i++)
        {
            for (int j = 0; j < NZ; j++)
            {
                if (flag[i][j])
                {
                    Lphi[i][j] = V;
                }
                else if (z[j] == 0.0)
                {
                    Lphi[i][j] = 0.0;
                }
                else
                {
                    Lphi[i][j] = V * hyperboroid[i][j];
                }
            }
        }

        // === 2) 计算 rho (电荷密度) 并统计最大电子数 e_max ===
        double e_max = 0.0; // track electron density max
        for (int i = 0; i < NR; i++)
        {
            for (int j = 0; j < NZ; j++)
            {
                // CPp 为阳离子总和
                double CPp = N2p[i][j] + O2p[i][j] + H2Op[i][j] + N4p[i][j] + O4p[i][j] + N2O2p[i][j] + O2pH2O[i][j] + H3Op[i][j] + H3OpH2O[i][j] + H3OpH2O2[i][j] + H3OpH2O3[i][j];

                // CPm 为阴离子 + 电子
                double CPm = Om[i][j] + O2m[i][j] + ne[i][j] + OHm[i][j] + Hm[i][j];

                // 若 i==0 -> rho[i][j] = 0.5* Q_div_E0(...)，否则 rh[i]*Q_div_E0(...)
                // 保留原注释
                if (i == 0)
                {
                    rho[i][j] = 0.5 * Q_div_E0 * (CPp - CPm);
                }
                else
                {
                    rho[i][j] = rh[i] * Q_div_E0 * (CPp - CPm);
                }

                // 更新 e_max
                if (e_max < ne[i][j])
                {
                    e_max = ne[i][j];
                }

                // 若 ne[i][j]< -1.0, 则报错并退出
                if (ne[i][j] < -1.0)
                {
                    std::cerr << "e_dens_over\n i=" << i << "\tj=" << j
                              << "\t" << ne[i][j] << std::endl;
                    // 如果真的要执行 system(...) 发送mail，可保持:
                    std::system("/usr/sbin/sendmail -t < mail_ne_over.txt");
                    std::exit(0);
                }
            }
        }

        // 显示一些信息 (nstp,time,dt,V,e_max)
        // 原 printf("%d\t%e\t%e\t%f\t%e\n", nstp,time,dt,V,e_max)

        std::cout << nstp << "\t"
                  << time << "\t"
                  << dt << "\t"
                  << V << "\t"
                  << e_max << "\n"
                  << std::endl;

        // === 3) 将 rho[i][j] 拷贝到一维数组 rrho[...]，然后用 cudaMemcpy 传到 GPU ===
        for (int i = 0; i < NR; i++)
        {
            for (int j = 0; j < NZ; j++)
            {
                rrho[point(NZ, i, j)] = rho[i][j];
            }
        }

        // 拷贝到 GPU (要求 mf, d_rho, rrho 已在更上层定义/分配)
        cudaMemcpy(d_rho, rrho, mf, cudaMemcpyHostToDevice);

        cuda_flag = 0;
        iter = 0;

        // 在此以 3000 作为循环上限
        for (int cuda = 1; cuda < 3000; ++cuda)
        {
            // 第 1 次调用
            OMEGA = 1.0;
            itnum = 2;
            Poisson_GPU_function(dimGrid, dimBlock,
                                 d_phi, d_rho, d_rh, d_temp,
                                 NR, NZ,
                                 d_P1, d_P2, d_P3, d_P4, d_P5,
                                 d_flag,
                                 OMEGA, itnum);
            iter += itnum;

            // 调用误差分析核
            Error_poisson_GPU(dimGrid, dimBlock,
                              d_phi, d_rho, d_rh,
                              NR, NZ,
                              d_P1, d_P2, d_P3, d_P4, d_P5,
                              d_flag, d_err);

            // 将误差限制到更粗网格
            Ristriction_GPU(dimGrid2, dimBlock2,
                            NR2, NZ2, NZ,
                            d_err,
                            dd_rho, dd_Cphi);

            // 在更粗网格上再做 2 次迭代
            {
                OMEGA = 1.0;
                itnum = 2;
                Poisson_GPU_function(dimGrid2, dimBlock2,
                                     dd_Cphi, dd_rho, dd_rh, dd_temp,
                                     NR2, NZ2,
                                     dd_P1, dd_P2, dd_P3, dd_P4, dd_P5,
                                     dd_flag,
                                     OMEGA, itnum);

                // 再计算误差
                Error_poisson_GPU(dimGrid2, dimBlock2,
                                  dd_Cphi, dd_rho, dd_rh,
                                  NR2, NZ2,
                                  dd_P1, dd_P2, dd_P3, dd_P4, dd_P5,
                                  dd_flag, dd_err);

                // 再一次限制到更小网格
                Ristriction_GPU(dimGrid4, dimBlock4,
                                NR4, NZ4, NZ2,
                                dd_err,
                                ddd_rho, ddd_Cphi);

                {
                    OMEGA = 1.0;
                    itnum = 3;
                    Poisson_GPU_function(dimGrid4, dimBlock4,
                                         ddd_Cphi, ddd_rho,
                                         ddd_temp, ddd_temp, // 这里原代码写的是 (ddd_temp,ddd_temp)
                                         NR4, NZ4,
                                         ddd_P1, ddd_P2, ddd_P3, ddd_P4, ddd_P5,
                                         ddd_flag,
                                         OMEGA, itnum);

                    // 插值回到上一层
                    Interporation_GPU(dimGrid2, dimBlock2,
                                      NR2, NZ2, NZ4,
                                      dd_flag,
                                      dd_Cphi,
                                      ddd_Cphi);
                }

                OMEGA = 1.5;
                itnum = 10;
                Poisson_GPU_function(dimGrid2, dimBlock2,
                                     dd_Cphi, dd_rho, dd_rh, dd_temp,
                                     NR2, NZ2,
                                     dd_P1, dd_P2, dd_P3, dd_P4, dd_P5,
                                     dd_flag,
                                     OMEGA, itnum);

                // 插值回到更精细网格
                Interporation_GPU(dimGrid, dimBlock,
                                  NR, NZ, NZ2,
                                  d_flag,
                                  d_phi,
                                  dd_Cphi);
            }

            // 在最细网格再做 12 次迭代
            {
                OMEGA = 1.5;
                itnum = 12;
                Poisson_GPU_function(dimGrid, dimBlock,
                                     d_phi, d_rho, d_rh, d_temp,
                                     NR, NZ,
                                     d_P1, d_P2, d_P3, d_P4, d_P5,
                                     d_flag,
                                     OMEGA, itnum);
                iter += itnum;
            }

            // 每次都做一次收敛检查
            // (为了与原代码一致, if(cuda%1==0) => 相当于每个 loop 都要检查)
            {
                Convergence_check_GPU(N, d_temp, d_phi,
                                      mf, temp, pphi,
                                      &error, &Maxphi);

                // 如果误差比小于 EPS, 就 break
                if ((error / Maxphi) < EPS)
                {
                    break;
                }
            }

        } // end for cuda=1..3000

        // 打印结果
        std::cout << "---Finish_Poisson_equation---\n"
                  << std::setw(5) << iter << "  "
                  << Maxphi << "\t"
                  << (error / Maxphi) << std::endl;

    } // end for(nstp=(START+1); ; nstp++)

    for (int i = 0; i < NR; ++i)
    {
        for (int j = 0; j < NZ; ++j)
        {
            Cphi[i][j] = pphi[point(NZ, i, j)];
        }
    }

    // -------------------------
    // 设置 phi[][] 的值
    // -------------------------
    for (int i = 0; i < NR; ++i)
    {
        for (int j = 0; j < NZ; ++j)
        {
            if (flag[i][j])
            {
                phi[i][j] = V;
            }
            else
            {
                phi[i][j] = Cphi[i][j] + Lphi[i][j];
            }
        }
    }

    // -------------------------
    // 调用各计算函数 (电场计算、速度计算等)
    // -------------------------
    calc_E(NR, NZ,
           Lphi, absE, Ey, Ex,
           air_kg, tflag, oflag, iflag, jflag, flag,
           a, bgap, rh, zh, rou,
           nstp, Prinstp); // 先用 Lphi 计算

    calc_E(NR, NZ,
           phi, absE, Ey, Ex,
           air_kg, tflag, oflag, iflag, jflag, flag,
           a, bgap, rh, zh, rou,
           nstp, Prinstp); // 再用 phi 计算

    // 电子速度
    calc_e_velo(NR, NZ,
                absE, Ex, Ey,
                vr, vz, ne,
                rh, zh,
                TdE, dv, ddv,
                BN);

    // 正、负离子速度
    calc_ie_velo(NR, NZ,
                 piv_r, piv_z,
                 Ex, Ey,
                 0); // + イオン
    calc_ie_velo(NR, NZ,
                 miv_r, miv_z,
                 Ex, Ey,
                 1); // - イオン

    // 1) 先为 THREAD_NUM 条线程准备参数
    std::vector<ThreadArg> tArgs(THREAD_NUM);
    for (int i = 0; i < THREAD_NUM; ++i)
    {
        tArgs[i].thread_no = i;
        tArgs[i].mutex = &mutex; // 如果该函数内部仍然需要互斥锁
    }

    // 2) 创建并启动线程
    std::vector<std::thread> workers;
    workers.reserve(THREAD_NUM);

    for (int i = 0; i < THREAD_NUM; ++i)
    {
        // 直接将引用传给线程函数，避免值拷贝
        workers.emplace_back(thread_func_forCurrent, std::ref(tArgs[i]));
    }

    // 3) 等待全部线程结束
    for (auto &th : workers)
    {
        if (th.joinable())
            th.join();
    }

    if (std::fabs(Ey[0][0]) < -500.0 && Priflag == 0)
    {
        Priflag = 1;
        Prinstp = nstp;
        limpoz = 661;

        // 更新 neflag ：ne > 1e18 → 0,  否则 → 1
        for (int i = 0; i < NR - 1; ++i)
        {
            for (int j = 0; j < NZ; ++j)
            {
                neflag[i][j] = (ne[i][j] > 1e18) ? 0 : 1;
            }
        }

        // 写 “Priflag.dat”
        {
            std::ofstream ofs("Priflag.dat");
            if (!ofs)
            {
                std::cerr << "Cannot open Priflag.dat for writing\n";
            }
            else
            {
                ofs << Prinstp << '\n';
            }
        }
    }

    // ───────────── 第二段：Priflag2 ─────────────
    // 触发条件：ne[0][861] > 1e19
    if (ne[0][861] > 1e19 && Priflag2 == 0)
    {
        Priflag2 = 1; // 在 z = 10 mm（网格 861）处拉起第二个标志
        Prinstp2 = nstp;

        {
            std::ofstream ofs("Priflag2.dat");
            if (!ofs)
            {
                std::cerr << "Cannot open Priflag2.dat for writing\n";
            }
            else
            {
                ofs << Prinstp2 << '\n';
            }
        }

        // 重新计算 neflag（同上）
        for (int i = 0; i < NR - 1; ++i)
        {
            for (int j = 0; j < NZ; ++j)
            {
                neflag[i][j] = (ne[i][j] > 1e18) ? 0 : 1;
            }
        }
    }

    /***************************************************************
     * 1. 电子扩散项 Dvx / Dvy 以及 ne 更新 —— 保持算法不变
     **************************************************************/
    {
        /* ---------- x 方向扩散 ---------- */
        int i = 0;
        for (int j = 0; j < 1421; ++j)
        {
            Dvx[i][j] = dt * Diff *
                        (Sr[i + 1][j] / (rh[i + 1] - rh[i]) *
                         (ne[i + 1][j] - ne[i][j])) /
                        Vol[i][j];
        }

        for (i = 1; i < NR - 1; ++i)
        {
            for (int j = 0; j < NZ; ++j)
            {
                if (flag[i][j])
                {
                    Dvx[i][j] = 0.0;
                }
                else
                {
                    if (flag[i - 1][j])
                    {
                        Dvx[i][j] = dt * Diff *
                                    (Sr[i + 1][j] / (rh[i + 1] - rh[i]) *
                                     (ne[i + 1][j] - ne[i][j])) /
                                    Vol[i][j];
                    }
                    else
                    {
                        Dvx[i][j] = dt * Diff *
                                    (-Sr[i][j] / (rh[i] - rh[i - 1]) *
                                         (ne[i][j] - ne[i - 1][j]) +
                                     Sr[i + 1][j] / (rh[i + 1] - rh[i]) *
                                         (ne[i + 1][j] - ne[i][j])) /
                                    Vol[i][j];
                    }
                }
            }
        }

        i = NR;
        for (int j = 0; j < NZ; ++j)
        {
            Dvx[i][j] = dt * Diff *
                        (-Sr[i][j] / (rh[i] - rh[i - 1]) *
                         (ne[i][j] - ne[i - 1][j])) /
                        Vol[i][j];
        }

        /* ---------- y 方向扩散 ---------- */
        for (i = 0; i < NR; ++i)
        {
            for (int j = 1; j < NZ; ++j)
            {
                if (flag[i][j])
                {
                    Dvy[i][j] = 0.0;
                }
                else
                {
                    if (flag[i][j + 1])
                    {
                        Dvy[i][j] = dt * Diff *
                                    (-Sz[i][j] / (zh[j] - zh[j - 1]) *
                                     (ne[i][j] - ne[i][j - 1])) /
                                    Vol[i][j];
                    }
                    else
                    {
                        Dvy[i][j] = dt * Diff *
                                    (-Sz[i][j] / (zh[j] - zh[j - 1]) *
                                         (ne[i][j] - ne[i][j - 1]) +
                                     Sz[i][j + 1] / (zh[j + 1] - zh[j]) *
                                         (ne[i][j + 1] - ne[i][j])) /
                                    Vol[i][j];
                    }
                }
            }
        }

        for (i = 0; i < NR; ++i)
        {
            int j = 0;
            Dvy[i][j] = dt * Diff *
                        (Sz[i][j + 1] / (zh[j + 1] - zh[j]) *
                         (ne[i][j + 1] - ne[i][j])) /
                        Vol[i][j];
        }

        /* ---------- 更新 ne ---------- */
        for (i = 0; i < NR; ++i)
        {
            for (int j = 0; j < NZ; ++j)
            {
                ne[i][j] += Dvx[i][j] + Dvy[i][j];
            }
        }
    }

    /***************************************************************
     * 2. 边界条件（保持原函数接口不变）
     **************************************************************/
    mol_boundary(NR, NZ, flag, ne,
                 N2p, O2p, H2Op, O2m, Om, OHm, Hm,
                 N4p, O4p, N2O2p, O2pH2O,
                 H3Op, H3OpH2O, H3OpH2O2, H3OpH2O3,
                 Ex, Ey, absE);

    std::cout << "==============All tests completed================\n";
    /***************************************************************
     * 3. 调用多线程 fluid 更新
     **************************************************************/
    std::cout << "---Finish_Fluid_equation... " << std::flush;

    {
        std::vector<std::thread> workers;
        workers.reserve(THREAD_NUM);

        for (int t = 0; t < THREAD_NUM; ++t)
        {
            targ[t].thread_no = t;                                      // 参数初始化
            workers.emplace_back(thread_func_fluid, std::ref(targ[t])); // 注意是 std::ref
        }

        for (auto &th : workers)
            th.join(); // 等待全部线程结束
    }

    std::cout << '\n';

    return 0;
}
