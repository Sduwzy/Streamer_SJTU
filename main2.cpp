#include <iostream>
#include <cmath>
#include <cstdio>

// 包含您各模块的头文件
#include "memory.h"          // vec(), mat(), imat(), free_vec(), free_mat()
#include "mesh_generator.h"  // mesh_generator(...)
#include "spline.h"          // spline, splint, sp6
#include "calcE.h"           // calc_E(...)
#include "get_constant.h"    // get_constant()

// 如果您需要别的整型也可加在这

#define BOLTZMANNFILE	"inputdata/Bolsig_Data/O2_20p.dat"
#define BN	1200	//ボルツマン方程式のデーター数
#define DATA 200 //テキトーに入れておく。反応数(C53とか)より多ければok
#define NUM_REACT 41 //C53+1
#define FILE_MESH_Z "inputdata/mesh_z.dat"
#define FILE_MESH_R "inputdata/mesh_r.dat"

// 全局指针 (extern double**) 在此真正定义
double **dc   = nullptr;
double **ddc  = nullptr;
double *TdE   = nullptr;
double *deV   = nullptr;
double *dv    = nullptr;
double *ddeV  = nullptr;
double *ddv   = nullptr;
double *alpTd = nullptr;
double *dalp  = nullptr;
double *ddalp = nullptr;
double *powerTd = nullptr;
double *dpower  = nullptr;
double *ddpower = nullptr;
int NR=4, NZ=5;
// ====== 头文件引用 ======
#include "memory.h"          // vec(), mat(), free_vec(), free_mat(), ...
#include "get_constant.h"    // get_constant()
#include "spline.h"          // spline(), splint(), sp6
#include "mesh_generator.h"  // mesh_generator(...)
#include "calcE.h"           // calc_E(...)

//============================================================
// main 函数：依次测试您已有的模块
//============================================================
int main()
{
    std::cout << "===== 1) Test memory.h =====\n";
    {
        // 测试 vec(...) / free_vec(...)
        double* testVec = vec(5);
        for(int i=0; i<5; i++){
            testVec[i] = i*0.1;
        }
        std::cout << "testVec: ";
        for(int i=0; i<5; i++){
            std::cout << testVec[i] << " ";
        }
        std::cout << std::endl;
        free_vec(testVec, 5);

        // 测试 mat(...) / free_mat(...)
        double** testMat = mat(3, 3);
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                testMat[i][j] = i*10 + j;
            }
        }
        std::cout << "testMat:\n";
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                std::cout << testMat[i][j] << " ";
            }
            std::cout << "\n";
        }
        free_mat(testMat, 3, 3);
    }

    std::cout << "\n===== 2) Test get_constant() =====\n";
    {
        // 调用 get_constant(), 里面会用到 BN,DATA,NUM_REACT,BOLTZMANNFILE 等
        // 并给 dc, ddv, ... 分配内存并读文件(若文件存在)
        std::cout << "Calling get_constant()...\n";
        get_constant();
        // 这里若文件 data 不存在则会报错，请根据实际路径
        std::cout << "Finished get_constant.\n";

        // (可选) 展示一下 dc, dv, etc. 里是否有数据
        //   例如：
        if(dc != nullptr){
            std::cout << "dc allocated with size DATA x BN = " 
                      << DATA << " x " << BN << "\n";
            // 打印 dc[1][1] 之类
            std::cout << "dc[1][1] = " << dc[1][1] << "\n";
        }
    }

    std::cout << "\n===== 3) Test spline (spline.h) =====\n";
    {
        int n=5;
        double* x  = vec(n);
        double* y  = vec(n);
        double* y2 = vec(n);

        // 简单测试：x=0,1,2,3,4; y= sin(x)
        for(int i=0; i<n; i++){
            x[i] = (double)i;
            y[i] = std::sin(x[i]);
        }

        spline(x, y, n, y2);
        double Xval = 2.5, Yval;
        splint(x, y, y2, n, Xval, &Yval);
        std::cout << "spline at x=2.5 => " << Yval 
                  << ", actual sin(2.5)=" << sin(2.5) << "\n";

        free_vec(x, n);
        free_vec(y, n);
        free_vec(y2, n);
    }

    std::cout << "\n===== 4) Test mesh_generator =====\n";
    {

        double* r   = vec(NR);
        double* z   = vec(NZ);
        double* rh  = vec(NR);
        double* zh  = vec(NZ);
        double* dr  = vec(NR);
        double* dz  = vec(NZ);
        double** Vol= mat(NR,NZ);
        double** Sr = mat(NR,NZ);
        double** Sz = mat(NR,NZ);

        // 需要 mesh_r.dat / mesh_z.dat
        // 此文件必须和 mesh_generator 中读取方式对应
        // 若没有, 可自己创建简单数据
        mesh_generator(FILE_MESH_R,FILE_MESH_Z,
                       NR,NZ,
                       r,z,
                       rh,zh,
                       dr,dz,
                       Vol,Sr,Sz);

        // 打印 r,z (演示)
        std::cout << "r[]:\n";
        for(int i=0; i<NR; i++){
            std::cout << i << ": " << r[i] << "\n";
        }
        std::cout << "z[]:\n";
        for(int j=0; j<NZ; j++){
            std::cout << j << ": " << z[j] << "\n";
        }

        free_vec(r,NR);
        free_vec(z,NZ);
        free_vec(rh,NR);
        free_vec(zh,NZ);
        free_vec(dr,NR);
        free_vec(dz,NZ);

        free_mat(Vol,NR,NZ);
        free_mat(Sr, NR,NZ);
        free_mat(Sz, NR,NZ);
    }

    std::cout << "\n===== 5) Test calcE =====\n";
    {
        // 仅演示：实际使用calcE时, 需要更多参数
        // 假设 NR=3, NZ=3

        int NR=3, NZ=3;
        double** phi   = mat(NR,NZ);
        double** absE  = mat(NR,NZ);
        double** Ex    = mat(NR,NZ);
        double** Ey    = mat(NR,NZ);
        double air_kg  = 28.966/1000.0;

        int** totuflag = imat(NR,NZ);
        int** otuflag  = imat(NR,NZ);
        int** iflag    = imat(NR,NZ);
        int** jflag    = imat(NR,NZ);
        int** flag     = imat(NR,NZ);

        // 简单初始化
        for(int i=0; i<NR; i++){
            for(int j=0; j<NZ; j++){
                phi[i][j] = (i==1 && j==1)? 100.0 : 0.0;
                absE[i][j]=0.0;
                Ex[i][j]  =0.0;
                Ey[i][j]  =0.0;
                totuflag[i][j]=0;
                otuflag[i][j]=0;
                iflag[i][j]=0;
                jflag[i][j]=0;
                flag[i][j]=0;
            }
        }

        // 如果calc_E需要 a,b,rhalf,zhalf,Mol等，请自行设置
        // 这里仅演示调用
        std::cout << "call calc_E(...) for demonstration.\n";
        // calc_E(NR,NZ,phi,absE,Ey,Ex, air_kg, totuflag, otuflag, iflag, jflag, flag,
        //         2.0,3.0, rhalf, zhalf, ???, ???, ???);

        // 打印Ex[1][1]
        std::cout << "Ex[1][1]=" << Ex[1][1] << "\n";

        free_mat(phi,  NR,NZ);
        free_mat(absE, NR,NZ);
        free_mat(Ex,   NR,NZ);
        free_mat(Ey,   NR,NZ);

        free_mat((double**)totuflag, NR,NZ);
        free_mat((double**)otuflag,  NR,NZ);
        free_mat((double**)iflag,    NR,NZ);
        free_mat((double**)jflag,    NR,NZ);
        free_mat((double**)flag,     NR,NZ);
    }

    std::cout << "\nAll tests done successfully.\n";
    return 0;
}