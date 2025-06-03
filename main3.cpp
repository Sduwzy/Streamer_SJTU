#include <iostream>
#include <cmath>
#include <cstdio>

// ===== 1) 引入头文件 =====
#include "memory.h"          // vec(), mat(), imat(), free_vec(), free_mat()
#include "spline.h"          // spline, splint, sp6
#include "calcE.h"           // calc_E(...)
#include "mesh_generator.h"  // mesh_generator(...)
#include "get_constant.h"    // get_constant()


// ===== 2) 定义与 get_constant(...) 相关宏及变量 =====
#define BOLTZMANNFILE  "inputdata/Bolsig_Data/O2_20p.dat"
#define FILE_MESH_Z "inputdata/mesh_z.dat"
#define FILE_MESH_R "inputdata/mesh_r.dat"

#define BN             1200    // ボルツマン方程式のデーター数
#define DATA           200     // reaction 数，需大于实际反应数
#define NUM_REACT      41      // (示例值)

// 在 get_constant() 中会用到的全局指针：
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


// ===== 3) store(...) / calcE(...) 所需的全局变量 =====

// 例如，这里给出一个简单网格大小:
int NR = 5;
int NZ = 5;

// 一些全局的 2D 数组指针 (store 和您原先的程序中会用到)：
double** absE  = nullptr;
double** Ex    = nullptr;
double** ne    = nullptr;
double** Cphi  = nullptr;
double** Lphi  = nullptr;
double** T     = nullptr;
double** p     = nullptr;
double** rou   = nullptr;

double** N2p   = nullptr;
double** N4p   = nullptr;
double** O2p   = nullptr;
double** O4p   = nullptr;
double** N2O2p = nullptr;
double** O2m   = nullptr;
double** Om    = nullptr;
double** H2Op  = nullptr;
double** OHm   = nullptr;
double** Hm    = nullptr;
double** O2pH2O= nullptr;
double** H3Op  = nullptr;
double** H3OpH2O  = nullptr;
double** H3OpH2O2 = nullptr;
double** H3OpH2O3 = nullptr;

// zh 数组供 store() 最后部分输出
double* zh     = nullptr;
int i, j;

// -----------------------------------------------------------------------
// 5) store(...) 函数的实现
//    (需要访问全局变量 absE, Ex, ne, N2p, Lphi, zh 等)
// -----------------------------------------------------------------------
void store(int Z, double time, double** Ey, char** particle, 
           double atmpa, double air_kg, double molV)
{
    FILE *fp;
    char filename[256];
    int skipnum = 1; // 简易演示

    if(Z % 1000 == 0)
    {
        // 1) absE
        std::sprintf(filename, "outputdata/2DE/2DE_%d.dat", Z);
        fp = std::fopen(filename, "w");
        if(!fp){
            std::cerr << "Cannot open " << filename << " for writing.\n";
            return;
        }
        for(i=0; i<NR; i+=skipnum){
            for(j=0; j<NZ; j+=skipnum){
                std::fprintf(fp, "%lf\n", absE[i][j]);
            }
        }
        std::fclose(fp);

        // 2) Ex
        std::sprintf(filename,"outputdata/2DEx/2DEx_%d.dat", Z);
        fp = std::fopen(filename,"w");
        if(!fp){
            std::cerr << "Cannot open " << filename << " for writing.\n";
            return;
        }
        for(i=0; i<NR; i+=skipnum){
            for(j=0; j<NZ; j+=skipnum){
                std::fprintf(fp, "%lf\n", Ex[i][j]);
            }
        }
        std::fclose(fp);

        // 3) Ey
        std::sprintf(filename,"outputdata/2DEy/2DEy_%d.dat", Z);
        fp = std::fopen(filename,"w");
        if(!fp){
            std::cerr << "Cannot open " << filename << " for writing.\n";
            return;
        }
        for(i=0; i<NR; i+=skipnum){
            for(j=0; j<NZ; j+=skipnum){
                std::fprintf(fp, "%lf\n", Ey[i][j]);
            }
        }
        std::fclose(fp);

        // 4) ne
        std::sprintf(filename,"outputdata/2Dne/ne_%d.dat", Z);
        fp = std::fopen(filename,"w");
        if(!fp){
            std::cerr << "Cannot open " << filename << " for writing.\n";
            return;
        }
        for(i=0; i<NR; i+=skipnum){
            for(j=0; j<NZ; j+=skipnum){
                std::fprintf(fp, "%e\n", ne[i][j]);
            }
        }
        std::fclose(fp);

        // 5) Cphi
        std::sprintf(filename,"outputdata/2DCphi/Cphi_%d.dat", Z);
        fp = std::fopen(filename,"w");
        if(!fp){
            std::cerr << "Cannot open " << filename << " for writing.\n";
            return;
        }
        for(i=0; i<NR; i+=skipnum){
            for(j=0; j<NZ; j+=skipnum){
                std::fprintf(fp, "%e\n", Cphi[i][j]);
            }
        }
        std::fclose(fp);

        // 6) Lphi
        std::sprintf(filename,"outputdata/2DLphi/Lphi_%d.dat", Z);
        fp = std::fopen(filename,"w");
        if(!fp){
            std::cerr << "Cannot open " << filename << " for writing.\n";
            return;
        }
        for(i=0; i<NR; i+=skipnum){
            for(j=0; j<NZ; j+=skipnum){
                std::fprintf(fp, "%e\n", Lphi[i][j]);
            }
        }
        std::fclose(fp);
    }

    // 下面的 if(0) 不会执行，只是保留原始代码以供参考
    if(0){
        // ... (输出 T, p, rou 等)
    }

    // 最后一个演示：ne_%d.dat
    // (注意这里硬编码了 j<1430，若 NZ < 1430 会越界)
    std::sprintf(filename,"outputdata/Ey_ne/ne_%d.dat", Z);
    fp = std::fopen(filename, "w");
    if(!fp){
        std::cerr << "Cannot open " << filename << " for writing.\n";
        return;
    }
    i=0;
    for(j=0; j<1430; j++){
        if(j >= NZ){
            // 若 j>=NZ，会越界，这里可以 break 或 continue
            break;
        }
        std::fprintf(fp, "%f\t", zh[j]*1000.0);
        std::fprintf(fp, "%f\t", Ey[i][j]);
        std::fprintf(fp, "%e\t", ne[i][j]);
        std::fprintf(fp, "%e\t", N2p[i][j]);
        std::fprintf(fp, "%e\t", N4p[i][j]);
        std::fprintf(fp, "%e\t", O2p[i][j]);
        std::fprintf(fp, "%e\t", O4p[i][j]);
        std::fprintf(fp, "%e\t", N2O2p[i][j]);
        std::fprintf(fp, "%e\t", O2m[i][j]);
        std::fprintf(fp, "%e\t", Om[i][j]);
        std::fprintf(fp, "%e\t", H2Op[i][j]);
        std::fprintf(fp, "%e\t", OHm[i][j]);
        std::fprintf(fp, "%e\t", Hm[i][j]);
        std::fprintf(fp, "%e\t", O2pH2O[i][j]);
        std::fprintf(fp, "%e\t", H3Op[i][j]);
        std::fprintf(fp, "%e\t", H3OpH2O[i][j]);
        std::fprintf(fp, "%e\t", H3OpH2O2[i][j]);
        std::fprintf(fp, "%e\n", H3OpH2O3[i][j]);
    }
    std::fclose(fp);

    std::cout << "store(...) done with Z=" << Z << std::endl;
}



// ====== 6) 主函数实现，用来测试上面的模块 ======
int main()
{
    // 1) 测试 memory.h
    std::cout << "===== 1) Test memory.h =====\n";
    {
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

    // 2) 测试 get_constant()
    std::cout << "\n===== 2) Test get_constant() =====\n";
    {
        std::cout << "Calling get_constant()...\n";
        get_constant(); // 读入碰撞截面或相关Bolsig文件
        std::cout << "Finished get_constant.\n";

        if(dc != nullptr){
            std::cout << "dc allocated with size DATA x BN = " 
                      << DATA << " x " << BN << "\n";
            std::cout << "dc[1][1] = " << dc[1][1] << "\n";
        }
    }

    // 3) 测试 spline()
    std::cout << "\n===== 3) Test spline() =====\n";
    {
        int n=5;
        double* x  = vec(n);
        double* y  = vec(n);
        double* y2 = vec(n);

        for(int i=0; i<n; i++){
            x[i] = (double)i;
            y[i] = std::sin(x[i]);
        }
        spline(x, y, n, y2);

        double Xval = 2.5, Yval=0.0;
        splint(x, y, y2, n, Xval, &Yval);
        std::cout << "splint at x=2.5 => " << Yval 
                  << ", actual sin(2.5)=" << std::sin(2.5) << "\n";

        free_vec(x, n);
        free_vec(y, n);
        free_vec(y2, n);
    }

    // 4) 测试 mesh_generator()
    std::cout << "\n===== 4) Test mesh_generator =====\n";
    {
        // 这里我们已经将 NR,NZ 设为全局 =5
        double* r   = vec(NR);
        double* z   = vec(NZ);
        double* rh2 = vec(NR);
        double* zh2 = vec(NZ);
        double* dr2 = vec(NR);
        double* dz2 = vec(NZ);
        double** Vol2= mat(NR,NZ);
        double** Sr2 = mat(NR,NZ);
        double** Sz2 = mat(NR,NZ);

        // 若您没有真实 mesh_r.dat / mesh_z.dat，可自己生成简单测试文件
        mesh_generator(FILE_MESH_R, FILE_MESH_Z,
                       NR,NZ, r,z, rh2,zh2, dr2,dz2,
                       Vol2,Sr2,Sz2);

        std::cout << "r[] from mesh:\n";
        for(int i=0; i<NR; i++){
            std::cout << i << ": " << r[i] << "\n";
        }
        std::cout << "z[] from mesh:\n";
        for(int j=0; j<NZ; j++){
            std::cout << j << ": " << z[j] << "\n";
        }

        // 释放
        free_vec(r,NR);
        free_vec(z,NZ);
        free_vec(rh2,NR);
        free_vec(zh2,NZ);
        free_vec(dr2,NR);
        free_vec(dz2,NZ);

        free_mat(Vol2,NR,NZ);
        free_mat(Sr2, NR,NZ);
        free_mat(Sz2, NR,NZ);
    }

    // 5) 分配全局数组, 测试calcE
    std::cout << "\n===== 5) Alloc & Test calcE + store =====\n";
    {
        // 分配全局数组: absE, Ex, ne, Cphi, Lphi, T, p, rou, ...
        absE  = mat(NR,NZ);
        Ex    = mat(NR,NZ);
        ne    = mat(NR,NZ);
        Cphi  = mat(NR,NZ);
        Lphi  = mat(NR,NZ);
        T     = mat(NR,NZ);
        p     = mat(NR,NZ);
        rou   = mat(NR,NZ);

        // 分配其他粒子矩阵
        N2p    = mat(NR,NZ);
        N4p    = mat(NR,NZ);
        O2p    = mat(NR,NZ);
        O4p    = mat(NR,NZ);
        N2O2p  = mat(NR,NZ);
        O2m    = mat(NR,NZ);
        Om     = mat(NR,NZ);
        H2Op   = mat(NR,NZ);
        OHm    = mat(NR,NZ);
        Hm     = mat(NR,NZ);
        O2pH2O = mat(NR,NZ);
        H3Op   = mat(NR,NZ);
        H3OpH2O  = mat(NR,NZ);
        H3OpH2O2 = mat(NR,NZ);
        H3OpH2O3 = mat(NR,NZ);

        zh = vec(NZ);

        // 初始化一些数据，便于store查看
        for(int i=0; i<NR; i++){
            for(int j=0; j<NZ; j++){
                absE[i][j] = 0.1*(i+j);
                Ex[i][j]   = 0.2*(i+j);
                ne[i][j]   = 1e12 + 1e6*(i+j);
                Cphi[i][j] = 100.0 + 0.5*(i+j);
                Lphi[i][j] = 200.0 + 0.3*(i+j);
                T[i][j]    = 300.0;
                p[i][j]    = 1.0e5;
                rou[i][j]  = 1.225;

                // 其他离子浓度
                N2p[i][j]    = 1.0*(i+j);
                N4p[i][j]    = 2.0*(i+j);
                O2p[i][j]    = 3.0*(i+j);
                O4p[i][j]    = 4.0*(i+j);
                N2O2p[i][j]  = 5.0*(i+j);
                O2m[i][j]    = 6.0*(i+j);
                Om[i][j]     = 7.0*(i+j);
                H2Op[i][j]   = 8.0*(i+j);
                OHm[i][j]    = 9.0*(i+j);
                Hm[i][j]     = 10.0*(i+j);
                O2pH2O[i][j] = 11.0*(i+j);
                H3Op[i][j]   = 12.0*(i+j);
                H3OpH2O[i][j]= 13.0*(i+j);
                H3OpH2O2[i][j]=14.0*(i+j);
                H3OpH2O3[i][j]=15.0*(i+j);
            }
        }
        for(int j=0; j<NZ; j++){
            zh[j] = 0.001*j; // 以 m 计
        }

        // 例如调用calcE(...) (伪代码,只演示)
        // int** totuflag = imat(NR,NZ); // ...
        // calc_E(NR,NZ, Cphi, absE, Ey, Ex,
        //        28.966/1000.0, totuflag, ...);

        // 最后: 调用 store(...) 函数输出到文件
        std::cout << "Call store(...) now...\n";

        // 先做个Ey分配/初始化
        double** Ey = mat(NR,NZ);
        for(int i=0; i<NR; i++){
            for(int j=0; j<NZ; j++){
                Ey[i][j] = 99.9*(i+j);
            }
        }

        // 假设没有真正的particle表，随便给个空指针
        char** dummyParticle = nullptr;

        // 调用 store
        store(1000, 0.123, Ey, dummyParticle, /*atmpa*/1.013e5, /*air_kg*/0.028966, /*molV*/0.0224);
        // 以上数字纯演示

        // 释放 Ey
        free_mat(Ey,NR,NZ);

        // 释放 dummyParticle 如有的话 (这里是nullptr, 无需释放)
        
        // ...结束
    }

    std::cout << "\nAll tests done successfully.\n";
    return 0;
}
