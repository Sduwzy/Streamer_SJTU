#include "get_constant.h"
#include <cstdio>   // fopen, fscanf, fclose
#include <cstdlib>  // exit
#include <cmath>    // 如果需要
#include <iostream> // 如果需要用 std::cout
// 如果您用到 "spline.h" 之类，也要包含
#include "spline.h"


#define BOLTZMANNFILE	"inputdata/Bolsig_Data/O2_20p.dat"
#define BN	1200	//ボルツマン方程式のデーター数
#define DATA 200 //テキトーに入れておく。反応数(C53とか)より多ければok
#define NUM_REACT 41 //C53+1

// 关键：需要对 main.cu 里定义的全局变量进行 extern 声明：
//extern int  AN, BN, DATA, NUM_REACT;
// 同理，如果 AN, etc. 在别处定义，也要 extern

// 如果 dc, ddc, TdE, deV, dv, 这些都是全局变量，也要 extern
extern double **dc, **ddc;
extern double *TdE, *deV, *dv;
extern double *ddeV, *ddv, *alpTd, *dalp, *ddalp;
extern double *powerTd, *dpower, *ddpower;



int AN=1200;




// 同理, spline(...) 函数若在 "spline.cpp" 中，需要能够链接到

void get_constant()
{
    int i, j;
    FILE *fp;

    // 分配 
    dc  = mat(DATA, BN);
    ddc = mat(DATA, BN);
    TdE = vec(BN);
    deV = vec(BN);
    dv  = vec(BN);
    ddeV= vec(BN);
    ddv = vec(BN);
    alpTd=vec(AN);
    dalp =vec(AN);
    ddalp=vec(AN);

    powerTd = vec(BN);
    dpower  = vec(BN);
    ddpower = vec(BN);

    // 1) 读 BOLTZMANNFILE
    fp = std::fopen(BOLTZMANNFILE, "r");
    if(!fp){
        std::perror("Cannot open BOLTZMANNFILE");
        std::exit(EXIT_FAILURE);
    }

    for(i=0; i<BN; i++){
        if(std::feof(fp)){
            std::fprintf(stderr, "Wrong Num of BN!!\n");
            std::exit(EXIT_FAILURE);
        }
        std::fscanf(fp,"%lf",&TdE[i]);
        std::fscanf(fp,"%lf",&deV[i]);
        std::fscanf(fp,"%lf",&dv[i]);
        for(j=1; j<NUM_REACT; j++){
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
    spline(TdE, dv,   BN, ddv);
    spline(TdE, deV,  BN, ddeV);
    for(i=1; i<NUM_REACT; i++){
        spline(TdE, dc[i], BN, ddc[i]);
    }
    // spline(alpTd, dalp, AN, ddalp); // 如果要

    // 3) 读 inputdata/PowerEdep.dat
    fp = std::fopen("inputdata/PowerEdep.dat", "r");
    if(!fp){
        std::perror("Cannot open inputdata/PowerEdep.dat");
        std::exit(EXIT_FAILURE);
    }
    for(i=0; i<BN; i++){
        std::fscanf(fp, "%lf", &powerTd[i]);
        std::fscanf(fp, "%lf", &dpower[i]);
    }
    std::fclose(fp);

    // 4) spline
    spline(powerTd, dpower, BN, ddpower);

    // 到此和原 get_constant() 相同
    std::cout << "get_constant() finished reading constants.\n";
}
