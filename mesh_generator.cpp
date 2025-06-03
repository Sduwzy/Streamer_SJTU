#include "mesh_generator.h"
#include <cstdio>   // fscanf, fopen, fclose
#include <cstdlib>  // NULL
#include <cmath>

void mesh_generator(const char *FILE_MESH_R, const char *FILE_MESH_Z,
                    int NR, int NZ,
                    double *r, double *z,
                    double *rh, double *zh,
                    double *dr, double *dz,
                    double **Vol, double **Sr, double **Sz)
{
    int i, j, dum;
    FILE *fp = nullptr;

    // 1) 读取 z 网格
    fp = std::fopen(FILE_MESH_Z, "r");
    if(!fp){
        std::perror("Failed to open FILE_MESH_Z");
        return; // 或者 throw std::runtime_error("Cannot open FILE_MESH_Z");
    }
    for(j=0; j<NZ; j++){
        // 原代码 fscanf(fp, "%d\t%le\n", &dum, &z[j]);
        // 这里用 %d %lf: le 对 double 是可行的
        std::fscanf(fp, "%d %lf", &dum, &z[j]);
    }
    std::fclose(fp);

    // 2) 读取 r 网格
    fp = std::fopen(FILE_MESH_R, "r");
    if(!fp){
        std::perror("Failed to open FILE_MESH_R");
        return;
    }
    for(i=0; i<NR; i++){
        std::fscanf(fp, "%d %lf", &dum, &r[i]);
    }
    std::fclose(fp);

    // 3) 计算中心点 rh[i] = (r[i] + r[i+1]) / 2, zh[j] = ...
    //   但原代码中写法对 i=NR-1,j=NZ-1 做了特殊处理
    //   这里保持逻辑一致

    // 先 for(i=0..NR-2), for(j=0..NZ-2)
    for(i=0; i<NR-1; i++){
        for(j=0; j<NZ-1; j++){
            rh[i] = (r[i] + r[i+1]) * 0.5; 
            zh[j] = (z[j] + z[j+1]) * 0.5;
        }
    }

    // j=NZ-1
    {
        j = NZ-1;
        for(i=0; i<NR-1; i++){
            rh[i] = (r[i] + r[i+1]) * 0.5;
            // zh[j] = (z[j] + z[j]) * 0.5 => 就是 z[j]
            // 可能 z[j] = z[j]? => 不变
            zh[j] = (z[j] + z[j]) * 0.5;
        }
    }

    // i=NR-1
    {
        i = NR-1;
        for(j=0; j<NZ-1; j++){
            rh[i] = (r[i] + r[i]) * 0.5; 
            zh[j] = (z[j] + z[j+1]) * 0.5;
        }
    }

    // i=NR-1, j=NZ-1
    {
        i=NR-1; j=NZ-1;
        rh[i] = (r[i]+r[i]) * 0.5;
        zh[j] = (z[j]+z[j]) * 0.5;
    }

    // 4) 计算 dr[i], dz[j]
    //   dr[i] = r[i+1] - r[i], except for the last i => repeat
    for(i=0; i<NR-1; i++){
        dr[i] = r[i+1] - r[i];
    }
    dr[NR-1] = dr[NR-2];

    for(j=0; j<NZ-1; j++){
        dz[j] = z[j+1] - z[j];
    }
    dz[NZ-1] = dz[NZ-2];

    // 5) 计算 Vol[i][j], Sr[i][j], Sz[i][j]
    //   原逻辑: 
    //   Vol[i][j] = dz[j]*(r[i]+r[i+1])*dr[i]/2.0
    //   Sr[i][j]  = r[i]*dz[j]
    //   Sz[i][j]  = (r[i]+r[i+1])*dr[i]/2.0
    //   并且 i=NR-1 要特殊 if
    for(i=0; i<NR; i++){
        for(j=0; j<NZ; j++){
            if(i == NR-1){
                // Vol: dz[j]*(r[i]+r[i])*dr[i]/2.0
                Vol[i][j] = dz[j]*(r[i]+r[i])*dr[i] * 0.5;
                Sr[i][j]  = r[i]*dz[j]; 
                Sz[i][j]  = (r[i]+r[i])*dr[i]*0.5;
            } else {
                Vol[i][j] = dz[j]*(r[i]+r[i+1])*dr[i]*0.5;
                Sr[i][j]  = r[i]*dz[j];
                Sz[i][j]  = (r[i]+r[i+1])*dr[i]*0.5;
            }
        }
    }
}
