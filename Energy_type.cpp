#include "Energy_type.h"

// convert_M: 将 O2, N2, H2O 各9项 (0~8) 合并至 M
void convert_M(double* O2, double* N2, double* H2O, double* M)
{
    int Q = 8;
    for(int i = 0; i <= Q; i++){
        M[i]                = O2[i];
        M[i + (Q+1)]        = N2[i];
        M[i + 2*(Q+1)]      = H2O[i];
    }
}

// convert_ONH2O: 将大数组M(长度3*(Q+1))拆分为 O2[], N2[], H2O[]
void convert_ONH2O(double* O2, double* N2, double* H2O, double* M)
{
    int Q = 8;
    for(int i = 0; i <= Q; i++){
        O2[i]  = M[i];
        N2[i]  = M[i + (Q+1)];
        H2O[i] = M[i + 2*(Q+1)];
    }
}

// vib_relaxation: 计算网格(i,j)处的振动弛豫过程
void vib_relaxation(
    int i, 
    int j,
    double* M,     
    double* dM,     
    double Oz,      
    double** kvt,
    double** kvv0, double** kvv1, double** kvv2, double** kvv3, double** kvv4,
    double** re_kvt,
    double** re_kvv0, double** re_kvv1, double** re_kvv2, double** re_kvv3, double** re_kvv4,
    double T,
    double** E
)
{
    // ---------------------------------------------------------------------
    // 这里直接贴上您之前的核心逻辑(并做适当 C++ 风格修饰)
    // 注意：若 i,j 并未用到，这里只是保留其形参签名
    // ---------------------------------------------------------------------

    int Q = 8;

    // 可能的常数
    double koz   = 3.2e-12;
    double kN2oz = (2.3e-13*std::exp(-1280.0 / T) + 2.7e-11*std::exp(-10840.0 / T)) * 1e-6;

    // 分配临时数组
    double* O2  = vec(20);
    double* N2  = vec(20);
    double* H2O = vec(20);

    double* O2VT     = vec(Q+5); 
    double* N2VT     = vec(Q+5); 
    double* H2OVT    = vec(Q+5); 
    double* O2OVT    = vec(Q+5);
    double* reO2VT   = vec(Q+5); 
    double* reN2VT   = vec(Q+5); 
    double* reH2OVT  = vec(Q+5); 
    double* reO2OVT  = vec(Q+5);

    double* O2N2VT   = vec(Q+5); 
    double* N2O2VT   = vec(Q+5); 
    double* reO2N2VT = vec(Q+5);
    double* reN2O2VT = vec(Q+5);

    double* N2OVT    = vec(Q+5); 
    double* reN2OVT  = vec(Q+5);

    double* O2VV      = vec(Q+5); 
    double* N2VV      = vec(Q+5); 
    double* O2N2VV    = vec(Q+5); 
    double* N2O2VV    = vec(Q+5);
    double* O2H2OVV   = vec(Q+5); 
    double* H2OO2VV   = vec(Q+5); 
    double* N2H2OVV   = vec(Q+5); 
    double* H2ON2VV   = vec(Q+5);

    double* reO2VV    = vec(Q+5); 
    double* reN2VV    = vec(Q+5); 
    double* reO2N2VV  = vec(Q+5);
    double* reN2O2VV  = vec(Q+5);
    double* reO2H2OVV = vec(Q+5); 
    double* reH2OO2VV = vec(Q+5);
    double* reN2H2OVV = vec(Q+5); 
    double* reH2ON2VV = vec(Q+5);

    // 全部清0
    for(int w=0; w<=Q; w++){
        O2VT[w]    = N2VT[w]    = H2OVT[w]   = O2OVT[w]   =
        reO2VT[w]  = reN2VT[w]  = reH2OVT[w] = reO2OVT[w] =
        O2N2VT[w]  = N2O2VT[w]  = reO2N2VT[w]= reN2O2VT[w]=
        N2OVT[w]   = reN2OVT[w] = 0.0;

        O2VV[w]      = N2VV[w]      = O2N2VV[w]    = N2O2VV[w]    =
        O2H2OVV[w]   = H2OO2VV[w]   = N2H2OVV[w]   = H2ON2VV[w]   =
        reO2VV[w]    = reN2VV[w]    = reO2N2VV[w]  = reN2O2VV[w]  =
        reO2H2OVV[w] = reH2OO2VV[w] = reN2H2OVV[w] = reH2ON2VV[w] =
        0.0;
    }
    for(int v=0; v<10; v++){
        O2[v] = N2[v] = H2O[v] = 0.0;
    }

    // 拆分 M => O2[],N2[],H2O[]
    convert_ONH2O(O2, N2, H2O, M);

    // ---------------------------
    // （下面照您之前的逻辑写）
    // ---------------------------

    // 1) V-T
    for(int n=0; n<=Q; n++){
        O2VT[n]  += kvt[0][n+1]*O2[n+1]*O2[0]  - kvt[0][n]*O2[n]*O2[0];
        N2VT[n]  += kvt[1][n+1]*N2[n+1]*N2[0]  - kvt[1][n]*N2[n]*N2[0];
        H2OVT[n] += kvt[3][n+1]*H2O[n+1]*H2O[0]- kvt[3][n]*H2O[n]*H2O[0];
        O2OVT[n] += kvt[5][n+1]*O2[n+1]*Oz     - kvt[5][n]*O2[n]*Oz;
    }
    // 2) O2-N2
    for(int n=0; n<=Q; n++){
        O2N2VT[n] += kvt[2][n+1]*O2[n+1]*N2[0] - kvt[2][n]*O2[n]*N2[0];
        N2O2VT[n] += kvt[4][n+1]*N2[n+1]*O2[0] - kvt[4][n]*N2[n]*O2[0];
    }
    // 3) N2 - O(=Oz)
    N2OVT[0] +=    kN2oz * N2[1]*Oz;
    N2OVT[1] +=  - kN2oz * N2[1]*Oz;

    // 4) V-V 过程 (O2-O2, N2-N2, O2-N2, ... H2O ... )
    for(int w=0; w<Q; w++){
        for(int v=0; v<Q; v++){
            // 依次对 O2VV, N2VV, O2N2VV, ...
            // 并加上 reO2VV, reN2VV, ...
            // (省略注释，只保留代码)
            O2VV[v+1] -= kvv0[w][v+1]*O2[v+1]*O2[w];
            O2VV[w]   -= kvv0[w][v+1]*O2[v+1]*O2[w];
            O2VV[v]   += kvv0[w][v+1]*O2[v+1]*O2[w];
            O2VV[w+1] += kvv0[w][v+1]*O2[v+1]*O2[w];

            reO2VV[v]   -= re_kvv0[w][v+1]*O2[v]*O2[w+1];
            reO2VV[w+1] -= re_kvv0[w][v+1]*O2[v]*O2[w+1];
            reO2VV[v+1] += re_kvv0[w][v+1]*O2[v]*O2[w+1];
            reO2VV[w]   += re_kvv0[w][v+1]*O2[v]*O2[w+1];

            N2VV[v+1] -= kvv1[w][v+1]*N2[v+1]*N2[w];
            N2VV[w]   -= kvv1[w][v+1]*N2[v+1]*N2[w];
            N2VV[v]   += kvv1[w][v+1]*N2[v+1]*N2[w];
            N2VV[w+1] += kvv1[w][v+1]*N2[v+1]*N2[w];

            reN2VV[v]   -= re_kvv1[w][v+1]*N2[v]*N2[w+1];
            reN2VV[w+1] -= re_kvv1[w][v+1]*N2[v]*N2[w+1];
            reN2VV[v+1] += re_kvv1[w][v+1]*N2[v]*N2[w+1];
            reN2VV[w]   += re_kvv1[w][v+1]*N2[v]*N2[w+1];

            N2O2VV[v+1] -= kvv2[w][v+1]*N2[v+1]*O2[w];
            O2N2VV[w]   -= kvv2[w][v+1]*N2[v+1]*O2[w];
            N2O2VV[v]   += kvv2[w][v+1]*N2[v+1]*O2[w];
            O2N2VV[w+1] += kvv2[w][v+1]*N2[v+1]*O2[w];

            reN2O2VV[v]   -= re_kvv2[w][v+1]*N2[v]*O2[w+1];
            reO2N2VV[w+1] -= re_kvv2[w][v+1]*N2[v]*O2[w+1];
            reN2O2VV[v+1] += re_kvv2[w][v+1]*N2[v]*O2[w+1];
            reO2N2VV[w]   += re_kvv2[w][v+1]*N2[v]*O2[w+1];

            O2H2OVV[v+1] -= kvv3[w][v+1]*O2[v+1]*H2O[w];
            H2OO2VV[w]   -= kvv3[w][v+1]*O2[v+1]*H2O[w];
            O2H2OVV[v]   += kvv3[w][v+1]*O2[v+1]*H2O[w];
            H2OO2VV[w+1] += kvv3[w][v+1]*O2[v+1]*H2O[w];

            reO2H2OVV[v]   -= re_kvv3[w][v+1]*O2[v]*H2O[w+1];
            reH2OO2VV[w+1] -= re_kvv3[w][v+1]*O2[v]*H2O[w+1];
            reO2H2OVV[v+1] += re_kvv3[w][v+1]*O2[v]*H2O[w+1];
            reH2OO2VV[w]   += re_kvv3[w][v+1]*O2[v]*H2O[w+1];

            N2H2OVV[v+1] -= kvv4[w][v+1]*N2[v+1]*H2O[w];
            H2ON2VV[w]   -= kvv4[w][v+1]*N2[v+1]*H2O[w];
            N2H2OVV[v]   += kvv4[w][v+1]*N2[v+1]*H2O[w];
            H2ON2VV[w+1] += kvv4[w][v+1]*N2[v+1]*H2O[w];

            reN2H2OVV[v]   -= re_kvv4[w][v+1]*N2[v]*H2O[w+1];
            reH2ON2VV[w+1] -= re_kvv4[w][v+1]*N2[v]*H2O[w+1];
            reN2H2OVV[v+1] += re_kvv4[w][v+1]*N2[v]*H2O[w+1];
            reH2ON2VV[w]   += re_kvv4[w][v+1]*N2[v]*H2O[w+1];
        }
    }

    // 5) reO2VT, reN2VT ...
    for(int n=0; n<=Q; n++){
        if(n == 0){
            reO2VT[n]  += -re_kvt[0][1]*O2[0]*O2[0];
            reN2VT[n]  += -re_kvt[1][1]*N2[0]*N2[0];
            reO2N2VT[n]+= -re_kvt[2][1]*O2[0]*N2[0];
            reN2O2VT[n]+= -re_kvt[4][1]*N2[0]*O2[0];
            reH2OVT[n] += -re_kvt[3][1]*H2O[0]*H2O[0];
            reO2OVT[n] += -re_kvt[5][1]*O2[0]*Oz;
        }
        else {
            reO2VT[n]  +=  re_kvt[0][n]*O2[n-1]*O2[0]   - re_kvt[0][n+1]*O2[n]*O2[0];
            reN2VT[n]  +=  re_kvt[1][n]*N2[n-1]*N2[0]   - re_kvt[1][n+1]*N2[n]*N2[0];
            reO2N2VT[n]+=  re_kvt[2][n]*O2[n-1]*N2[0]   - re_kvt[2][n+1]*O2[n]*N2[0];
            reN2O2VT[n]+=  re_kvt[4][n]*N2[n-1]*O2[0]   - re_kvt[1][n+1]*N2[n]*O2[0];
            reH2OVT[n] +=  re_kvt[3][n]*H2O[n-1]*H2O[0] - re_kvt[3][n+1]*H2O[n]*H2O[0];
            reO2OVT[n] +=  re_kvt[5][n]*O2[n-1]*Oz      - re_kvt[5][n+1]*O2[n]*Oz;
        }
    }
    reN2OVT[0] += -kN2oz*std::exp(-(E[1][1]-E[1][0]) / T)*N2[0]*Oz;
    reN2OVT[1] +=  kN2oz*std::exp(-(E[1][1]-E[1][0]) / T)*N2[0]*Oz;

    // 6) H2O(v1)/H2O(v3)
    double H2Ov1VT = 2.2e-17 * M[20] * H2O[0]; // H2O(v1)+H2O(0)=>2H2O(0)
    double H2Ov3VT = 2.2e-17 * M[21] * H2O[0]; // H2O(v3)+H2O(0)=>2H2O(0)

    // 7) 写入 dM
    for(int n=0; n<=Q; n++){
        // O2
        dM[n] = (O2VT[n]  + reO2VT[n])
              + (O2OVT[n] + reO2OVT[n])
              + (O2N2VT[n]+ reO2N2VT[n])
              + (O2VV[n]   + reO2VV[n])
              + (O2N2VV[n] + reO2N2VV[n])
              + (O2H2OVV[n]+ reO2H2OVV[n]);

        // N2
        dM[n + (Q+1)] = (N2VT[n]   + reN2VT[n])
                      + (N2OVT[n]  + reN2OVT[n])
                      + (N2O2VT[n] + reN2O2VT[n])
                      + (N2VV[n]   + reN2VV[n])
                      + (N2O2VV[n] + reN2O2VV[n])
                      + (N2H2OVV[n]+ reN2H2OVV[n]);

        // H2O
        dM[n + 2*(Q+1)] = (H2OVT[n] + reH2OVT[n])
                        + (H2OO2VV[n]+ reH2OO2VV[n])
                        + (H2ON2VV[n]+ reH2ON2VV[n]);
    }
    dM[18] += (H2Ov1VT + H2Ov3VT);
    dM[20]  = -H2Ov1VT;
    dM[21]  = -H2Ov3VT;

    // 释放内存
    free_vec(O2, 20);   free_vec(N2, 20);   free_vec(H2O, 20);
    free_vec(O2VT, Q+5);    free_vec(N2VT, Q+5);     free_vec(H2OVT, Q+5); 
    free_vec(O2OVT, Q+5);   free_vec(reO2VT, Q+5);   free_vec(reN2VT, Q+5);
    free_vec(reH2OVT, Q+5); free_vec(reO2OVT, Q+5);  free_vec(O2N2VT, Q+5);
    free_vec(N2O2VT, Q+5);  free_vec(reO2N2VT, Q+5); free_vec(reN2O2VT, Q+5);
    free_vec(N2OVT, Q+5);   free_vec(reN2OVT, Q+5);
    free_vec(O2VV, Q+5);    free_vec(N2VV, Q+5);    free_vec(O2N2VV, Q+5); 
    free_vec(N2O2VV, Q+5);  free_vec(O2H2OVV, Q+5); free_vec(H2OO2VV, Q+5);
    free_vec(N2H2OVV, Q+5); free_vec(H2ON2VV, Q+5);
    free_vec(reO2VV, Q+5);  free_vec(reN2VV, Q+5);  free_vec(reO2N2VV, Q+5);
    free_vec(reN2O2VV, Q+5); free_vec(reO2H2OVV, Q+5);
    free_vec(reH2OO2VV, Q+5); free_vec(reN2H2OVV, Q+5); free_vec(reH2ON2VV, Q+5);
}

// Read_constant: 读取相关文件，以填充 dE[][], E[][], c1[], c2[], c3[]
void Read_constant(double** dE, double** E, double* c1, double* c2, double* c3)
{
    // 与您之前的逻辑一致
    int i, j;
    char filename[256];
    std::FILE* fp = nullptr;

    // 先清 0
    for(i=0; i<5; i++){
        for(j=0; j<30; j++){
            dE[i][j] = 0.0;
        }
    }

    // 1) dEv_O2.dat
    fp = std::fopen("inputdata/dEv_O2.dat", "r");
    if(!fp){
        std::printf("Error: cannot open inputdata/dEv_O2.dat\n");
        std::exit(1);
    }
    for(i=1; i<=10; i++){
        std::fscanf(fp, "%lf\n", &dE[0][i]);
    }
    std::fclose(fp);

    // 2) dEv_N2.dat
    fp = std::fopen("inputdata/dEv_N2.dat", "r");
    if(!fp){
        std::printf("Error: cannot open inputdata/dEv_N2.dat\n");
        std::exit(1);
    }
    for(i=1; i<=10; i++){
        std::fscanf(fp, "%lf\n", &dE[1][i]);
    }
    std::fclose(fp);

    // 3) dEv_H2O.dat
    fp = std::fopen("inputdata/dEv_H2O.dat", "r");
    if(!fp){
        std::printf("Error: cannot open inputdata/dEv_H2O.dat\n");
        std::exit(1);
    }
    for(i=1; i<=10; i++){
        std::fscanf(fp, "%lf\n", &dE[3][i]);
    }
    std::fclose(fp);

    // 4) sum_dE(m-1)_{i}.dat
    //   若需要跳过 i=2，则在循环里 if(i==2) continue;
    for(i=0; i<=3; i++){
        // 在原代码中 i==2 可能被跳过; 视您需求决定
        if(i == 2) {
            continue; // 如果您确实要跳过
        }
        std::sprintf(filename, "inputdata/sum_dE(m-1)_%d.dat", i);
        fp = std::fopen(filename, "r");
        if(!fp){
            std::printf("Error: cannot open %s\n", filename);
            std::exit(1);
        }
        for(j=0; j<=10; j++){
            std::fscanf(fp, "%lf\n", &E[i][j]);
        }
        std::fclose(fp);
    }

    // 5) c1.dat
    fp = std::fopen("inputdata/c1.dat","r");
    if(!fp){
        std::printf("Error: cannot open c1.dat\n");
        std::exit(1);
    }
    for(i=1; i<=5; i++){
        std::fscanf(fp, "%lf\n", &c1[i]);
    }
    std::fclose(fp);

    // 6) c2.dat
    fp = std::fopen("inputdata/c2.dat","r");
    if(!fp){
        std::printf("Error: cannot open c2.dat\n");
        std::exit(1);
    }
    for(i=1; i<=5; i++){
        std::fscanf(fp, "%lf\n", &c2[i]);
    }
    std::fclose(fp);

    // 7) c3.dat
    fp = std::fopen("inputdata/c3.dat","r");
    if(!fp){
        std::printf("Error: cannot open c3.dat\n");
        std::exit(1);
    }
    for(i=1; i<=5; i++){
        std::fscanf(fp, "%lf\n", &c3[i]);
    }
    std::fclose(fp);
}
