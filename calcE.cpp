#include "calcE.h"    // ����ͷ�ļ�����֤����ƥ��
#include <cmath>      // sqrt, pow
#include <cstdlib>    // malloc, free
#include <iostream>   // �����Ҫ�ں������õ� std::cout

// ���û��ȫ�ֵĶ�ά���亯������������Զ���
static double** allocate2DDouble(int rows, int cols) {
    double** arr = (double**)std::malloc(rows * sizeof(double*));
    for(int i = 0; i < rows; i++){
        arr[i] = (double*)std::malloc(cols * sizeof(double));
        for(int j = 0; j < cols; j++){
            arr[i][j] = 0.0;
        }
    }
    return arr;
}
static void free2DDouble(double** arr, int rows) {
    if(!arr) return;
    for(int i = 0; i < rows; i++){
        if(arr[i]) std::free(arr[i]);
    }
    std::free(arr);
}

//================= calc_E ������ʵ�� ===================//
void calc_E(int NR,int NZ,
            double **phi,   // ��λ
            double **absE,  // �糡��С
            double **Ey,    // y����糡
            double **Ex,    // x����糡
            double air_kg,
            int **totuflag,
            int **otuflag,
            int **iflag,
            int **jflag,
            int **flag,
            double a,
            double b,
            double *rhalf,
            double *zhalf,
            double **Mol,
            int num,
            int Prinstp)
{
    int i, j;
    double Avo = 6.022141e+23;  // ����٤���޳���

    // ���ԭ�� mat(NR,NZ)
    double** LEx = allocate2DDouble(NR, NZ);
    double** mol = allocate2DDouble(NR, NZ);

    // �� kg/m^3 ת��Ϊ (������)/m^3
    for(i=0; i<NR; i++){
        for(j=0; j<NZ; j++){
            mol[i][j] = Mol[i][j] / (air_kg / Avo);
        }
    }

    // ���� Ex �� LEx
    for(i=0; i<NR; i++){
        for(j=0; j<NZ; j++){
            LEx[i][j] = Ex[i][j];
        }
    }

    // ��ѭ�������� Ex, Ey
    for(i = 1; i < NR; i++){
        for(j = 1; j < NZ; j++){
            if(jflag[i][j-1] == 1) {
                double hzm = b * std::sqrt(1.0 + std::pow(rhalf[i]/a, 2)) - zhalf[j-1];
                Ex[i][j] = -1e+21 * (phi[i][j] - phi[i-1][j]) 
                           / ((rhalf[i] - rhalf[i-1]) * mol[i][j]);
                Ey[i][j] = -1e+21 * (phi[i][j] - phi[i][j-1])
                           / (hzm * mol[i][j]);

            } else if(iflag[i][j] == 1 && j != NZ-1) {
                double hrm = rhalf[i] - a * std::sqrt(std::pow(zhalf[j]/b,2) - 1.0);
                Ex[i][j] = -1e+21 * (phi[i][j] - phi[i-1][j]) 
                           / (hrm * mol[i][j]);
                Ey[i][j] = -1e+21 * (phi[i][j] - phi[i][j-1])
                           / ((zhalf[j] - zhalf[j-1]) * mol[i][j]);

            } else if(otuflag[i][j] == 1) {
                double hrm = rhalf[i] - a * std::sqrt(std::pow(zhalf[j]/b,2) - 1.0);
                Ex[i][j] = -1e+21 * (phi[i][j] - phi[i-1][j]) 
                           / (hrm * mol[i][j]);
                Ey[i][j] = -1e+21 * (phi[i][j] - phi[i][j-1]) 
                           / ((zhalf[j] - zhalf[j-1]) * mol[i][j]);

            } else {
                // ��ͨ���
                Ex[i][j] = -1e+21 * (phi[i][j] - phi[i-1][j]) 
                           / ((rhalf[i] - rhalf[i-1]) * mol[i][j]);
                Ey[i][j] = -1e+21 * (phi[i][j] - phi[i][j-1]) 
                           / ((zhalf[j] - zhalf[j-1]) * mol[i][j]);
            }

            // otuflag[i][j-1] �ٴ����� Ey
            if(otuflag[i][j-1]) {
                double hzm = b * std::sqrt(1.0 + std::pow(rhalf[i]/a,2)) - zhalf[j-1];
                Ey[i][j] = -1e+21 * (phi[i][j] - phi[i][j-1]) 
                           / (hzm * mol[i][j]);
            }
        }
    }

    // i=0
    {
        i=0;
        for(j=1; j<NZ; j++){
            if(flag[i][j] == 0){
                Ey[i][j] = -1e+21*(phi[i][j] - phi[i][j-1]) 
                           / ((zhalf[j] - zhalf[j-1]) * mol[i][j]);
            }
            if(jflag[i][j-1] == 1){
                double hzm = b * std::sqrt(1.0 + std::pow(rhalf[i]/a,2)) - zhalf[j-1];
                Ey[i][j] = -1e+21*(phi[i][j] - phi[i][j-1]) 
                           / (hzm * mol[i][j]);
            }
            Ex[i][j] = 0.0;
        }
    }

    // j=0
    {
        int j=0;
        for(i=1; i<NR; i++){
            Ex[i][j] = -1e+21*(phi[i][j] - phi[i-1][j]) 
                       / ((rhalf[i] - rhalf[i-1]) * mol[i][j]);
            Ey[i][j] = -1e+21*(phi[i][j] - 0.0)
                       / ((zhalf[j] - 0.0) * mol[i][j]);
        }
    }

    // i=0, j=0
    {
        i=0; int j=0;
        Ex[i][j] = 0.0;
        Ey[i][j] = -1e+21*(phi[i][j] - 0.0)
                   / ((zhalf[j] - 0.0) * mol[i][j]);
    }

    // �� Ex ����������
    for(i=0; i<NR; i++){
        for(j=0; j<NZ; j++){
            double zmm = zhalf[j]*1000 - 12.7;
            double sigm = (2.0*zmm/(1+std::fabs(2.0*zmm)) + 1.0)*0.5;
            Ex[i][j] = Ex[i][j] - 0.5 * sigm * LEx[i][j];
        }
    }

    // ���� absE
    double Ezlim = -30.0;
    for(i=0; i<NR-1; i++){
        for(j=0; j<NZ-1; j++){
            double aExEy = std::sqrt(Ex[i][j]*Ex[i][j] + Ey[i][j]*Ey[i][j]);
            if(aExEy == 0.0){
                absE[i][j] = 0.0;
            } else {
                double aEx = (Ex[i][j] + Ex[i+1][j]) * 0.5;
                double aEy;
                if(Ey[i][j+1] == 0.0){
                    aEy = Ey[i][j];
                } else {
                    if(Ey[i][j] < Ey[i][j+1] + Ezlim && num < Prinstp + 10000){
                        aEy = Ey[i][j];
                    } else {
                        aEy = (Ey[i][j] + Ey[i][j+1]) * 0.5;
                    }
                }
                absE[i][j] = std::sqrt(aEx*aEx + aEy*aEy);
            }
        }
    }

    // �ұ� i=NR-1
    {
        i=NR-1;
        for(int j=0; j<NZ-1; j++){
            double aEx = Ex[i][j];
            double aEy = (Ey[i][j] + Ey[i][j+1])*0.5;
            absE[i][j] = std::sqrt(aEx*aEx + aEy*aEy);
        }
    }
    // ���� j=NZ-1
    {
        int j=NZ-1;
        for(i=0; i<NR-1; i++){
            double aEx = (Ex[i][j] + Ex[i+1][j])*0.5;
            double aEy = Ey[i][j];
            absE[i][j] = std::sqrt(aEx*aEx + aEy*aEy);
        }
    }
    // ���Ͻ�
    {
        i=NR-1; int j=NZ-1;
        double aEx = Ex[i][j];
        double aEy = Ey[i][j];
        absE[i][j] = std::sqrt(aEx*aEx + aEy*aEy);
    }

    // �������Ҫ�ü�����糡�����ڴ˼��߼�

    // �ͷ���ʱ��ά����
    free2DDouble(mol, NR);
    free2DDouble(LEx, NR);
}


//================= ��ѡ��һ���򵥵Ĳ��� main() =================//
// �粻��Ҫ����ɾ����ע�͵�
#ifdef TEST_CALCE_MAIN
int main()
{
    std::cout << "=== Testing calc_E ===\n";

    int NR=3, NZ=3;
    // ��������� phi, absE, Ex, Ey, Mol ��
    // ���� totuflag, otuflag, iflag, jflag, flag ���򵥳�ʼ��
    // ���� calc_E(...), Ȼ������鿴

    std::cout << "Done.\n";
    return 0;
}
#endif

