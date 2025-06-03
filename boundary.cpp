#include "boundary.h"
#include <cstdlib>   // malloc/free
#include <iostream>  // ���������Ҫ�õ� cout�����Ա���
// Ҳ���԰������⼸����ά���亯�����ڱ𴦣�
// ����ֻ���� boundary ��������

// ----------- boundary ����ʵ�� -----------
void negative_boundary_condition(int NR,int NZ,double **ne,int **flag)
{
    for(int i = 0; i < NR; i++){
        for(int j = 0; j < NZ; j++){
            if(flag[i][j]) {
                ne[i][j] = 0.0;
            }
        }
    }
}

void positive_boundary_condition(int num_x,int num_y,double **ion,int **flag)
{
    for(int i = 0; i < num_x; i++){
        for(int j = 0; j < num_y; j++){
            if(flag[i][j]) {
                ion[i][j] = 0.0;
            }
        }
    }
}

void call_flag(int num_x,int num_y,int **flag,int **totuflag,int **otuflag,int **iflag_array,int **jflag_array)
{
    for(int i = 1; i < num_x; i++){
        for(int j = 0; j < num_y; j++){
            // ע��ʵ��ʹ��ʱҪ��Խ��
            if( flag[i-1][j+1] == 1 && flag[i-1][j] == 0 && flag[i][j+1] == 0 ){
                totuflag[i][j] = 1;
            }
            else if( flag[i][j] == 0 && flag[i-1][j] == 1 && flag[i][j+1] == 1 ){
                otuflag[i][j] = 1;
            }
            else if( flag[i][j+1] == 1 && flag[i][j] == 0 ){
                jflag_array[i][j] = 1;
            }
            else if( flag[i-1][j] == 1 && flag[i][j] == 0 ){
                iflag_array[i][j] = 1;
            }
            else {
                totuflag[i][j] = 0;
                otuflag[i][j]   = 0;
                jflag_array[i][j] = 0;
                iflag_array[i][j] = 0;
            }
        }
    }

    // i=0 �߽�
    {
        int i = 0;
        for(int j = 0; j < num_y; j++){
            if(flag[i][j+1] == 1 && flag[i][j] == 0){
                jflag_array[i][j] = 1;
            }
        }
    }
}

void mol_boundary(
    int NR,int NZ,int **flag,
    double **ne,double **N2p,double **O2p, double **H2Op,
    double **O2m,double **Om,double **OHm,double **Hm,
    double **N4p,double **O4p,double **N2O2p,
    double **O2pH2O, double **H3Op, double **H3OpH2O,
    double **H3OpH2O2, double **H3OpH2O3,
    double **Ex,double **Ey,double **absE)
{
    // �ϲ�: j = NZ - 1
    {
        int j = NZ - 1;
        for(int i = 0; i < NR; i++){
            if(!flag[i][j]){
                ne[i][j]   = ne[i][j-1];
                Om[i][j]   = Om[i][j-1];
                O2m[i][j]  = O2m[i][j-1];
                OHm[i][j]  = OHm[i][j-1];
                Hm[i][j]   = Hm[i][j-1];

                N2p[i][j]   = N2p[i][j-1];
                O2p[i][j]   = O2p[i][j-1];
                H2Op[i][j]  = H2Op[i][j-1];
                N4p[i][j]   = N4p[i][j-1];
                O4p[i][j]   = O4p[i][j-1];
                N2O2p[i][j] = N2O2p[i][j-1];
                O2pH2O[i][j] = O2pH2O[i][j-1];
                H3Op[i][j]   = H3Op[i][j-1];
                H3OpH2O[i][j] = H3OpH2O[i][j-1];
                H3OpH2O2[i][j] = H3OpH2O2[i][j-1];
                H3OpH2O3[i][j] = H3OpH2O3[i][j-1];

                absE[i][j] = 0.0;
                Ex[i][j]   = 0.0;
                Ey[i][j]   = 0.0;
            }
        }
    }

    // �Ҳ�: i = NR - 1
    {
        int i = NR - 1;
        for(int j = 0; j < NZ; j++){
            ne[i][j]   = ne[i-1][j];
            Om[i][j]   = Om[i-1][j];
            O2m[i][j]  = O2m[i-1][j];
            OHm[i][j]  = OHm[i-1][j];
            Hm[i][j]   = Hm[i-1][j];

            N2p[i][j]   = N2p[i-1][j];
            O2p[i][j]   = O2p[i-1][j];
            H2Op[i][j]  = H2Op[i-1][j];
            N4p[i][j]   = N4p[i-1][j];
            O4p[i][j]   = O4p[i-1][j];
            N2O2p[i][j] = N2O2p[i-1][j];
            O2pH2O[i][j] = O2pH2O[i-1][j];
            H3Op[i][j]   = H3Op[i-1][j];
            H3OpH2O[i][j] = H3OpH2O[i-1][j];
            H3OpH2O2[i][j] = H3OpH2O2[i-1][j];
            H3OpH2O3[i][j] = H3OpH2O3[i-1][j];
        }
    }

    // �����Ҫ���� j=0��i=0 �������߽磬�����в��䡣
}

