#ifndef CALC_E_H
#define CALC_E_H

// ��������ڲ���Ҫ�õ�һЩ��׼�������������������� include <cmath>, <cstdlib> �ȣ�
// ����ͷ�ļ�ֻ�������Ļ�������ֻд��������������ʵ���� calcE.cpp

// ����һ������/�ͷŶ�ά����ĺ���(����Ҫ), Ҳ���ڱ�ͳһ����
//double** allocate2DDouble(int rows, int cols);
//void free2DDouble(double** arr, int rows);

// ���� calc_E(...) ����
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
            int Prinstp);

#endif // CALC_E_H

