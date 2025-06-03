#ifndef MEMORY_H
#define MEMORY_H

#include <cstddef> // size_t
#include <cstdio>  // printf
#include <cstdlib> // malloc/free

// ���ַ��亯������
double* vec(int n);
float* fvec(int n);
int* ivec(int n);

double** mat(int n1, int n2);
float** fmat(int n1, int n2);
int** imat(int n1, int n2);
char** cmat(int n1, int n2);

char* cvec(int n);

// �����ͷź�������
void free_vec(double* v, int n);
void free_mat(double** m, int n1, int n2);
void free_cvec(char* v, int n);
void free_ivec(int* v, int n);

// �������Ҫ�� fvec/fmat/cmat ��Ҳ��Ӷ�Ӧ free_fvec / free_fmat �ȣ�������д

#endif // MEMORY_H

