#ifndef MEMORY_H
#define MEMORY_H

#include <cstddef> // size_t
#include <cstdio>  // printf
#include <cstdlib> // malloc/free

// 各种分配函数声明
double* vec(int n);
float* fvec(int n);
int* ivec(int n);

double** mat(int n1, int n2);
float** fmat(int n1, int n2);
int** imat(int n1, int n2);
char** cmat(int n1, int n2);

char* cvec(int n);

// 各种释放函数声明
void free_vec(double* v, int n);
void free_mat(double** m, int n1, int n2);
void free_cvec(char* v, int n);
void free_ivec(int* v, int n);

// 如果您需要把 fvec/fmat/cmat 等也添加对应 free_fvec / free_fmat 等，可以再写

#endif // MEMORY_H

