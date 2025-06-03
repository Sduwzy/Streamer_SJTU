#include "memory.h"
#include <cstdio>   // 如果函数内要用 printf
#include <cstdlib>  // malloc/free
#include <cmath>    // 如果有用到 math 函数

//------------------- 分配函数 -------------------//

double* vec(int n) {
    // 分配 n+1 个 double
    double* v = (double *)std::malloc((n+1) * sizeof(double));
    if(!v) {
        std::printf("allocation error in vector()\n");
        // 或者可以考虑throw异常或 std::abort()
    }
    return v;
}

float* fvec(int n){
    float* v = (float*)std::malloc((n+1)*sizeof(float));
    if(!v) {
        std::printf("allocation error in fvec()\n");
    }
    return v;
}

int* ivec(int n){
    int* v = (int*)std::malloc((n+1)*sizeof(int));
    if(!v) {
        std::printf("allocation error in ivec()\n");
    }
    return v;
}

// mat 分配 n1+1 行, 每行 n2+1 列
double** mat(int n1, int n2) {
    double** m = (double**)std::malloc((n1+1)*sizeof(double*));
    if(!m) {
        std::printf("allocation error1 in matrix()\n");
    }
    // m[0] 分配 n1*(n2+1) 个 double
    m[0] = (double*)std::malloc(n1*(n2+1)*sizeof(double));
    if(!m[0]) {
        std::printf("allocation error2 in matrix()\n");
    }
    // 后续行指针指向 m[0] 中的相应位置
    for(int i=0; i<n1; i++){
        m[i+1] = m[i] + (n2+1);
    }
    return m;
}

float** fmat(int n1, int n2){
    float** m = (float**)std::malloc((n1+1)*sizeof(float*));
    if(!m) {
        std::printf("allocation error1 in fmat()\n");
    }
    m[0] = (float*)std::malloc(n1*(n2+1)*sizeof(float));
    if(!m[0]) {
        std::printf("allocation error2 in fmat()\n");
    }
    for(int i=0; i<n1; i++){
        m[i+1] = m[i] + (n2+1);
    }
    return m;
}

int** imat(int n1, int n2){
    int** m = (int**)std::malloc((n1+1)*sizeof(int*));
    if(!m) {
        std::printf("allocation error1 in imat()\n");
    }
    m[0] = (int*)std::malloc(n1*(n2+1)*sizeof(int));
    if(!m[0]) {
        std::printf("allocation error2 in imat()\n");
    }
    for(int i=0; i<n1; i++){
        m[i+1] = m[i] + (n2+1);
    }
    return m;
}

char** cmat(int n1, int n2){
    char** m = (char**)std::malloc((n1+1)*sizeof(char*));
    if(!m) {
        std::printf("allocation error1 in cmat()\n");
    }
    m[0] = (char*)std::malloc(n1*(n2+1)*sizeof(char));
    if(!m[0]) {
        std::printf("allocation error2 in cmat()\n");
    }
    for(int i=0; i<n1; i++){
        m[i+1] = m[i] + (n2+1);
    }
    return m;
}

char* cvec(int n){
    char* v = (char*)std::malloc((n+1)*sizeof(char));
    if(!v) {
        std::printf("allocation error in cvec()\n");
    }
    return v;
}


//------------------- 释放函数 -------------------//

void free_vec(double* v, int n){
    (void)n;  // 如果不需要用 n，也可以直接不写
    std::free(v);
}

void free_mat(double** m, int n1, int n2){
    (void)n1; 
    (void)n2; 
    // m[0] 分配了 n1*(n2+1)
    std::free(m[0]);
    std::free(m);
}

void free_cvec(char* v, int n){
    (void)n; 
    std::free(v);
}

void free_ivec(int* v, int n){
    (void)n; 
    std::free(v);
}

// 如果您还需要 free_fmat, free_imat, free_cmat等，可以类似写

