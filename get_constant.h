#ifndef GET_CONSTANT_H
#define GET_CONSTANT_H

#ifdef __cplusplus
extern "C" {
#endif

// 如果您使用的是 CUDA .cu + nvcc, 可能需要 extern "C"。
// 如果您使用纯 C++，不需要 extern "C"。可自行去掉。
// 如果您只是把 main.cu 改成 .cpp 也可直接写普通声明。

// 声明：从文件读取常量/数组等，并初始化 spline
void get_constant(void);

#ifdef __cplusplus
}
#endif

#endif // GET_CONSTANT_H
