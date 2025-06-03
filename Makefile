##############################################################################
# Makefile (使用 NVCC 编译链接所有源文件的示例)
##############################################################################

# 配置 NVCC 及编译选项
NVCC      = /usr/local/cuda/bin/nvcc
#NVCCFLAGS = -O2 -std=c++11
NVCCFLAGS += -allow-unsupported-compiler -Xcompiler -w

# 默认目标
all: main

# --------------------- 各个 .o 的生成规则 ---------------------
# 无论文件里有没有用到CUDA API, 用NVCC编译也可正常处理
boundary.o: boundary.cpp boundary.h
	$(NVCC) $(NVCCFLAGS) -c boundary.cpp -o boundary.o

Reaction.o: Reaction.cpp Reaction.h
	$(CXX) $(CXXFLAGS) -c Reaction.cpp -o Reaction.o

poisson.o: poisson.cpp poisson.h memory.h
	$(CXX) $(CXXFLAGS) -c poisson.cpp -o poisson.o

first_q.o: first_q.cpp first_q.h
	$(CXX) $(CXXFLAGS) -c first_q.cpp -o first_q.o

Energy_type.o: Energy_type.cpp Energy_type.h
	$(CXX) $(CXXFLAGS) -c Energy_type.cpp -o Energy_type.o

memory.o: memory.cpp memory.h
	$(NVCC) $(NVCCFLAGS) -c memory.cpp -o memory.o

spline.o: spline.cpp spline.h memory.h
	$(NVCC) $(NVCCFLAGS) -c spline.cpp -o spline.o

mesh_generator.o: mesh_generator.cpp mesh_generator.h
	$(NVCC) $(NVCCFLAGS) -c mesh_generator.cpp -o mesh_generator.o

calcE.o: calcE.cpp calcE.h
	$(NVCC) $(NVCCFLAGS) -c calcE.cpp -o calcE.o

calc_velo.o: calc_velo.cpp calc_velo.h
	$(NVCC) $(NVCCFLAGS) -c calc_velo.cpp -o calc_velo.o

Red_Black_SOR_Kernel.o: Red_Black_SOR_Kernel.cu Red_Black_SOR_Kernel.h
	$(NVCC) $(NVCCFLAGS) -c Red_Black_SOR_Kernel.cu -o Red_Black_SOR_Kernel.o

superbee.o: superbee.cpp superbee.h
	$(CXX) $(CXXFLAGS) -c superbee.cpp -o superbee.o

Helmholtz_Kernel.o: Helmholtz_Kernel.cu Helmholtz_Kernel.h
	$(NVCC) $(NVCC_FLAGS) -c Helmholtz_Kernel.cu -o Helmholtz_Kernel.o


Multigrid_kernel.o: Multigrid_kernel.cu Multigrid_kernel.h
	$(NVCC) $(NVCCFLAGS) -c Multigrid_kernel.cu -o Multigrid_kernel.o

main.o: main.cpp memory.h spline.h mesh_generator.h calcE.h boundary.h Energy_type.h Reaction.h first_q.h poisson.h Red_Black_SOR_Kernel.h Helmholtz_Kernel.h Multigrid_kernel.h calc_velo.h
	$(NVCC) $(NVCCFLAGS) -c main.cpp -o main.o

# --------------------- 最终的链接阶段 ---------------------
# 使用 NVCC 来进行链接, 自动链接到 CUDA runtime

main: memory.o spline.o mesh_generator.o calcE.o main.o boundary.o Energy_type.o Reaction.o first_q.o poisson.o Red_Black_SOR_Kernel.o Helmholtz_Kernel.o Multigrid_kernel.o calc_velo.o superbee.o
	$(NVCC) $(NVCCFLAGS) memory.o spline.o mesh_generator.o calcE.o boundary.o Energy_type.o Reaction.o first_q.o poisson.o Red_Black_SOR_Kernel.o Helmholtz_Kernel.o Multigrid_kernel.o calc_velo.o superbee.o main.o -o main

# --------------------- 清理规则 ---------------------

clean:
	rm -f *.o main
