FC = gfortran
FFLAGS = -Wall -O2 -Wno-unused-variable -ffree-line-length-256  -fmax-errors=1 -Wno-uninitialized -fopenmp
SERIAL_FLAGS = -Wall -O2 -Wno-unused-variable -ffree-line-length-256  -fmax-errors=1 -Wno-uninitialized
FILES = global_var.f90 vector.f90 file_io.f90 initialize.f90 force.f90 integrate.f90 main.f90
# 如果希望进行单步调试, 应该使用-O0强制编译器不优化代码,
# 在FFLAGS中加入-cpp可以启用GNU编译器中的define宏定义

target:
	${FC} ${FILES} -o $@ ${FFLAGS} 

serial:
	${FC} ${FILES} -o target_serial ${SERIAL_FLAGS}

.PHONY: clean

clean:
	rm *.mod  ../xyz/trajectory_*