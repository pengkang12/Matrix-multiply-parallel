# Scalepack-matrix
# 1. using mkl
# 2.download mkl library from Intel

# 3. install mkl library to your computer

source  /opt/intel/compilers_and_libraries_2020/linux/mkl/bin/mklvars.sh intel64

pay attention: your need to modify address

# 4. compiler and link mkl library  

using this link to check your link and compile parameter.  
https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

you can use mpicc or gcc
```
export MKLROOT=/home/peng/intel/mkl

mpicc matrix-multiplication.c  -DMKL_ILP64 -m64 -I${MKLROOT}/include -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
```

# 5. Runing program 
```
mpicc -n 1 ./a.out
```
