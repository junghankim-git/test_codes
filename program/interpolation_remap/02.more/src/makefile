FC = mpif90
#FCFLAGS = -O2 -ffree-line-length-none
FCFLAGS = -O2 -ffree-line-length-none -fcheck=all

NETCDF     = /usr/lib/x86_64-linux-gnu
NETCDFLIB  = -L$(NETCDF) -lnetcdff -lnetcdf -I$(NETCDF)/../../include
LAPACK_DIR = /usr/lib/x86_64-linux-gnu/lapack
BLAS_DIR   = /usr/lib/x86_64-linux-gnu/blas
LAPACKLIB  =  -L$(LAPACK_DIR) -llapack -L$(BLAS_DIR) -lblas

VPATH = ../../../../share

SRCS = kinds.F90 constant.F90 statistics.F90 gen_test_func.F90 \
       lagrange_interp.F90 piecewise_remap.F90 spline_remap.F90 io_module.F90
OBJS = $(SRCS:.F90=.o)
SRC = main_test.F90
EXE = run.exe

.SUFFIXES:
.SUFFIXES:.F90 .o

build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) $^ -o $(EXE) $(LAPACKLIB) $(NETCDFLIB)
	cp $(EXE) ../

clean:
	rm -rf *.mod *.o $(EXE)

.F90.o:
	$(FC) $(FCFLAGS) -c $< $(LAPACKLIB) $(NETCDFLIB)
