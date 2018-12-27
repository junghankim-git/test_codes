FC = mpif90
#FCFLAGS = -O2 -ffree-line-length-none
FCFLAGS = -O2 -ffree-line-length-none -fcheck=all

LAPACKLIB = -L/app/compilers/gcc/4.7.1/applib1/LAPACK/3.4.1 -llapack -lrefblas

VPATH = /home/jhkim/work/share

SRCS = kinds.F90 constant.F90 logger.F90 io_ascii.F90 statistics.F90 matrix_control.F90 quadrature.F90 \
       gen_test_func.F90 lagrange_interp.F90 piecewise_remap.F90 spline_remap.F90
OBJS = $(SRCS:.F90=.o)
SRC = main_test.F90
EXE = run.exe

.SUFFIXES:
.SUFFIXES:.F90 .o

build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) $^ -o $(EXE) $(LAPACKLIB)

clean:
	rm -rf *.mod *.o a.out $(EXE)

.F90.o:
	$(FC) $(FCFLAGS) -c $< $(LAPACKLIB)
