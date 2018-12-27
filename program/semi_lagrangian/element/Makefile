FC = gfortran
FCFLAGS = -ffree-line-length-none -fcheck=all

LAPACKLIB = -L/app/compilers/gcc/4.7.1/applib1/LAPACK/3.4.1 -llapack -lrefblas

VPATH = /home/jhkim/work/share

SRCS = kinds.F90 logger.F90 io_ascii.F90 controls.F90 quadrature.F90 element.F90  \
       lagrange_interp.F90 piecewise_remap.F90 spline_remap.F90 semi_lagrangian.F90
OBJS = $(SRCS:.F90=.o)
SRC = main_test.F90
EXE = run.exe

.SUFFIXES:
.SUFFIXES:.F90 .o

build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) $^ -o $(EXE) $(LAPACKLIB)

clean:
	rm -rf *.mod *.o $(EXE)

.F90.o:
	$(FC) $(FCFLAGS) -c $< $(LAPACKLIB)
