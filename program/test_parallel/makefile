FC = mpif90
FCFLAGS = -ffree-line-length-none

LAPACKLIB = -L/app/compilers/gcc/4.7.1/applib1/LAPACK/3.4.1 -llapack -lrefblas

VPATH = /home/jhkim/work/share

SRCS = kinds.F90 parallel.F90 par_env.F90 par_mpi.F90 par_func.F90 repro_sum.F90
OBJS = $(SRCS:.F90=.o)
SRC = main_test.F90

.SUFFIXES:
.SUFFIXES:.F90 .o

build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) $^ $(LAPACKLIB)

clean:
	rm -rf *.mod *.o a.out

.F90.o:
	$(FC) $(FCFLAGS) -c $< $(LAPACKLIB)
