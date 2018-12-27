FC = gfortran
FCFLAGS = -O2
LAPACKLIB = -L/app/compilers/gcc/4.7.1/applib1/LAPACK/3.4.1 -llapack -lrefblas

VPATH = ./ /home/jhkim/work/share

SRCS = kinds.F90 matrix_control.F90 lapack.F90 pseudo_inverse.F90
OBJS = $(SRCS:.F90=.o)
SRC = main_program.F90

.SUFFIXES:
.SUFFIXES:.F90 .o

build: $(SRC) $(OBJS)
	$(FC) $^ $(LAPACKLIB)

clean:
	rm -rf *.mod *.o a.out


.F90.o:
	$(FC) $(FCFLAGS) -c $^
