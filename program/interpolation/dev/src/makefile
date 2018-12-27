FC = gfortran
FCFLAGS = -ffree-line-length-none -fcheck=all

LAPACKLIB = -L/app/compilers/gcc/4.7.1/applib1/LAPACK/3.4.1 -llapack -lrefblas
NETCDF = /app/compilers/gcc/4.7.1/applib1/NETCDF4/4.1.3
NETCDFLIB = -L$(NETCDF)/lib -lnetcdff -lnetcdf -I$(NETCDF)/include

VPATH = /home/jhkim/work/share

SRCS = kinds.F90 logger.F90 io_ascii.F90 quadrature.F90 gen_test_func.F90 \
       newton_interp.F90 lagrange_interp.F90 spline_interp.F90 io_module.F90
OBJS = $(SRCS:.F90=.o)
SRC  = main_test.F90
EXE  = run.exe

.SUFFIXES:
.SUFFIXES:.F90 .o

build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) $^ -o $(EXE) $(LAPACKLIB) $(NETCDFLIB)
	cp $(EXE) ../

clean:
	rm -rf *.mod *.o a.out $(EXE)

.F90.o:
	$(FC) $(FCFLAGS) -c $< $(LAPACKLIB) $(NETCDFLIB)
