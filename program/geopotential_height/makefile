FC = gfortran
FCFLAGS = -O2 -ffree-line-length-none -fcheck=all

LAPACKLIB = -L/app/compilers/gcc/4.7.1/applib1/LAPACK/3.4.1 -llapack -lrefblas

VPATH = ./ /home/jhkim/work/share

SRCS = kinds.F90 constant.F90 std_atmosphere.F90 logger.F90 statistics.F90 matrix_control.F90 std_atmosphere.F90 read_hybrid_file.F90 cubic_spline_interp.F90 spline_interp.F90 hybrid_coord.F90 geopotential.F90
OBJS = $(SRCS:.F90=.o)

SRC  = main_program.F90
EXEC = run.exe


.SUFFIXES:
.SUFFIXES:.F90 .o

exe: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) -o $(EXEC) $^ $(LAPACKLIB)

clean:
	rm -rf *.mod *.o a.out *.exe

.F90.o:
	$(FC) $(FCFLAGS) -c $< $(LAPACKLIB)
