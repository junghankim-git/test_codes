FC      = mpif90
FCFLAGS = -O2 -ffree-line-length-none -fcheck=all -fbacktrace
#FCFLAGS = -O2 -ffree-line-length-none
#FCFLAGS = -O2 -mcmodel=large

VPATH = /home/jhkim/work/share

SRCS = kinds.F90 statistics.F90 matrix_control.F90
OBJS = $(SRCS:.F90=.o)
MODS = $(SRCS:.F90=.mod)

SRC = main_program.F90
OBJ = $(SRC:.F90=.o)

EXE = run.exe

LAPACK   = /app/compilers/gcc/4.7.1/applib1/LAPACK/3.4.1
#LAPACK   = /app/compilers/intel/2013.0.028/applib1/LAPACK
NETCDF   = /app/compilers/gcc/4.7.1/applib1/NETCDF4/4.1.3
EXT_LIB  = -L$(NETCDF)/lib -lnetcdff -lnetcdf -I$(NETCDF)/include -L$(LAPACK) -llapack -lrefblas

.SUFFIXES:
.SUFFIXES:.F90 .o

build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) -o $(EXE) $^ $(EXT_LIB)
	cp $(EXE) ../

clean:
	rm -rf $(OBJS) $(MODS) $(EXE)

run:
	mpirun -n 4 ./$(EXE)


.F90.o:
	$(FC) $(FCFLAGS) -c $< $(EXT_LIB)
