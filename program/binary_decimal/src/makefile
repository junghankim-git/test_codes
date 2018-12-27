#------------------------------------
# Compiler Environment
#------------------------------------
FC    = gfortran
MPIFC = mpif90
CPP   = cpp
AR    = ar
MV    = mv
CP    = cp
RM    = rm

FCFLAGS := -O2 -ffree-line-length-none -fcheck=all 

NETCDF = /app/compilers/gcc/4.7.1/applib1/NETCDF4/4.1.3
NETCDF_INC = $(NETCDF)/include
NETCDF_LIB = $(NETCDF)/lib
EXTERN_INC = -I$(NETCDF_INC)
EXTERN_LIB = -L$(NETCDF_LIB) -lnetcdff -lnetcdf

VPATH = ./src /home/jhkim/work/share

# Auxiliary
SRCS = kinds.F90 binary_decimal.F90
OBJS = $(SRCS:.F90=.o)
MODS = $(shell echo $(SRCS:.F90=.mod) | tr A-Z a-z)

# Main
SRC  = main_program.F90
EXE  = run.exe

.SUFFIXES:
.SUFFIXES: .F90 .o


build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) -o $(EXE) $^ $(EXTERN_INC) $(EXTERN_LIB)
	cp run.exe ../

clean:
	rm -rf $(OBJS) $(MODS) $(EXE) *.exe ../*.exe

.F90.o:
	$(FC) $(FCFLAGS) -c $< -o $@ $(EXTERN_INC) $(EXTERN_LIB)
