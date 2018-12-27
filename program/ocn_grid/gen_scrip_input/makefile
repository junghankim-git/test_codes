# Compiler
FC=gfortran
FCFLAGS=-O2 -ffree-line-length-none -fcheck=all

# External library
NETCDF_DIR=/app/compilers/gcc/4.7.1/applib1/NETCDF4/4.1.3
NETCDF_INC=$(NETCDF_DIR)/include
NETCDF_LIB=$(NETCDF_DIR)/lib
EXTERNAL=-L$(NETCDF_LIB) -lnetcdff -lnetcdf -I$(NETCDF_INC)


SRCS=kinds.F90 scrip_input.F90 scrip_hycom.F90
OBJS=$(SRCS:.F90=.o)
MODS=$(shell echo $(OBJS:.o=.mod) | tr A-Z a-z)

SRC = main_program.F90
EXE = run.exe

build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) $^ -o $(EXE) $(EXTERNAL)

clean:
	rm -f $(OBJS) $(MODS) $(EXE)

.SUFFIXES:
.SUFFIXES:.F90 .o

.F90.o:
	$(FC) $(FCFLAGS) -c $< $(EXTERNAL)
