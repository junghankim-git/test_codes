COMPILER = gnu

ifeq ($(COMPILER),gnu)
FC      = mpif90
FCFLAGS = -O2 -ffree-line-length-none -fcheck=all
endif
ifeq ($(COMPILER),intel)
FC      = mpif90
FCFLAGS = -O2 -free -check all
endif

VPATH = /home/jhkim/work/share
INC_PATH = -I$(VPATH)

SRCS = kinds.F90 lagrangian_particle.F90
OBJS = $(SRCS:.F90=.o)
MODS = $(SRCS:.F90=.mod)

SRC  = main_program.F90
EXE  = run.exe

ifeq ($(COMPILER),gnu)
NETCDF_DIR = /app/compilers/gcc/4.7.1/applib1/NETCDF4/4.1.3
PNETCDF_DIR = /app/compilers/gcc/4.7.1/mpi/mvapich2/1.8.1/applib2/PNETCDF/1.6.1
endif
ifeq ($(COMPILER),intel)
NETCDF_DIR = /app/compilers/intel/2013.0.028/applib1/NETCDF4/4.1.3
PNETCDF_DIR = /app/compilers/intel/2013.0.028/mpi/mvapich2/1.8.1/applib2/PNETCDF/1.6.1
endif
NETCDF_INC = $(NETCDF_DIR)/include
NETCDF_LIB = $(NETCDF_DIR)/lib
PNETCDF_INC = $(PNETCDF_DIR)/include
PNETCDF_LIB = $(PNETCDF_DIR)/lib

EXT_LIBS   = -I$(NETCDF_INC) -L$(NETCDF_LIB) -lnetcdff -lnetcdf -I$(PNETCDF_INC) -L$(PNETCDF_LIB) -lpnetcdf

.SUFFIXES:
.SUFFIXES:.F90 .o


build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) $^ -o $(EXE) $(INC_PATH) $(EXT_LIBS)
	cp $(EXE) ../

clean:
	rm -f $(OBJS) $(MODS) $(EXE) *.nc
	rm -f ../$(EXE)

.F90.o:
	$(FC) $(FCFLAGS) -c $< $(INC_PATH) $(EXT_LIBS)
