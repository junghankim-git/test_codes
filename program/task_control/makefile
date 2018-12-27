FC      = mpif90
FCFLAGS = -O2 -ffree-line-length-none -fcheck=all
#FC      = ifort
#FCFLAGS = -O2 -free -check=all

VPATH = /home/jhkim/work/share
INC_PATH = -I$(VPATH)

SRCS = kinds.F90 multi_io.F90 task_control.F90
OBJS = $(SRCS:.F90=.o)
MODS = $(SRCS:.F90=.mod)

SRC  = main_program.F90
EXE  = run.exe

NETCDF_DIR = /app/compilers/gcc/4.7.1/applib1/NETCDF4/4.1.3
NETCDF_INC = $(NETCDF_DIR)/include
NETCDF_LIB = $(NETCDF_DIR)/lib

EXT_LIBS   = -I$(NETCDF_INC) -L$(NETCDF_LIB) -lnetcdff -lnetcdf

.SUFFIXES:
.SUFFIXES:.F90 .o


build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) $^ -o $(EXE) $(INC_PATH) $(EXT_LIBS)

clean:
	rm -f $(OBJS) $(MODS) $(EXE) *.nc

.F90.o:
	$(FC) $(FCFLAGS) -c $< $(INC_PATH) $(EXT_LIBS)
