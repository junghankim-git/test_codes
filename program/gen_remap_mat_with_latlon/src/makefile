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


# NetCDF
NETCDF_DIR     = /app/compilers/gcc/4.7.1/applib1/NETCDF4/4.1.3
NETCDF_INC_DIR = $(NETCDF_DIR)/include
NETCDF_LIB_DIR = $(NETCDF_DIR)/lib
LAPACK_DIR     = /app/compilers/gcc/4.7.1/applib1/LAPACK/3.4.1

EXTERN_INC  = -I$(NETCDF_INC_DIR)
EXTERN_LIB  = -L$(NETCDF_LIB_DIR) -lnetcdff -lnetcdf -L$(LAPACK_DIR) -llapack -lrefblas


VPATH = /home/jhkim/work/share ./src

# Auxiliary
SRCS = kinds.F90 logger.F90 buffer.F90 mesh_control.F90 quadrature.F90 coordinates.F90 \
       space_filling_curve.F90 statistics.F90 matrix_control.F90 bilinear.F90          \
       scrip_input.F90 cubed_sphere.F90 remap_matrix.F90 remap_netcdf.F90 remapper.F90
OBJS = $(SRCS:.F90=.o)
MODS = $(shell echo $(SRCS:.F90=.mod) | tr A-Z a-z)

# Main
SRC  = main_program.F90
EXE  = run.exe

.SUFFIXES:
.SUFFIXES: .F90 .o


build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) -o $(EXE) $^ $(EXTERN_INC) $(EXTERN_LIB)

check: scrip_check.F90 $(OBJS)
	$(FC) $(FCFLAGS) -o check.exe $^ $(EXTERN_INC) $(EXTERN_LIB)

clean:
	rm -rf $(OBJS) $(MODS) $(EXE)

.F90.o:
	$(FC) $(FCFLAGS) -c $< -o $@ $(EXTERN_INC) $(EXTERN_LIB)
