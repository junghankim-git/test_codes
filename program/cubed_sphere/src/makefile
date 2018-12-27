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
NETCDF_DIR     = /usr/lib/x86_64-linux-gnu
NETCDF_INC_DIR = $(NETCDF_DIR)/../../include
NETCDF_LIB_DIR = $(NETCDF_DIR)
LAPACK_DIR     = /usr/lib/x86_64-linux-gnu/lapack
BLAS_DIR       = /usr/lib/x86_64-linux-gnu/blas

EXTERN_INC  = -I$(NETCDF_INC_DIR)
EXTERN_LIB  = -L$(NETCDF_LIB_DIR) -lnetcdff -lnetcdf -L$(LAPACK_DIR) -llapack -L$(BLAS_DIR) -lblas


VPATH = ./src ../../../share

# Auxiliary
SRCS = kinds.F90 logger.F90 mesh_control.F90 quadrature.F90 coordinates.F90 space_filling_curve.F90 statistics.F90 matrix_control.F90 bilinear.F90 scrip_input.F90 cubed_sphere.F90
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

rotation1: test_rotation1.F90 $(OBJS)
	$(FC) $(FCFLAGS) -o rotation1.exe $^ $(EXTERN_INC) $(EXTERN_LIB)
	cp rotation1.exe ../

rotation2: test_rotation2.F90 $(OBJS)
	$(FC) $(FCFLAGS) -o rotation2.exe $^ $(EXTERN_INC) $(EXTERN_LIB)
	cp rotation2.exe ../

rotation3: test_rotation3.F90 $(OBJS)
	$(FC) $(FCFLAGS) -o rotation3.exe $^ $(EXTERN_INC) $(EXTERN_LIB)
	cp rotation3.exe ../

rotation4: test_rotation4.F90 $(OBJS)
	$(FC) $(FCFLAGS) -o rotation4.exe $^ $(EXTERN_INC) $(EXTERN_LIB)
	cp rotation4.exe ../

inverse: test_inverse.F90 $(OBJS)
	$(FC) $(FCFLAGS) -o inverse.exe $^ $(EXTERN_INC) $(EXTERN_LIB)
	cp inverse.exe ../

convert: test_convert_coord.F90 $(OBJS)
	$(FC) $(FCFLAGS) -o convert.exe $^ $(EXTERN_INC) $(EXTERN_LIB)
	cp convert.exe ../

clean:
	rm -rf $(OBJS) $(MODS) $(EXE) *.exe ../*.exe



.F90.o:
	$(FC) $(FCFLAGS) -c $< -o $@ $(EXTERN_INC) $(EXTERN_LIB)
