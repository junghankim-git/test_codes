FC=gfortran
FCFLAGS=-O2 -ffree-line-length-none -fcheck=all

#NETCDF_DIR=/opt/kiaps/apps/netcdf/4.1.3/gcc/4.5.3
NETCDF_DIR=/app/compilers/gcc/4.7.1/applib1/NETCDF4/4.1.3
NETCDF_INC=$(NETCDF_DIR)/include
NETCDF_LIB=$(NETCDF_DIR)/lib
EXTERNAL=-L$(NETCDF_LIB) -lnetcdff -lnetcdf -I$(NETCDF_INC)


SRCS = ResHis2Post.F90
OBJS = $(SRCS:.F90=.o)
EXE  = genpost.exe

SRC = GenPost.F90


exe: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) $^ -o $(EXE) $(EXTERNAL)

clean:
	rm -rf *.mod *.o *.exe

.SUFFIXES:
.SUFFIXES:.F90 .o

.F90.o:
	$(FC) $(FCFLAGS) -c $< $(EXTERNAL)
