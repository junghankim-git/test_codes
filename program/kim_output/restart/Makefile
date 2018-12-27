FC      = gfortran
FCFLAGS = -O2

SRCS    = CheckRestart.F90
OBJS    = $(SRCS:.F90=.o)
EXE     = run.exe

NETCDF_DIR=/app/compilers/gcc/4.7.1/applib1/NETCDF4/4.1.3
NETCDF_INC=$(NETCDF_DIR)/include
NETCDF_LIB=$(NETCDF_DIR)/lib
EXTERNAL=-L$(NETCDF_LIB) -lnetcdff -lnetcdf -I$(NETCDF_INC)


build: $(OBJS)
	$(FC) $(FCFLAGS) TestMain.F90 $^ -o $(EXE) $(EXTERNAL)


clean:
	rm -rf *.mod *.o *.exe


.SUFFIXES:
.SUFFIXES:.F90 .o

.F90.o:
	$(FC) $(FCFLAGS) -c $< $(EXTERNAL)
