FC = mpif90
FCFLAGS = -O2 -ffree-line-length-none -fsecond-underscore -fcheck=all

VPATH = /home/jhkim/work/share
VER=4.7.1
NETCDF_DIR=/app/compilers/gcc/$(VER)/mpi/mvapich2/1.8.1/applib2/NETCDF4/4.1.3
EXT_LIBS=-L$(NETCDF_DIR)/lib -lnetcdff -lnetcdf -I$(NETCDF_DIR)/include


SRCS = kinds.F90 parallel.F90 utilities.F90
OBJS = $(SRCS:.F90=.o)

SRC  = main_program.F90
EXEC = run.exe


.SUFFIXES:
.SUFFIXES:.F90 .o

exe: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) -o $(EXEC) $^
	cp $(EXEC) ../

clean:
	rm -rf *.mod *.o *.exe

.F90.o:
	$(FC) $(FCFLAGS) -c $<
