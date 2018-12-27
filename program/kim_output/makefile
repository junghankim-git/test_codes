FC=gfortran
#FCFLAGS=-O2 -ffree-line-length-none -fcheck=all
FCFLAGS=-O2 -ffree-line-length-none

#NETCDF_DIR=/opt/kiaps/apps/netcdf/4.1.3/gcc/4.5.3
NETCDF_DIR=/app/compilers/gcc/4.7.1/applib1/NETCDF4/4.1.3
NETCDF_INC=$(NETCDF_DIR)/include
NETCDF_LIB=$(NETCDF_DIR)/lib
EXTERNAL=-L$(NETCDF_LIB) -lnetcdff -lnetcdf -I$(NETCDF_INC)

VPATH = ./src /home/jhkim/work/share
SRCS = kinds.F90 calendar.F90 time.F90 logger.F90 statistics.F90 list.F90 opt_parser.F90 \
       kim_std_vars.F90 kim_std_output.F90
OBJS = $(SRCS:.F90=.o)
#MODS = $(SRCS:.F90=.o | tr A-Z a-z)

EXE = KIMStdOut.exe
SRC = main_program.F90

.SUFFIXES:
.SUFFIXES:.F90 .o

exe: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) $^ -o $(EXE) $(EXTERNAL) -I/home/jhkim/work/share

gen: main_gen_db.F90 $(OBJS)
	$(FC) $(FCFLAGS) $^ -o gen_db.exe $(EXTERNAL)

mod: test_kim.F90 $(OBJS)
	$(FC) $(FCFLAGS) $^ -o test_kim.exe $(EXTERNAL)

clean:
	rm -rf *.mod *.o *.exe

.F90.o:
	$(FC) $(FCFLAGS) -c $< $(EXTERNAL) -I/home/jhkim/work/share
