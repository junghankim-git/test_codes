FC       = gfortran
FCFLAGS  = -O2

VPATH = /home/jhkim/work/share

SRCS = kinds.F90 io_ascii.F90 cubic_spline_interp.F90
OBJS = $(SRCS:.F90=.o)
SRC  = main_program.F90
EXE  = run.exe

.SUFFIXES:
.SUFFIXES:.F90 .o


exe: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) -o $(EXE) $^

clean:
	rm -rf *.mod *.o $(EXE)

.F90.o:
	$(FC) $(FCFLAGS) -c $^
