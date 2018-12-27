FC = gfortran
FCFLAGS = -ffree-line-length-none

VPATH = /home/jhkim/work/share

SRCS = kinds.F90 mesh_control.F90 space_filling_curve.F90
OBJS = $(SRCS:.F90=.o)
SRC = main_program.F90
EXE = run.exe

.SUFFIXES:
.SUFFIXES:.F90 .o

build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) $^ -o $(EXE)

clean:
	rm -rf *.mod *.o $(EXE)

.F90.o:
	$(FC) $(FCFLAGS) -c $<
