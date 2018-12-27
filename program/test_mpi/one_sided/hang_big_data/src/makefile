FC = mpif90
FCFLAGS = -O2 -ffree-line-length-none -mcmodel=large -fcheck=all -fbacktrace
#FC = mpiifort
#FCFLAGS = -O2 -free

#VPATH = /home/jhkim/work/share

SRCS = kinds.F90 parallel.F90 infra.F90
OBJS = $(SRCS:.F90=.o)
SRC = main_program.F90
EXE = run.exe

.SUFFIXES:
.SUFFIXES:.F90 .o

build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) -o $(EXE) $^
	cp $(EXE) ../

clean:
	rm -rf *.mod *.o $(EXE) ../$(EXE)

.F90.o:
	$(FC) $(FCFLAGS) -c $<
