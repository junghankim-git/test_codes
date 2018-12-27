FC = gfortran
FCFLAGS = -O2

VPATH = /home/jhkim/work/share

SRCS = kinds.F90
OBJS = $(SRCS:.F90=.o)
MODS = $(shell echo $(SRCS:.F90=.mod) | tr A-Z a-z)


SRC = test_lagrange.F90
EXE = run.exe

.SUFFIXES:
.SUFFIXES:.F90 .o


exe: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) -o $(EXE) $^

clean:
	rm -rf *.mod *.o a.out

.F90.o:
	$(FC) $(FCFALGS) -c $<
