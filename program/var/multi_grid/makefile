FC      = mpif90
FCFLAGS = -O2 -ffree-line-length-none -fcheck=all -fbacktrace

SRCS = mg_comm.F90
OBJS = $(SRCS:.F90=.o)
MODS = $(SRCS:.F90=.mod)

SRC = main_program.F90
OBJ = $(SRC:.F90=.o)

EXE = run.exe

.SUFFIXES:
.SUFFIXES:.F90 .o

build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) -o $(EXE) $^

clean:
	rm -rf $(OBJS) $(MODS) $(EXE)

run:
	mpirun -n 4 ./$(EXE)


.F90.o:
	$(FC) $(FCFLAGS) -c $<
