CC      = mpicc
CCFLAGS = -O2

SRCS = class_fdm.cpp
OBJS = $(SRCS:.cpp=.o)

SRC  = main_program.cpp
OBJ  = $(SRC:.cpp=.o)
EXE  = run.exe

.SUFFIXES:
.SUFFIXES:.cpp .o


build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) -o $(EXE) $^

clean:
	rm -rf $(OBJS) $(EXE)

run:
	mpirun -n 4 ./$(EXE)


.cpp.o:
	$(FC) $(FCFLAGS) -c $<
