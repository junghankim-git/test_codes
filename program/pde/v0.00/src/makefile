
FC      = g++
MPIFC   = mpic++
FCFLAGS = -O2 -fbounds-check -fstack-check

SRCS    = parallel.cpp
OBJS    = $(SRCS:.cpp=.o)

SRC     = main_program.cpp
OBJ     = $(SRC:.cpp=.o)
EXE     = run.exe

INCLUDE = -I.

.SUFFIXES:
.SUFFIXES:.cpp .o


build: $(SRC) $(OBJS)
	$(MPIFC) $(FCFLAGS) -o $(EXE) $^ $(INCLUDE)

clean:
	rm -rf $(OBJS) $(EXE)

test:
	$(MPIFC) $(FCFLAGS) test.cpp

.cpp.o:
	$(MPIFC) $(FCFLAGS) -c $< $(INCLUDE)
