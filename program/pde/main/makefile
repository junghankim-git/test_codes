include makefile.in

SRC_DIR = src

SRC     = main_program.cpp
OBJ     = $(SRC:.cpp=.o)
EXE     = run.exe

INCLUDE = -I. -I$(SRC_DIR)

.SUFFIXES:
.SUFFIXES:.cpp .o

build: $(SRC)
	make -C src
	$(MPIFC) $(FCFLAGS) -o $(EXE) $^ -L./ -l$(LIBNAME) $(INCLUDE)

clean:
	make -C $(SRC_DIR) $@
	rm -rf $(EXE) lib$(LIBNAME).a

.cpp.o:
	$(MPIFC) $(FCFLAGS) -c $< $(INCLUDE)
