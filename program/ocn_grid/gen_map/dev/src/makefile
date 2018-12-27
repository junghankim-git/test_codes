
include ../makefile.in

VPATH = /home/jhkim/work/share

SRCS  = kinds.F90 ocn_domain.F90
OBJS  = $(SRCS:.F90=.o)
MODS  = $(shell echo $(SRCS:.F90=.mod) | tr A-Z a-z)
SRC   = main_program.F90
OBJ   = $(SRC:.F90=.o)
EXE   = run.exe

.SUFFIXES:
.SUFFIXES: .F90 .o

build: $(OBJS) $(OBJ)
	$(FC) $(FCFLAGS) -o $(EXE)  $^ $(EXTERN_LIB) $(EXTERN_INC)
	cp $(EXE) ../

clean:
	rm -f $(EXE) $(CEXE) *.o *.mod

.F90.o:
	$(FC) $(FCFLAGS) -c $< -o $@ $(EXTERN_INC) $(EXTERN_LIB)
