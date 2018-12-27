
include ../makefile.in

VPATH = /home/jhkim/work/share

SRCS  = kinds.F90 buffer.F90 remap_matrix.F90 remapper.F90
OBJS  = $(SRCS:.F90=.o)
MODS  = $(shell echo $(SRCS:.F90=.mod) | tr A-Z a-z)
SRC   = remap_main.F90
OBJ   = $(SRC:.F90=.o)
EXE   = remap.exe
CSRC  = remap_check_main.F90
COBJ  = $(CSRC:.F90=.o)
CEXE  = checker.exe

.SUFFIXES:
.SUFFIXES: .F90 .o

build: $(OBJS) $(OBJ)
	$(FC) $(FCFLAGS) -o $(EXE)  $^ $(EXTERN_LIB) $(EXTERN_INC)
	cp $(EXE) ../

check: $(OBJS) $(COBJ)
	$(FC) $(FCFLAGS) -o $(CEXE) $^ $(EXTERN_LIB) $(EXTERN_INC)
	cp $(CEXE) ../

convert: ll_convert.F90 kinds.o
	$(FC) $(FCFLAGS) -o convert.exe $^ $(EXTERN_LIB) $(EXTERN_INC)
	cp convert.exe ../

clean:
	rm -f $(EXE) $(CEXE) *.o *.mod

.F90.o:
	$(FC) $(FCFLAGS) -c $< -o $@ $(EXTERN_INC) $(EXTERN_LIB)
