
FC=gfortran
MPIFC=mpif90
FCFALGS=-O2 -ffree-line-length-none

SRCS=Infra.F90 KIO.F90 Computation.F90 Model.F90
OBJS=$(SRCS:.F90=.o)
MODS=$(shell echo $(SRCS:.F90=.mod) | tr A-Z a-z)

.SUFFIXES:
.SUFFIXES: .F90 .o

all: $(OBJS)
	$(MPIFC) $(FCFLAGS) test.F90 $^

clean:
	rm -rf $(OBJS) $(MODS) a.out

run:
	mpirun -n 5 ./a.out

.F90.o:
	$(MPIFC) -c $(FCFALGS) $<
