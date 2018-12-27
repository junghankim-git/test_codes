FC=gfortran
FCFLAGS=-O2 -fcheck=all

all:
	$(FC) $(FCFLAGS) -c buffer.F90
	$(FC) $(FCFLAGS) TestMain.F90 buffer.o -I./


clean:
	rm -rf a.out
