FC = gfortran
FCFLAGS = -ffree-line-length-none
LAPACKLIB = -L/app/compilers/gcc/4.7.1/applib1/LAPACK/3.4.1 -llapack -lrefblas
SHAREDDIR = /home/jhkim/Study/Library/Shared/etc
LIBDIR = $(SHAREDDIR)/lib
INCDIR = $(SHAREDDIR)/include

MYLIB = -L$(LIBDIR) -lshared -I$(INCDIR)

SRC = TestMain.F90

exe: $(SRC)
	make -C /home/jhkim/Study/Library/Shared
	$(FC) $(FCFLAGS) $^ $(MYLIB) $(LAPACKLIB)

clean:
	rm -rf *.mod *.o a.out

