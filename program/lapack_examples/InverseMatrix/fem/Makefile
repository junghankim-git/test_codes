FC = gfortran
SHAREDDIR=/home/jhkim/Study/NumericalMethod/Shared/etc
LIBDIR = $(SHAREDDIR)/lib
INCDIR = $(SHAREDDIR)/include
MYLIB  = -L$(LIBDIR) -lshared -I$(INCDIR)

EXTLIB = -L/app/compilers/gcc/4.7.1/applib1/LAPACK/3.4.1 -llapack -lrefblas


SRC = InverseMatrix.F90
GMRESDIR = /home/jhkim/Study/NumericalMethod/Main/LinearEquation/GMRES
GMRESSRC = GMRES_Solver.F90
GMRESOBJ = GMRES_Solver.o

SUMSRC = ReproSum.F90
SUMOBJ = ReproSum.o

#exe: $(SRC) $(GMRESOBJ)
#exe: $(SRC) $(GMRESOBJ) $(SUMOBJ)
exe: $(SRC)
	$(FC) $^ $(MYLIB) $(EXTLIB)

$(GMRESOBJ): $(GMRESDIR)/$(GMRESSRC)
	$(FC) -c $^

$(SUMOBJ): $(SUMSRC)
	$(FC) -c $^

clean:
	rm -rf *.o *.mod a.out
