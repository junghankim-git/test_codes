FC    = gfortran

#SHAREDDIR=/home/jhkim/Study/NumericalMethod/Shared/etc
#LIBDIR = $(SHAREDDIR)/lib
#INCDIR = $(SHAREDDIR)/include
#MYLIB  = -L$(LIBDIR) -lshared -I$(INCDIR)

EXTLIB = -L/app/compilers/gcc/4.7.1/applib1/LAPACK/3.4.1 -llapack -lrefblas

SRCDIR=shyoon
OBJS  = RealSPH.o

SRC   = SingularVectorDecompse.F90

exe: $(SRC) $(OBJS)
	$(FC) $^ $(EXTLIB) #$(MYLIB)

clean:
	rm -rf a.out *.o *.mod

RealSPH.o: $(SRCDIR)/RealSPH.F90
	$(FC) -c $^
