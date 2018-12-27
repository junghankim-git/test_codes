#FC = gfortran
FCFLAGS = $(FCFLAGS) -fbounds-check

SHAREDDIR = /home/jhkim/Study/Library/Shared/etc
include $(SHAREDDIR)/../Makefile.default
LIBDIR = $(SHAREDDIR)/lib
INCDIR = $(SHAREDDIR)/include
MYLIB = -L$(LIBDIR) -lshared -I$(INCDIR)

SRC = ReadNetCDF.F90

exe: $(SRC)
	$(MPIFC) $(FCFLAGS) $^ $(MYLIB) $(EXTINC) $(EXTLIB)



clean:
	rm -rf *.mod *.o a.out

