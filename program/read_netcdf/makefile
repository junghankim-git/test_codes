

SRC = ReadNetCDF_10step.F90

exe: $(SRC)
	$(FC) $(FCFLAGS) $^ $(MYLIB) $(EXTINC) $(EXTLIB)



clean:
	rm -rf *.mod *.o a.out

