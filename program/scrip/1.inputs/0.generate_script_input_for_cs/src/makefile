include ../../../makefile.in


#SRCS = kinds.F90 khkim2scrip.F90 hycom_grid.F90
#OBJS = $(SRCS:.F90=.o)
#MODS = $(SRCS:.F90=.mod)

#EXE = gen_input.exe
#SRC = generate_scrip_input.F90
EXE = run.exe
SRC = main_program.F90

exe: $(SRC)
	make -C ../../../
	$(FC) $(FCFLAGS) $^ -o $(EXE) $(ACV) $(EXTERNAL)
	cp $(EXE) ../

convert: convert_to_radian.F90
	make -C ../../../
	$(FC) $(FCFLAGS) $^ -o convert.exe $(ACV) $(EXTERNAL)
	cp convert.exe ../

check: check_scrip_input.F90
	make -C ../../../
	$(FC) $(FCFLAGS) $^ -o check.exe $(ACV) $(EXTERNAL)
	cp check.exe ../

clean:
	rm -rf *.mod *.o *.exe ../*.exe

.SUFFIXES:
.SUFFIXES:.F90 .o

.F90.o:
	$(FC) $(FCFLAGS) -c $< $(ACV) $(EXTERNAL)
