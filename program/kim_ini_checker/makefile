FC=ftn
FCFLAGS=-O2 -free -no-wrap-margin -fp-model precise -check bounds

VPATH = . ../../share/

SRCS = kinds.F90 chk_kim_ini.F90
OBJS = $(SRCS:.F90=.o)
MODS = $(SRCS:.F90=.mod)

SRC = main_program.F90
EXE = run.exe

EXETERN = -lnetcdff -lnetcdf

.SUFFIXES:
.SUFFIXES: .F90 .o


build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) -o $(EXE) $^ $(EXTERN)


clean:
	rm -rf $(OBJS) $(MODS) $(EXE)


.F90.o:
	$(FC) $(FCFLAGS) -c $^ $(EXTERN)
