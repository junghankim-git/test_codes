FC      = mpif90
FCFLAGS = -O2 -ffree-line-length-none -fcheck=all
#FC      = ifort
#FCFLAGS = -O2 -free -check=all

VPATH = /home/jhkim/work/share

SRCS = kinds.F90
OBJS = $(SRCS:.F90=.o)
MODS = $(SRCS:.F90=.mod)

SRC  = main.F90
EXE  = run.exe

.SUFFIXES:
.SUFFIXES:.F90 .o


build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) $^ -o $(EXE)

clean:
	rm -f $(OBJS) $(MODS) $(EXE) *.nc

.F90.o:
	$(FC) $(FCFLAGS) -c $<
