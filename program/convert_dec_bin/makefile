FC = gfortran
FCFLAGS = -O2

VPATH = /home/jhkim/work/share

SRCS = binary_decimal.F90
OBJS = $(SRCS:.F90=.o)

SRC  = main_program.F90
EXEC = run.exe

.SUFFIXES:
.SUFFIXES:.F90 .o

build: $(SRC) $(OBJS)
	$(FC) $(FCFLAGS) -o $(EXEC) $^ -I./

clean:
	rm -rf *.o *.mod $(EXEC)

.F90.o:
	$(FC) $(FCFLAGS) -c $<
