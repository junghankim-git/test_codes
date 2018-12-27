FC = gfortran
CC = gcc


SRC = TestMain.F90
LIB = -L/app/compilers/gcc/4.7.1/applib1/FFTW3/3.3.2/lib -lfftw3
INC = -I/app/compilers/gcc/4.7.1/applib1/FFTW3/3.3.2/include


exe: $(SRC)
	$(FC) $^ $(LIB) $(INC)
