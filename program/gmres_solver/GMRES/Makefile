FC = gfortran
LIB = -L/app/compilers/gcc/4.7.1/applib1/LAPACK/3.4.1 -llapack -lrefblas


#SRC = LeastSquares.F90
#SRC = LinearEqu.F90
SRC = GMRES_Solver.F90
OBJ = GMRES_Solver.o


$(OBJ): $(SRC)
	$(FC) -c $(SRC) $(LIB)
