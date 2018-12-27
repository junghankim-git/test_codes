
include ./makefile.default

SRCS = kinds.F90 constant.F90 parallel.F90 par_env.F90 par_mpi.F90 par_func.F90 logger.F90    \
       opt_parser.F90 controls.F90 quadrature.F90 element.F90 lapack.F90 matrix_control.F90   \
       repro_sum.F90 std_atmosphere.F90 cubic_spline_interp.F90 polynomial.F90                \
       polynomial_lagrange.F90 piecewise_method.F90 pseudo_inverse.F90 semi_lagrangian.F90    \
       spline_method.F90 statistics.F90 

# todo: atmos_remap.F90 

OBJS = $(SRCS:.F90=.o)
MODS = $(SRCS:.F90=.mod)

.SUFFIXES:
.SUFFIXES:.F90 .o


build: $(OBJS)
	mkdir -p $(ETCLIB) $(ETCINC)
	ar cr libshared.a $^

install:
	mv $(MODS) $(ETCINC)
	cp libshared.a $(ETCLIB)

clean:
	rm -f $(OBJS) $(MODS) libshared.a

distclean: clean
	$(RM) -rf $(LIB_DIR)/etc


.F90.o:
	$(MPIFC) $(FCFLAGS) -c $^ $(EXTLIB) $(EXTINC)
