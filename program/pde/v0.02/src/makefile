include ../makefile.in

VPATH   = $(SRC_DIR)

SRCS    = carray.cpp parallel.cpp par_func.cpp
OBJS    = $(SRCS:.cpp=.o)

INCLUDE = -I.

.SUFFIXES:
.SUFFIXES:.cpp .o


build: $(OBJS)
	ar -rcs lib$(LIBNAME).a $^
	cp lib$(LIBNAME).a ../

clean:
	rm -rf $(OBJS) lib$(LIBNAME).a

.cpp.o:
	$(MPIFC) $(FCFLAGS) -c $< $(INCLUDE)
