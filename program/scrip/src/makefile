include ../makefile.in

VPATH = /home/jhkim/work/share

SRCS=kinds.F90 coordinates.F90 scrip_input.F90 scrip_khkim.F90 scrip_hycom.F90 buffer.F90 \
     remap_netcdf.F90 remap_matrix.F90 remapper.F90
OBJS=$(SRCS:.F90=.o)
MODS=$(shell echo $(OBJS:.o=.mod) | tr A-Z a-z)

build: builddir $(OBJS)
#	make -C ./ builddir
#	make -C ./ $(OBJS)
	ar -cr lib$(LIBNAME).a $(OBJS)
	cp $(MODS) $(BUILDDIR)/include
	cp lib$(LIBNAME).a $(BUILDDIR)/lib

clean:
	rm -f $(OBJS) $(MODS) lib$(LIBNAME).a

distclean: clean
	rm -f $(BUILDDIR)/lib/lib$(LIBNAME).a $(addprefix $(BUILDDIR)/include/, $(MODS))

builddir:
	mkdir -p $(BUILDDIR)/lib
	mkdir -p $(BUILDDIR)/include

.SUFFIXES:
.SUFFIXES:.F90 .o

.F90.o:
	$(FC) $(FCFLAGS) -c $< $(EXTERNAL)
