include make.inc

LIBS_EXT = -lblas -llapack -lgsl -lgslcblas -lfftw3

SRC = mesh.f90 output.f90 \
coeff.f90 \
fftw3.f90 cf.f90 poisson.f90 \
epot.f90 states.f90 energy.f90 \
ipot.f90  perturbation.f90  qdf.f90  sf.f90 \
cg.f90  excitations.f90  gs.f90  read_wfs.f90 \
expokit.f90  hpsi.f90  propagate.f90  rho.f90  td.f90  vxc.f90 \
exponential.f90  init_wfs.f90  propagator.f90    scf.f90  write_wfs.f90 \
expokit.f90 hpsi.f90 init_wfs.f90 \
test_hartree.f90 test_exp.f90 test_laplacian.f90


C_SRC = gslwrappers.c

OBJ = $(SRC:.f90=.o) $(SRC:.f=.o) $(C_SRC:.c=.o)

#
# Suffix rule for Fortran 90
#
%.mod :
	@if [! -f $@ ]; then \
		rm $(*F).o; \
		fi
	$(MAKE) $<

%.o : %.f90
	$(F90) $(F90_OPTS) -c -o $(*F).o $<

#
# Fortran 77 sources
# supress warning
.SUFFIXES: .o .f
.f.o:
	$(F77) -c $(F77_OPTS) $<

#
# C source
#
.SUFFIXES: .o .c
.c.o:
	$(CC) -c $(CC_OPTS) $<

# Libraries
lib: $(OBJ)
	ar rcs libmain.a *.o


#
# Main executable
#

qdmain: lib
	$(F90) $(F90_OPTS) qdmain.f90 -o qdmain.x libmain.a $(LIBS_EXT)

clean:
	rm -rfv *.o *.mod *.x libmain.a


