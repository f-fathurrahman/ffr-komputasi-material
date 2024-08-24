include make.inc

LIBS_EXT = $(LIBS)


SRC = \
m_qd3d.f90 \
xc.f90 \
3dqdot.f90

C_SRC = 

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
main: lib 
	$(F90) $(F90_OPTS) main_qd3d.f90 -o $(EXE_MAIN) libmain.a $(LIBS_EXT)

# does not delete *.x files
clean:
	rm -rf *.o *.mod libmain.a

# also delete *.x files
cleanx:
	rm -rf *.o *.mod libmain.a *.x
