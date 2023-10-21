MAKE = make
AR = ar

F90 = gfortran
F90_OPTS = -Wall -O3
F90_OPTS_NOWARN = -cpp -O3

CC = gcc -cpp -I./
CC_OPTS = 

LIBS =

SRC = m_constants.f90 m_random.f90 m_position.f90

OBJ = $(SRC:.f90=.o)

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


# Libraries
lib: $(OBJ)
	ar rcs libmain.a *.o

# Main executable
main: lib
	$(F90) $(F90_OPTS) main_Hmc.f90 -o main.x libmain.a $(LIBS)

clean:
	rm -rfv *.o *.mod *.x libmain.a


