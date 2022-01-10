# Makefile for Car-Parrinello Program

COMPILER = gfortran
LDOPTS = -O
COMPOPTS = -c -O  
LIBS = -L/usr/lib/ -lfftw -lblas -llapack 
OBJECTS	= fft.o globals.o utilities.o excorr.o pseudo.o grids.o energy.o force.o main.o

%.o : %.f90
	$(COMPILER) $(COMPOPTS) $<

%.o : %.f
	$(COMPILER) $(COMPOPTS) -I/usr/local/include $<

carpar: $(OBJECTS)
	$(COMPILER) $(OBJECTS) $(LDOPTS) -o carpar.x -L/usr/lib/ -lblas -llapack $(LIBS) 

clean: 
	rm -f *.o *.d *.mod

cleanall:
	rm -f *.o *.d *.mod carpar.x last.dat work*.*
