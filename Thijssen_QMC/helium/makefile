# Makefile to compile Quantum Monte Carlo simulation for
# Computational Physics

# COMPHY directory path:
COMPHY	=	$(HOME)/comphy
LIBS	=	-L$(COMPHY)/lib -lran
FFLAGS	=	-O -c
LDFLAGS	=	-O
FC	=	gfortran

helium: helium.o
	$(FC) -o helium $(LDFLAGS) helium.o $(LIBS)

helium.o: helium.f globHe
	$(FC) $(FFLAGS) helium.f
clean:
	rm -f *.o helium
