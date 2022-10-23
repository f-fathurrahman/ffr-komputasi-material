# Makefile to compile Quantum Monte Carlo simulation for
# Computational Physics

# COMPHY directory path:
COMPHY	=	$(HOME)/comphy
LIBS	=	-L$(COMPHY)/lib -lran
FFLAGS	=	-O -c
LDFLAGS	=	-O
FC	=	gfortran

hydrogen: hydrogen.o
	$(FC) -o hydrogen $(LDFLAGS) hydrogen.o $(LIBS)

hydrogen.o: hydrogen.f globHydr
	$(FC) $(FFLAGS) hydrogen.f
clean:
	rm -f *.o hydrogen
