# Makefile to compile Monte Carlo simulation for
# Computational Physics

# COMPHY directory path:
COMPHY	=	$(HOME)/comphy
# XPSDIR is an environment variable containing the path to the "lib" subdir
# of the XPS directory
XPSLIBS	= 	-L$(XPSDIR)/lib -L/usr/X11R6/lib -lxps -lX11
LIBS	=	-L$(COMPHY)/lib  -lran $(XPSLIBS)
FFLAGS	=	-O -c -DPlot
# The Plot flag enables plotting using the XPS package.
# This package can also be found in the computational
# physics program directory
# To switch this option off, use -UPlot
# In that case, $(XPSLIBS) can be omitted above!!!!!
LDFLAGS	=	-O 
FC	=	gfortran

mc: mc.o
	$(FC) -o mc $(LDFLAGS) mc.o $(LIBS)

mc.o: mc.F globmc.f
	$(FC) $(FFLAGS) mc.F
clean:
	rm -f *.o mc
