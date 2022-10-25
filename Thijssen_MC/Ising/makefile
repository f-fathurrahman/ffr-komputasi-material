COMPHY	=	$(HOME)/comphy
FC	=	gfortran
XPSLIBS	=	-L$(XPSDIR)/lib -lxps -L/usr/X11R6/lib -lX11 
XLIBDIR	=	/usr/X11R6/lib
FFLAGS	= 	-O -DPlot
# The Plot flag enables plotting using the XPS package.
# This package can also be found in the computational
# physics program directory
# To switch this option off, use -UPlot and remove $(XPSLIBS) below
LIBS = -L$(COMPHY)/lib/ -lran $(XPSLIBS)

IsMC: IsMC.F globIsMC
	$(FC) $(FFLAGS)  -o IsMC IsMC.F $(LIBS)


clean:
	rm -f *.o IsMC
