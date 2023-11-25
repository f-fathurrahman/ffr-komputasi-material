# COMPHY: path of COMPHY directory
COMPHY	=	$(HOME)/comphy
FC	=	gfortran
LIBS	=	
FCOPTS	=	-O -c
LDOPTS	=	-O

%o:%.f
	$(FC) $(FCOPTS) $<

elec: elec_md.f elec.glob
	$(FC) $(LDOPTS) elec_md.f -o elec $(LIBS)

nucl: nucl_md.f nucl.glob
	$(FC) $(LDOPTS) nucl_md.f -o nucl $(LIBS)

clean:
	rm -f *.o elec nucl
