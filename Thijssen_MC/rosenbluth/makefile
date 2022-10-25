COMPHY	=	$(HOME)/comphy
FC	=	gfortran
FFLAGS	= 	-O -c
LDFLAGS = 	-O
LIBS    =	-L$(COMPHY)/lib -lran

%.o: %.f90
	$(FC) $(FFLAGS) $<

rosenbluth: rosenbluth.o
	$(FC) $(LDFLAGS)  -o rosenbluth rosenbluth.o $(LIBS)

recur: recur.o
	$(FC) $(LDFLAGS)  -o recur recur.o $(LIBS)

perm: perm.o
	$(FC) $(LDFLAGS)  -o perm perm.o $(LIBS)

clean:
	rm -f *.o recur rosenbluth perm
