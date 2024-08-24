EXE=pw_dft.x
FC=gfortran
FCFLAGS=-O2
LINKFLAGS=-lblas -llapack

FFT_DIR=lib

OBJS=pw_dft.o gvector.o pw_small.o pw.o ps.o fft_fftpack_data.o \
band.o density.o dos.o fermi.o gs.o hamiltonian_small.o hpsi.o \
nl_ps.o phfac.o potential.o scf_dft.o solve_small.o steepest_descent.o \
strucf.o sw.o subs.o sym.o iflip.o lda_xc.o linalg.o order.o

FFT_OBJS=cfftb1.o cfftf.o fft_fftpack.o passb3.o passb.o passf4.o \
cfftb.o cffti1.o fft_fftpack_fast.o passb4.o passf2.o passf5.o \
cfftf1.o cffti.o passb2.o passb5.o passf3.o passf.o param.o

$(EXE): $(FFT_OBJS) $(OBJS)
	$(FC) -o $(EXE) $(FFT_OBJS) $(OBJS) $(LINKFLAGS)

pw_dft.o: pw_dft.f90 pw.o
	$(FC) $(FCFLAGS) -c pw_dft.f90

gvector.o: gvector.f90 iflip.o order.o linalg.o
	$(FC) $(FCFLAGS) -c gvector.f90

iflip.o: iflip.f90 
	$(FC) $(FCFLAGS) -c iflip.f90
	
lda_xc.o: lda_xc.f90
	$(FC) $(FCFLAGS) -c lda_xc.f90

linalg.o: linalg.f90
	$(FC) $(FCFLAGS) -c linalg.f90
			
order.o: order.f90
	$(FC) $(FCFLAGS) -c order.f90	       
	
pw_small.o: pw_small.f90 gvector.o iflip.o
	$(FC) $(FCFLAGS) -c pw_small.f90
	
pw.o: pw.f90 gvector.o pw_small.o ps.o fft_fftpack_data.o
	$(FC) $(FCFLAGS) -c pw.f90	

ps.o: ps.f90 gvector.o pw_small.o
	$(FC) $(FCFLAGS) -c ps.f90

fft_fftpack_data.o: fft_fftpack_data.f90
	$(FC) $(FCFLAGS) -c fft_fftpack_data.f90
	
band.o: band.f90 fft_fftpack_data.o gvector.o pw.o pw_small.o ps.o
	$(FC) $(FCFLAGS) -c band.f90	

density.o: density.f90 fft_fftpack_data.o gvector.o pw.o ps.o 
	$(FC) $(FCFLAGS) -c density.f90
	
dos.o: dos.f90 gvector.o pw.o ps.o 
	$(FC) $(FCFLAGS) -c dos.f90	
	
fermi.o: fermi.f90 gvector.o pw.o ps.o linalg.o
	$(FC) $(FCFLAGS) -c fermi.f90	
	
gs.o: gs.f90 pw.o 
	$(FC) $(FCFLAGS) -c gs.f90	
	
hamiltonian_small.o: hamiltonian_small.f90 fft_fftpack_data.o gvector.o pw.o pw_small.o ps.o linalg.o
	$(FC) $(FCFLAGS) -c hamiltonian_small.f90	
	
hpsi.o: hpsi.f90 fft_fftpack_data.o gvector.o pw.o ps.o 
	$(FC) $(FCFLAGS) -c hpsi.f90		
	
nl_ps.o: nl_ps.f90 gvector.o pw.o ps.o iflip.o
	$(FC) $(FCFLAGS) -c nl_ps.f90			
	
phfac.o: phfac.f90 gvector.o pw.o ps.o iflip.o
	$(FC) $(FCFLAGS) -c phfac.f90	
	
potential.o: potential.f90 fft_fftpack_data.o gvector.o pw.o ps.o lda_xc.o
	$(FC) $(FCFLAGS) -c potential.f90	
	
scf_dft.o: scf_dft.f90 fft_fftpack_data.o gvector.o pw.o pw_small.o ps.o 
	$(FC) $(FCFLAGS) -c scf_dft.f90		
	
solve_small.o: solve_small.f90 fft_fftpack_data.o gvector.o pw.o pw_small.o ps.o 
	$(FC) $(FCFLAGS) -c solve_small.f90	
	
steepest_descent.o: steepest_descent.f90 fft_fftpack_data.o gvector.o pw.o pw_small.o ps.o 
	$(FC) $(FCFLAGS) -c steepest_descent.f90	
		
strucf.o: strucf.f90 gvector.o ps.o iflip.o
	$(FC) $(FCFLAGS) -c strucf.f90		

sw.o: sw.f90 fft_fftpack_data.o gvector.o pw.o pw_small.o ps.o 
	$(FC) $(FCFLAGS) -c sw.f90	
				
subs.o: subs.f90
	$(FC) $(FCFLAGS) -c subs.f90

sym.o: sym.f90 fft_fftpack_data.o gvector.o
	$(FC) $(FCFLAGS) -c sym.f90	
	
	
param.o: $(FFT_DIR)/param.f90
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/param.f90
	
cfftb1.o: $(FFT_DIR)/cfftb1.f
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/cfftb1.f

cfftf.o: $(FFT_DIR)/cfftf.f
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/cfftf.f

fft_fftpack.o: $(FFT_DIR)/fft_fftpack.f90 param.o fft_fftpack_data.o
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/fft_fftpack.f90
	
passb3.o: $(FFT_DIR)/passb3.f
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/passb3.f

passb.o: $(FFT_DIR)/passb.f
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/passb.f

passf4.o: $(FFT_DIR)/passf4.f
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/passf4.f

cfftb.o: $(FFT_DIR)/cfftb.f
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/cfftb.f

cffti1.o: $(FFT_DIR)/cffti1.f
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/cffti1.f

fft_fftpack_fast.o: $(FFT_DIR)/fft_fftpack_fast.f
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/fft_fftpack_fast.f

passb4.o: $(FFT_DIR)/passb4.f
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/passb4.f 

passf2.o: $(FFT_DIR)/passf2.f
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/passf2.f 

passf5.o: $(FFT_DIR)/passf5.f
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/passf5.f

cfftf1.o: $(FFT_DIR)/cfftf1.f
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/cfftf1.f 

cffti.o: $(FFT_DIR)/cffti.f
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/cffti.f 

passb2.o: $(FFT_DIR)/passb2.f
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/passb2.f 

passb5.o: $(FFT_DIR)/passb5.f
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/passb5.f 

passf3.o: $(FFT_DIR)/passf3.f
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/passf3.f

passf.o: $(FFT_DIR)/passf.f
	$(FC) $(FCFLAGS) -c $(FFT_DIR)/passf.f
	
clean:
	-rm -rf *.o *.mod $(FFT_DIR)/*.o $(FFT_DIR)/*.mod
cleaner:
	-rm -rf *.o *.mod $(FFT_DIR)/*.o $(FFT_DIR)/*.mod $(EXE)
	
