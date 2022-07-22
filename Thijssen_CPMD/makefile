# Makefile for Car-Parrinello Program

COMPILER = gfortran
LDOPTS = -O
COMPOPTS = -c -O -Wall
LIBS = -L/usr/lib/ -lfftw -lblas -llapack 
OBJECTS	= fft.o globals.o \
calc_dens_and_coeffs_r.o calc_orbitals.o  get_E_ovrl.o   \
init_coeffs.o  nonloc_pp.o store_optimal.o \
calc_facs.o check_const_energy.o get_E_self.o init_grids.o periodic_boundary.o  total_energy.o \
calc_F_ES.o core_dens.o get_index_pp.o init_params.o     print_ions.o  Vxc.o \
calc_F_local.o epsilon_xc.o  get_last.o  init_pp_params.o  rattle.o \
calc_F_nonlocal.o get_optimal.o  init_sol.o run_carpar.o \
calc_F_ovrl.o fill_dyn_grids.o  get_pp_filename.o  inner_prod.o      short_pp.o \
calc_ion_force.o fill_static_grids.o  local_pp.o   simulate.o \
calc_orb_force.o G2_short.o  gram_schmidt.o  main.o  store_last.o


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
