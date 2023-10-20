rm *.o *.a *.mod -v
gfortran -c -O3 -Wall m_highlevel.f90
gfortran -c -O3 -Wall m_midlevel.f90
#gfortran -c -O3 -Wall m_random_H2mc.f90
#gfortran -c -O3 -Wall m_orbital.f90
#gfortran -c -O3 -Wall m_jastrow.f90
#ar rcs libmain.a *.o
#gfortran -O3 H2mc_final.f90 libmain.o
