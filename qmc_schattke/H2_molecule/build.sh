gfortran -c -O3 m_highlevel.f90
gfortran -c -O3 m_midlevel.f90
gfortran -c -O3 m_random_H2mc.f90
gfortran -c -O3 m_orbital.f90
gfortran -c -O3 m_jastrow.f90
gfortran -O3 H2mc_final.f90 m_*.o
