rm -v *.mod *.o # clean up
gfortran -c -Wall -O3 m_variables_Limc.f90
gfortran -c -Wall -O3 -Wno-unused-dummy-argument m_random_Li_atom.f90
gfortran -c -Wall -O3 m_orbital_Li_HF.f90
gfortran -c -Wall -O3 m_jastrow.f90
gfortran -c -Wall -O3 m_determinant.f90
gfortran -c -Wall -O3 m_observables.f90
gfortran -c -Wall -O3 m_output.f90
gfortran -Wall -O3 main_Li_HF.f90 *.o

