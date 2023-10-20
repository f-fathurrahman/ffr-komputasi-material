rm -v *.mod *.o # clean up
gfortran -c -Wall -O3 m_variables_cubicqd_cluster.f90
gfortran -c -Wall -O3 -Wno-unused-dummy-argument m_random_cubicqd_cluster.f90
gfortran -c -Wall -O3 m_orbital_cubicqd_cluster.f90
gfortran -c -Wall -O3 -Wno-maybe-uninitialized m_determinant_cubicqd_cluster.f90
gfortran -c -Wall -O3 m_jastrow_cubicqd_cluster.f90
gfortran -c -Wall -O3 m_observables_cubicqd_cluster.f90
gfortran -c -Wall -O3 m_output_cubicqd_cluster.f90
ar rcs libmain.a *.o
gfortran -Wall -O3 cubicqd_cluster.f90 libmain.a