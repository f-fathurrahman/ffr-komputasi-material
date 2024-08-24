gfortran -O3 -c m_qd2d.f90
gfortran -O3 -c xc.f90
gfortran -O3 2dqdot.f90 xc.o m_qd2d.o