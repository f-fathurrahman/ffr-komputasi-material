SUBROUTINE check_const_energy(NZCoeffsDot, R_ionDot)
  use globals
  IMPLICIT NONE
  complex(8), INTENT(IN) :: NZCoeffsDot(N_orbitals, NoOfPW), &
                                R_ionDot(N_ion, 3)
  complex(8) :: E
  INTEGER        :: N

  E = CMPLX(0.D0)
  DO N=1, N_ion
    E = E + SUM(R_ionDot(N,:)*CONJG(R_ionDot(N,:)))/(2*Ions(N)%Mass) 
  END DO
  E = E + SUM(NZCoeffsDot*CONJG(NZCoeffsDot))/(2*mu)
  print '(A23 D15.3 )', 'Constant Energy check:', DBLE(E)
  print *
END SUBROUTINE