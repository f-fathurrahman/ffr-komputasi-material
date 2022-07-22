SUBROUTINE store_optimal(NZCoeffs, NZCoeffsDot)
  use globals
  IMPLICIT NONE
  complex(8), INTENT(IN) :: NZCoeffs(N_orbitals, NoOfPW)
  complex(8), INTENT(IN) :: NZCoeffsDot(N_orbitals, NoOfPW)

  OPEN ( 10, FILE='opt.dat')
  WRITE (10,*) NZCoeffs
  WRITE (10,*) NZCoeffsDot
  CLOSE( 10 )
END SUBROUTINE