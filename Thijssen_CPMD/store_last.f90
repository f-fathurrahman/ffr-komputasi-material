SUBROUTINE store_last(NZCoeffs,NZCoeffsDot, R_IonDot, Time)
  use globals
  IMPLICIT NONE
  complex(8), INTENT(IN) :: NZCoeffs(N_orbitals, NoOfPW)
  complex(8), INTENT(IN) :: NZCoeffsDot(N_orbitals, NoOfPW)
  complex(8), INTENT(IN) :: R_IonDot(N_ion, 3)
  real(8), INTENT(IN) :: Time
  
  OPEN ( 10, FILE='last.dat')
  WRITE (10,*) Time
  WRITE (10,*) NZCoeffs
  WRITE (10,*) NZCoeffsDot
  WRITE (10,*) Ions
  WRITE (10,*) R_IonDot
  CLOSE( 10 )
END SUBROUTINE


  

