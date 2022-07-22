SUBROUTINE get_last(NZCoeffs,NZCoeffsDot, R_IonDot, Time)
  use globals, only: NoOfPW, Ions, N_orbitals, N_ion
  IMPLICIT NONE
  complex(8), INTENT(OUT) :: NZCoeffs(N_orbitals, NoOfPW)
  complex(8), INTENT(OUT) :: NZCoeffsDot(N_orbitals, NoOfPW)
  complex(8), INTENT(OUT) :: R_IonDot(N_ion, 3)
  real(8), INTENT(OUT) :: Time
         
  LOGICAL       :: Exists
  INQUIRE(FILE='last.dat', EXIST=Exists)
  IF (Exists) THEN
    OPEN ( 10, FILE='last.dat')
    READ (10,*) Time
    READ (10,*) NZCoeffs
    READ (10,*) NZCoeffsDot
    READ (10,*) Ions
    READ (10,*) R_IonDot
    CLOSE(10)
  END IF 
END SUBROUTINE