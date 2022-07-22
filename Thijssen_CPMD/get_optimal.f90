SUBROUTINE get_optimal(NZCoeffs, NZCoeffsDot)
  use globals
  IMPLICIT NONE
  LOGICAL       :: Exists
  complex(8), INTENT(OUT) :: NZCoeffs(N_orbitals, NoOfPW)
  complex(8), INTENT(OUT) :: NZCoeffsDot(N_orbitals, NoOfPW)
  
  INQUIRE (FILE='opt.dat', EXIST=Exists)
  IF (Exists) THEN
    OPEN ( 10, FILE='opt.dat')
    READ (10,*) NZCoeffs
    READ (10,*) NZCoeffsDot
    CLOSE(10)
  END IF
END SUBROUTINE