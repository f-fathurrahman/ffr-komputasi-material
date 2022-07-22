!*************************************************************
!****                  InitCoeffs                         ****
!*************************************************************
SUBROUTINE init_coeffs(NZCoeffs, NZCoeffsDot, R_ionDot, Time)
  use globals
  IMPLICIT NONE
  LOGICAL                      :: Opt_Exists, Mov_Exists
  complex(8), INTENT(OUT)  :: NZCoeffs(N_orbitals, NoOfPW)
  complex(8), INTENT(OUT)  :: NZCoeffsDot(N_orbitals, NoOfPW)
  complex(8), INTENT(OUT)  :: R_ionDot(N_ion, 3)
  real(8), INTENT(OUT):: Time
  
  INQUIRE(FILE='last.dat', EXIST=Mov_Exists)
  INQUIRE(FILE='opt.dat',EXIST=Opt_Exists)
  IF( Mov_Exists ) THEN
    CALL get_last(NZCoeffs, NZCoeffsDot, R_ionDot,  Time)
  ELSEIF( Opt_Exists ) THEN
    CALL get_optimal(NZCoeffs, NZCoeffsDot)
    Time = 0.D0
  ELSE
    CALL Init_Sol(NZCoeffs, N_Orbitals, NoOfPW)
  ENDIF
END SUBROUTINE