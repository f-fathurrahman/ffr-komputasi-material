!*************************************************************
!****                  InitCoeffs                         ****
!*************************************************************
SUBROUTINE init_coeffs(NZCoeffs, NZCoeffsDot, R_ionDot, Time)
  USE globals
  IMPLICIT NONE
  LOGICAL                      :: Opt_Exists, Mov_Exists
  COMPLEX(8), INTENT(OUT)  :: NZCoeffs(N_orbitals, NoOfPW)
  COMPLEX(8), INTENT(OUT)  :: NZCoeffsDot(N_orbitals, NoOfPW)
  COMPLEX(8), INTENT(OUT)  :: R_ionDot(N_ion, 3)
  REAL(8), INTENT(OUT):: Time
  
  INQUIRE(FILE='last.dat', EXIST=Mov_Exists)
  INQUIRE(FILE='opt.dat',EXIST=Opt_Exists)
  !
  IF( Mov_Exists ) THEN
    WRITE(*,*) 'Read ionic velocities from file'
    CALL get_last(NZCoeffs, NZCoeffsDot, R_ionDot,  Time)
  !
  ELSEIF( Opt_Exists ) THEN
    WRITE(*,*) 'Read orbitals from file: '
    CALL get_optimal(NZCoeffs, NZCoeffsDot)
    Time = 0.D0
  !
  ELSE
    WRITE(*,*) 'Initialize orbitals from scratch'
    CALL init_sol(NZCoeffs, N_Orbitals, NoOfPW)
  ENDIF

END SUBROUTINE

