!*************************************************************
!****                  InitSol                            ****
!*************************************************************
! Starting solution, Gaussian distribution
SUBROUTINE init_sol(NZCoeffs, N_orbitals, NoOfPW)
  use globals, only: G2Grid
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NoOfPW, N_orbitals
  complex(8), INTENT(OUT) :: NZCoeffs(N_orbitals, NoOfPW)
  INTEGER :: Iorb, IIndex
  real(8) :: Norm, Alpha=1.5D0, G2, X, Y
  
  NZCoeffs = CMPLX(0.D0)
  DO IIndex = 1, NoOfPW
    G2 = G2Grid(IIndex)
    DO Iorb = 1, N_orbitals
      Norm = EXP(-Alpha*G2)
      CALL Random_Number(X)
      CALL Random_Number(Y)
      NZCoeffs(Iorb, IIndex) = CMPLX(Norm*X, Norm*Y)
    END DO
  END DO
  CALL Gram_Schmidt(NZCoeffs, N_orbitals, NoOfPW)
  
END SUBROUTINE