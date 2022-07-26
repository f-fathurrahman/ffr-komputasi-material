!*************************************************************
!****                  InitSol                            ****
!*************************************************************
! Starting solution, Gaussian distribution
SUBROUTINE init_sol(NZCoeffs, N_orbitals, NoOfPW)
  USE globals, only: G2Grid
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NoOfPW, N_orbitals
  COMPLEX(8), INTENT(OUT) :: NZCoeffs(N_orbitals, NoOfPW)
  INTEGER :: Iorb, IIndex
  REAL(8) :: nrm, Alpha=1.5D0, G2, X, Y
  INTEGER, ALLOCATABLE :: rnd_state(:)
  INTEGER :: state_size

  CALL random_seed(size=state_size)
  ALLOCATE(rnd_state(state_size))

  rnd_state = 12345
  CALL random_seed(put=rnd_state)

  NZCoeffs = CMPLX(0.D0, kind=8)
  DO IIndex = 1, NoOfPW
    G2 = G2Grid(IIndex)
    DO Iorb = 1, N_orbitals
      nrm = EXP(-Alpha*G2)
      CALL Random_Number(X)
      CALL Random_Number(Y)
      NZCoeffs(Iorb, IIndex) = CMPLX(nrm*X, nrm*Y, kind=8)
    ENDDO
  ENDDO
  CALL Gram_Schmidt(NZCoeffs, N_orbitals, NoOfPW)
  
  DEALLOCATE(rnd_state)

END SUBROUTINE
