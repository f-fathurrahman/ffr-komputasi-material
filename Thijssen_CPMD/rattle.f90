SUBROUTINE Rattle(NZCoeffs, OldNZCoeffs, NZCoeffsdot)
  use globals
  IMPLICIT NONE
  complex(8), INTENT(INOUT) :: NZCoeffs(N_Orbitals, NoOfPW), &
                                   NZCoeffsdot(N_Orbitals, NoOfPW), &
                                   OldNZCoeffs(N_Orbitals, NoOfPW)

  INTEGER :: N, J, K
  real(8) :: Norm
  complex(8), ALLOCATABLE :: I(:,:), A(:,:), B(:,:), &
                                 X(:,:), Correction(:,:)
  IF (.NOT.ALLOCATED(A)) THEN
    ALLOCATE (A(N_Orbitals, N_Orbitals))
    ALLOCATE (B(N_Orbitals, N_Orbitals))
    ALLOCATE (X(N_Orbitals, N_Orbitals))
    ALLOCATE (Correction(N_Orbitals, N_Orbitals))
    ALLOCATE (I(N_Orbitals, N_Orbitals))
  END IF

  I = CMPLX(0.D0)
  DO N = 1, N_Orbitals 
    I(N,N) = CMPLX(1.D0)
  END DO

  A = MATMUL(CONJG(NZCoeffs),TRANSPOSE(NZCoeffs))
  B = MATMUL(CONJG(OldNZCoeffs),TRANSPOSE(NZCoeffs))

  X = 0.5D0*(I-A)

  DO 
    Correction = I-A - MATMUL(TRANSPOSE(CONJG(B)), X) - MATMUL(X, B) - &
               MATMUL(X,X)
    X = X + 0.5*Correction
    Norm = SUM(CONJG(Correction)*Correction)
    IF (Norm<1.D-10) EXIT  
  END DO 
  NZCoeffs = NZCoeffs + MATMUL(CONJG(X),OldNZCoeffs)
  

END SUBROUTINE Rattle