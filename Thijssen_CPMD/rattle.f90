!----------------------------------------------------
SUBROUTINE Rattle(NZCoeffs, OldNZCoeffs, NZCoeffsdot)
!----------------------------------------------------
  USE globals
  IMPLICIT NONE
  COMPLEX(8), INTENT(INOUT) :: NZCoeffs(N_Orbitals, NoOfPW), &
                               NZCoeffsdot(N_Orbitals, NoOfPW), &
                               OldNZCoeffs(N_Orbitals, NoOfPW)
  ! XXX: NZCoeffsdot is not used?
  INTEGER :: N
  REAL(8) :: nrm
  COMPLEX(8), ALLOCATABLE :: I(:,:), A(:,:), B(:,:), &
                             X(:,:), Correction(:,:)

  IF( .NOT. ALLOCATED(A) ) THEN
    ALLOCATE( A(N_Orbitals, N_Orbitals) )
    ALLOCATE( B(N_Orbitals, N_Orbitals) )
    ALLOCATE( X(N_Orbitals, N_Orbitals) )
    ALLOCATE( Correction(N_Orbitals, N_Orbitals) )
    ALLOCATE( I(N_Orbitals, N_Orbitals) )
  ENDIF

  I = CMPLX(0.D0, kind=8)
  DO N = 1, N_Orbitals 
    I(N,N) = CMPLX(1.D0, kind=8)
  ENDDO

  A = MATMUL(CONJG(NZCoeffs),TRANSPOSE(NZCoeffs))
  B = MATMUL(CONJG(OldNZCoeffs),TRANSPOSE(NZCoeffs))

  ! Initial value
  X = 0.5D0*(I - A)

  DO 
    Correction = I - A - MATMUL(TRANSPOSE(CONJG(B)), X) - MATMUL(X, B) - MATMUL(X,X)
    X = X + 0.5*Correction
    nrm = REAL(SUM( CONJG(Correction)*Correction ), kind=8)
    IF( nrm < 1.D-10 ) EXIT  
  ENDDO 
  NZCoeffs = NZCoeffs + MATMUL(CONJG(X),OldNZCoeffs)

END SUBROUTINE

