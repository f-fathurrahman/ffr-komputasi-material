!*************************************************************
!****                  Calc_Orbitals                      ****
!*************************************************************
! Solve 'Equation of motion' for the wavefunction coefficients
! using velocity-Verlet in combination with Rattle Algorithm. 
! Ions are restricted not to move!
  
!-------------------------
SUBROUTINE Calc_Orbitals()
!-------------------------
  use globals
  IMPLICIT NONE
  INTEGER                     :: Iter
  COMPLEX(8), ALLOCATABLE :: NZCoeffs(:,:), NZCoeffsDot(:,:), &
                             OrbForce(:,:), OldNZCoeffs(:,:), &
                             Y(:,:), R_IonDot(:,:)
  COMPLEX(8) :: E, Eold
  REAL(8) :: Time, TimeStep

  ALLOCATE(OrbForce(N_Orbitals, NoOfPW))
  ALLOCATE(NZCoeffs(N_Orbitals, NoOfPW))
  ALLOCATE(NZCoeffsDot(N_Orbitals, NoOfPW))
  ALLOCATE(R_ionDot(N_ion, 3))
  ALLOCATE(OldNZCoeffs(N_Orbitals, NoOfPW))
  ALLOCATE(Y(N_Orbitals, N_Orbitals))

  NZCoeffsDot = CMPLX(0.D0)
  CALL Init_Coeffs(NZCoeffs,NZCoeffsDot,R_ionDot, Time)
  TimeStep = TimeStepOrt
  OldNZCoeffs = NZCoeffs
  CALL Rattle(NZCoeffs, OldNZCoeffs, NZCoeffsDot)
  CALL Calc_Orb_Force(NZCoeffs, OrbForce)
  
  PrintOut = .TRUE.
  Eold = CMPLX(1.D0)
  E = Eold
  
  DO Iter = 1, MaxIter
    OldNZCoeffs=NZCoeffs
    NZCoeffsDot=NZCoeffsDot+TimeStep*OrbForce/2
    NZCoeffs = NZCoeffs + TimeStep*NZCoeffsDot
    CALL Rattle(NZCoeffs, OldNZCoeffs, NZCoeffsDot)
    CALL Calc_Orb_Force(NZCoeffs, OrbForce)
    NZCoeffsDot=NZCoeffsDot+TimeStep*OrbForce/2
    Y =  MATMUL(CONJG(NZCoeffsDot),TRANSPOSE(NZCoeffs))
    Y = -0.5D0*(CONJG(Y) + TRANSPOSE(Y))
    NZCoeffsdot = NZCoeffsdot + MATMUL(Y,NZCoeffs)
    NZCoeffsdot = RedFac*NZCoeffsdot
    !
    IF(PrintOut) THEN
      Eold = E
      CALL Total_Energy(NZCoeffs, E)
      PrintOut = .FALSE.
    ENDIF
    !
    IF( MOD(Iter,10) == 0 ) PrintOut = .TRUE.
    !
    IF( ABS(Eold-E) < Ediff ) THEN
      CALL Store_Optimal(NZCoeffs, NZCoeffsDot)
      OPEN(12, FILE='Energy.dat', POSITION='APPEND')
      WRITE(12, '(1x,2F18.10)') Ions(1)%R_I(1), DBLE(E)
      CLOSE(12)
      EXIT
    ELSEIF(PrintOut) THEN
      WRITE(*,*) 'Still converging orbitals before starting dynamics...'
    ENDIF  
  ENDDO
END SUBROUTINE Calc_Orbitals