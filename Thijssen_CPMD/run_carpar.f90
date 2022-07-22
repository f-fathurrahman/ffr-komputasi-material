!*************************************************************
!****                  Run_CarPar                         ****
!*************************************************************
! Solve "Equations of Motion" for the wavefunction coefficients
! using velocity-Verlet in combination with the RATTLE algorithm
! The ions are allowed to move now
SUBROUTINE Run_CarPar()
  use globals
  IMPLICIT NONE
  INTEGER                     :: N, Iter 
  complex(8), ALLOCATABLE :: NZCoeffs(:,:), NZCoeffsDot(:,:),&
                                 OrbForce(:,:), OldNZCoeffs(:,:), &
                                 Y(:,:), IonForce(:,:), R_ionDot(:,:)
  COMPLEX(8) :: E, Eold
  REAL(8) :: Time, TimeStep
  ALLOCATE (OrbForce(N_orbitals, NoOfPW))
  ALLOCATE (IonForce(N_ion, 3))
  ALLOCATE (NZCoeffs(N_orbitals, NoOfPW))
  ALLOCATE (NZCoeffsDot(N_orbitals, NoOfPW))
  ALLOCATE (R_ionDot(N_ion, 3))
  ALLOCATE (OldNZCoeffs(N_orbitals, NoOfPW))
  ALLOCATE (Y(N_orbitals, N_orbitals))

  TimeStep = TimeStepCP
  NZCoeffsDot = CMPLX(0.D0)
  OldNZCoeffs = NZCoeffs 
  CALL init_coeffs(NZCoeffs, NZCoeffsDot, R_ionDot,  Time)
  CALL rattle(NZCoeffs, OldNZCoeffs, NZCoeffsDot)
  CALL calc_orb_Force(NZCoeffs, OrbForce)
  R_ionDot = CMPLX(0.D0)
  IonForce = CMPLX(0.D0) 
  PrintOut = .TRUE.
  CALL Calc_Ion_Force(NZCoeffs, IonForce)
  DO Iter = 1, MaxIter
    DO N = 1, N_ion
      R_ionDot(N,:) = R_ionDot(N,:) + &
                   TimeStep*IonForce(N,:)/(2*Ions(N)%Mass)
      Ions(N)%R_I(:) = Ions(N)%R_I(:) + TimeStep*R_ionDot(N,:)
      CALL Periodic_Boundary()
      CALL Fill_Dyn_Grids()
    END DO !N
    Time = Time + TimeStep
    OldNZCoeffs=NZCoeffs
    NZCoeffsDot=NZCoeffsDot+TimeStep*OrbForce/(2*mu)
    NZCoeffs = NZCoeffs + TimeStep*NZCoeffsDot
    CALL Rattle(NZCoeffs, OldNZCoeffs, NZCoeffsDot)
    CALL Calc_Ion_Force(NZCoeffs, IonForce)
    CALL Calc_Orb_Force(NZCoeffs, OrbForce)
    DO N = 1, N_ion
        R_ionDot(N,:) = R_ionDot(N,:) + &
                       TimeStep*IonForce(N,:)/(2*Ions(N)%Mass)
    END DO !N
    NZCoeffsDot=NZCoeffsDot+TimeStep*OrbForce/(2*mu)
    Y =  MATMUL(CONJG(NZCoeffsDot),TRANSPOSE(NZCoeffs))
    Y = -0.5D0*(CONJG(Y) + TRANSPOSE(Y))
    NZCoeffsdot = NZCoeffsdot + MATMUL(Y,NZCoeffs)
    IF (PrintOut) THEN
      PRINT '(A23 F15.3)', 'Time:', Time
      CALL Total_Energy(NZCoeffs, E)
      CALL Check_Const_Energy (NZCoeffsdot, R_ionDot)
      CALL Print_Ions()
      PrintOut = .FALSE.
    END IF

    IF (MOD(Iter,10)==0) THEN
      PrintOut = .TRUE.
      CALL Store_Last(NZCoeffs,NZCoeffsDot, R_IonDot, Time)
    END IF
  END DO !Iter
    
END SUBROUTINE