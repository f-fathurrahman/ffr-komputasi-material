PROGRAM CarPar

USE globals
USE utilities
USE pseudo
USE grids
USE force
USE excorr

CALL InitParams()
CALL Simulate()


CONTAINS 

!*************************************************************
!****                  InitParams                         ****
!*************************************************************
!*** Parameters are read from file InCP and initialised   ****
!*************************************************************

SUBROUTINE InitParams()
  IMPLICIT NONE
  REAL(8) :: Energy_CutOff
  INTEGER :: I, I2, I3, I5, N

  PI = 4.D0*ATAN(1.D0)
  Im = CMPLX(0.D0,1.D0)
    
    OPEN (8, File="InCP")
    READ (8, *) OnlyStatic
    READ (8, *) BoxL
    Omega = BoxL**3
    READ (8, *) Energy_CutOff
    GMax = SQRT(2.D0*Energy_Cutoff)
    GMax = GMax*BoxL*0.5D0/PI
    I2 = Ceiling(log(4*Gmax)/Log(2.D0))
    I2 = 2**I2
    I3 = Ceiling(log(4*Gmax)/Log(3.D0))
    I3 = 3**I3
    I5 = Ceiling(log(4*Gmax)/Log(5.D0))
    I5 = 5**I5
    GridSize = MIN(I2,I3,I5)
    MaxBas = INT(GMax)
    ALLOCATE (Density_K(0:GridSize-1,0:GridSize-1,0:GridSize-1))
    ALLOCATE (Density_R(0:GridSize-1,0:GridSize-1,0:GridSize-1))
    CALL FFT_Plan(GridSize)
    print *, 'Per direction, the linear index of PWs runs from' 
    print *, 0, '  to ', GridSize-1
    READ (8,*) mu
    READ (8,*) TimeStepOrt      
    READ (8,*) TimeStepCP
    READ (8,*) Ediff
    READ (8,*) MaxIter
    READ (8,*) RedFac
    READ (8, *) No_OF_DIFF_IONS
    print *, 'there are ',  No_OF_DIFF_IONS, ' different ions'
    ALLOCATE (PP_Params(No_OF_DIFF_IONS))
    PP_Params(:)%AtomNum = 0
    READ (8, *) N_ion
    ALLOCATE (Ions(N_ion))
    DO I=1, N_ion
      READ (8,*) Ions(I)%AtomNum, &  ! Atomic Number
                 Ions(I)%Mass, &     ! Mass of Ion
                 Ions(I)%R_I(1), Ions(I)%R_I(2), Ions(I)%R_I(3) ! Positions
      Ions(I)%Mass = Ions(I)%Mass*1836.D0
      CALL Init_PP_Params(Ions(I)%AtomNum)
    END DO
    print *, 'No of ions', N_ion
    print *, 'Ion data', Ions  
    CALL InitGrids()

    WRITE(*,*) 'After InitGrids()'

    READ(8, *) N_electron
    READ(8, *) N_Orbitals
    ALLOCATE(FillFac(N_orbitals))
    DO N = 1, N_orbitals
      READ(8,*) FillFac(N)
    ENDDO
    CLOSE(8)

    WRITE(*,*) 'End of InitParams'

  END SUBROUTINE InitParams

 
!*************************************************************
!****                  Simulate                           ****
!*************************************************************
  SUBROUTINE Simulate()
    IMPLICIT NONE
    CALL Calc_Orbitals()
    IF( .NOT. OnlyStatic ) CALL Run_CarPar()
  END SUBROUTINE Simulate

!*************************************************************
!****                  Calc_Orbitals                      ****
!*************************************************************
  ! Solve 'Equation of motion' for the wavefunction coefficients
  ! using velocity-Verlet in combination with Rattle Algorithm. 
  ! Ions are restricted not to move!
  
!-------------------------
SUBROUTINE Calc_Orbitals()
!-------------------------
  IMPLICIT NONE
  INTEGER                     :: Iter
  LOGICAL                     :: PrintOut
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
  CALL InitCoeffs(NZCoeffs,NZCoeffsDot,R_ionDot, Time)
  TimeStep = TimeStepOrt
  OldNZCoeffs = NZCoeffs
  CALL Rattle(NZCoeffs, OldNZCoeffs, NZCoeffsDot)
  CALL Calc_OrbForce(NZCoeffs, OrbForce)
  
  PrintOut = .TRUE.
  Eold = CMPLX(1.D0)
  E = Eold
  
  DO Iter = 1, MaxIter
    OldNZCoeffs=NZCoeffs
    NZCoeffsDot=NZCoeffsDot+TimeStep*OrbForce/2
    NZCoeffs = NZCoeffs + TimeStep*NZCoeffsDot
    CALL Rattle(NZCoeffs, OldNZCoeffs, NZCoeffsDot)
    CALL Calc_OrbForce(NZCoeffs, OrbForce)
    NZCoeffsDot=NZCoeffsDot+TimeStep*OrbForce/2
    Y =  MATMUL(CONJG(NZCoeffsDot),TRANSPOSE(NZCoeffs))
    Y = -0.5D0*(CONJG(Y) + TRANSPOSE(Y))
    NZCoeffsdot = NZCoeffsdot + MATMUL(Y,NZCoeffs)
    NZCoeffsdot = RedFac*NZCoeffsdot
    IF (PrintOut) THEN
       Eold = E
       CALL Total_Energy(NZCoeffs, E)
       PrintOut = .FALSE.
    END IF
    IF (MOD(Iter,10)==0) PrintOut = .TRUE.
    IF (ABS(Eold-E)<Ediff) THEN
        CALL StoreOptimal(NZCoeffs, NZCoeffsDot)
        OPEN (12, FILE='Energy.dat', POSITION='APPEND')
        write(12, *) Ions(1)%R_I(1), DBLE(E)
        CLOSE (12)
        EXIT
    ELSE IF (PrintOut) THEN
      print * , 'Still converging orbitals before starting dynamics...'
    END IF  
 END DO
 END SUBROUTINE Calc_Orbitals
  
  
  
  
!*************************************************************
!****                  Run_CarPar                         ****
!*************************************************************
! Solve "Equations of Motion" for the wavefunction coefficients
! using velocity-Verlet in combination with the RATTLE algorithm
! The ions are allowed to move now


  SUBROUTINE Run_CarPar()

  IMPLICIT NONE
  INTEGER                     :: N, Iter 
  DOUBLE COMPLEX, ALLOCATABLE :: NZCoeffs(:,:), NZCoeffsDot(:,:),&
                                 OrbForce(:,:), OldNZCoeffs(:,:), &
                                 Y(:,:), IonForce(:,:), R_ionDot(:,:)
  DOUBLE COMPLEX              :: E, Eold
  DOUBLE PRECISION            :: Time, TimeStep
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
  CALL InitCoeffs(NZCoeffs, NZCoeffsDot, R_ionDot,  Time)
  CALL Rattle(NZCoeffs, OldNZCoeffs, NZCoeffsDot)
  CALL Calc_OrbForce(NZCoeffs, OrbForce)
  R_ionDot = CMPLX(0.D0)
  IonForce = CMPLX(0.D0) 
  PrintOut = .TRUE.
  CALL Calc_IonForce(NZCoeffs, IonForce)
  DO Iter = 1, MaxIter
    DO N = 1, N_ion
      R_ionDot(N,:) = R_ionDot(N,:) + &
                   TimeStep*IonForce(N,:)/(2*Ions(N)%Mass)
      Ions(N)%R_I(:) = Ions(N)%R_I(:) + TimeStep*R_ionDot(N,:)
      CALL PeriodicBoundary()
      CALL FillDynGrids()
    END DO !N
    Time = Time + TimeStep
    OldNZCoeffs=NZCoeffs
    NZCoeffsDot=NZCoeffsDot+TimeStep*OrbForce/(2*mu)
    NZCoeffs = NZCoeffs + TimeStep*NZCoeffsDot
    CALL Rattle(NZCoeffs, OldNZCoeffs, NZCoeffsDot)
    CALL Calc_IonForce(NZCoeffs, IonForce)
    CALL Calc_OrbForce(NZCoeffs, OrbForce)
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
      CALL Check_ConstEnergy (NZCoeffsdot, R_ionDot)
      CALL PrintIons 
      PrintOut = .FALSE.
    END IF

    IF (MOD(Iter,10)==0) THEN
      PrintOut = .TRUE.
      CALL StoreLast(NZCoeffs,NZCoeffsDot, R_IonDot, Time)
    END IF
  END DO !Iter
    
  END SUBROUTINE Run_CarPar
  



  
!*************************************************************
!****                  InitCoeffs                         ****
!*************************************************************
  
  SUBROUTINE InitCoeffs(NZCoeffs, NZCoeffsDot, R_ionDot, Time)
  IMPLICIT NONE
  LOGICAL                      :: Opt_Exists, Mov_Exists
  DOUBLE COMPLEX, INTENT(OUT)  :: NZCoeffs(N_orbitals, NoOfPW)
  DOUBLE COMPLEX, INTENT(OUT)  :: NZCoeffsDot(N_orbitals, NoOfPW)
  DOUBLE COMPLEX, INTENT(OUT)  :: R_ionDot(N_ion, 3)
  DOUBLE PRECISION, INTENT(OUT):: Time
  
  INQUIRE(FILE='last.dat', EXIST=Mov_Exists)
  INQUIRE(FILE='opt.dat',EXIST=Opt_Exists)
  IF (Mov_Exists) THEN
     CALL GetLast(NZCoeffs, NZCoeffsDot, R_ionDot,  Time)
  ELSE IF (Opt_Exists) THEN
     CALL GetOptimal(NZCoeffs, NZCoeffsDot)
     Time = 0.D0
  ELSE
     CALL InitSol(NZCoeffs, N_Orbitals, NoOfPW)
  END IF
  END SUBROUTINE InitCoeffs
  




  
!*************************************************************
!****                  InitSol                            ****
!*************************************************************
! Starting solution, Gaussian distribution
  
  SUBROUTINE InitSol(NZCoeffs, N_orbitals, NoOfPW)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NoOfPW, N_orbitals
  DOUBLE COMPLEX, INTENT(OUT) :: NZCoeffs(N_orbitals, NoOfPW)
  INTEGER :: Iorb, IIndex
  DOUBLE PRECISION :: Norm, Alpha=1.5D0, G2, X, Y
  
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
  
  END SUBROUTINE InitSol
          
          
          

END PROGRAM CarPar
    
