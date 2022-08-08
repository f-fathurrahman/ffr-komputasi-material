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
  use energy_components, only: print_energy_components
  IMPLICIT NONE
  INTEGER  :: Iter
  COMPLEX(8), ALLOCATABLE :: NZCoeffs(:,:), NZCoeffsDot(:,:), &
                             OrbForce(:,:), OldNZCoeffs(:,:), &
                             Y(:,:), R_IonDot(:,:)
  COMPLEX(8) :: E, Eold
  REAL(8) :: Time, TimeStep
  REAL(8) :: delta_E

  ALLOCATE(OrbForce(N_Orbitals, NoOfPW))
  ALLOCATE(NZCoeffs(N_Orbitals, NoOfPW))
  ALLOCATE(NZCoeffsDot(N_Orbitals, NoOfPW))
  ALLOCATE(R_ionDot(N_ion, 3))
  ALLOCATE(OldNZCoeffs(N_Orbitals, NoOfPW))
  ALLOCATE(Y(N_Orbitals, N_Orbitals))

  NZCoeffsDot = CMPLX(0.D0, kind=8)
  !
  CALL init_coeffs(NZCoeffs, NZCoeffsDot, R_ionDot, Time)
  !
  TimeStep = TimeStepOrt
  OldNZCoeffs = NZCoeffs
  !
  WRITE(*,*) 'Before rattle: sum NZCoeffs = ', sum(NZCoeffs)
  CALL rattle(NZCoeffs, OldNZCoeffs, NZCoeffsDot)
  WRITE(*,*) 'After rattle: sum NZCoeffs = ', sum(NZCoeffs)
  !STOP 'ffr 39'
  !
  CALL calc_orb_force(NZCoeffs, OrbForce)
  
  PrintOut = .TRUE.
  Eold = CMPLX(1.D0, kind=8)
  E = Eold
  
  DO Iter = 1, MaxIter
    OldNZCoeffs = NZCoeffs ! save old orbitals
    NZCoeffsDot = NZCoeffsDot + TimeStep*OrbForce/2 ! orbital velocity
    NZCoeffs = NZCoeffs + TimeStep*NZCoeffsDot ! update orbitals
    !
    CALL rattle(NZCoeffs, OldNZCoeffs, NZCoeffsDot) ! apply constraint
    CALL calc_orb_force(NZCoeffs, OrbForce) ! calc gradient
    !
    ! Electron velocities, mass is taken to be one
    NZCoeffsDot = NZCoeffsDot + TimeStep*OrbForce/2
    Y =  MATMUL(CONJG(NZCoeffsDot),TRANSPOSE(NZCoeffs))
    Y = -0.5D0*(CONJG(Y) + TRANSPOSE(Y))
    NZCoeffsdot = NZCoeffsdot + MATMUL(Y,NZCoeffs)
    NZCoeffsdot = RedFac*NZCoeffsdot
    !
    IF(PrintOut) THEN
      Eold = E
      CALL Total_Energy(NZCoeffs, E)
      !PrintOut = .FALSE.
    ENDIF

    delta_E = abs(Eold - E)
    write(*,'(1x,A,I8,F18.10,ES18.10)') 'calc_orbitals: ', Iter, real(E), delta_E

    !
    !IF( MOD(Iter,10) == 0 ) PrintOut = .TRUE.
    !
    
    IF( ABS(Eold-E) < Ediff ) THEN
      !
      write(*,'(1x,A,ES18.10)') 'calc_orbitals: CONVERGED, diffE = ', abs(Eold-E) 
      !
      call print_energy_components()
      !
      CALL Store_Optimal(NZCoeffs, NZCoeffsDot)
      !
      OPEN(12, FILE='Energy.dat', POSITION='APPEND')
      WRITE(12, '(1x,2F18.10)') Ions(1)%R_I(1), DBLE(E)
      CLOSE(12)
      EXIT
    ENDIF

  ENDDO
END SUBROUTINE Calc_Orbitals
