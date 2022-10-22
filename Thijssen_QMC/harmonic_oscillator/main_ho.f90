PROGRAM HarmOscVMC
! *****************************************************************
! * This program estimates the ground state energy of the         *
! * harmonic oscillator using a variational quantum Monta Carlo   *
! * procedure.                                                    *
! *                                                               *
! * Program written by Jos Thijssen	                         *
! * Spring 1999  february 2006                                    *
! * This program is described in section 12.2.2 of the book       *
! * "Computational Physics" by Jos Thijssen,                      *
! * Cambridge University Press 2007                               *
! *****************************************************************

  IMPLICIT NONE 
  REAL(8) :: Alpha, Energy
  INTEGER :: InitStep, WalkNum, StepNum

  CALL InitWalkers(Alpha, InitStep, StepNum, WalkNum)

  CALL Metropolis(Energy, Alpha, InitStep, StepNum, WalkNum)
  
  WRITE(*,*) 'Energy = ', Energy
  WRITE(*,*) 'analytical energy = ', 0.5*alpha + 1/(8*alpha)
  WRITE(*,*) 'var = ', (1-4*alpha**2)**2/(32*alpha**2)
  
  CLOSE(11)
  CLOSE(12)
  CLOSE(8)
END PROGRAM



SUBROUTINE InitWalkers(Alpha, InitStep, StepNum, WalkNum)
  ! Initialise the parameters of the procedure and initialse the 
  ! random walkers to random initial positions. 
  USE m_global_ho
  IMPLICIT NONE 
  INTEGER :: K, InitStep, StepNum, WalkNum
  REAL(8) :: RealRand, Alpha

  WRITE(*,*) 'Give nr of steps and nr of equilibration steps'
  READ(*,*) StepNum, InitStep
  
  WRITE(*,*) 'Give nr of Walkers'
  READ(*,*) WalkNum
  
  IF( WalkNum > MaxWalkNum) THEN
    WRITE(*,*) 'This number should be smaller than', MaxWalkNum
    STOP
  ENDIF 
      
  WRITE(*,*) 'Give exponential parameter alpha'
  READ(*,*) Alpha
  
  CALL InitRand(4537)
  DO K = 1, WalkNum
    Walker(K) = RealRand() - 0.5D0
  ENDDO
  
  OPEN(11, File='energy.dat')
  OPEN(12, File='variance.dat')
  OPEN(8, File='ground_state.dat')
END SUBROUTINE 


real(8) FUNCTION TransProb(X, NewX, Alpha)
! Calculate the transition probability for a walker with coordinate 
! X stepping to NewX, with variational parameter Alpha
  use m_global_ho
  implicit none 
  real(8) :: Alpha, X, NewX

  TransProb = EXP(-2*Alpha*(NewX**2-X**2))
END FUNCTION


SUBROUTINE CalcLocal(ELocal, X, Alpha)
  ! Calculate the local energy for a walker at position X
  use m_global_ho
  implicit none
  real(8) :: ELocal, Alpha, X
  ELocal = Alpha + X**2*(0.5D0 - 2.D0*Alpha**2)
END SUBROUTINE
      


SUBROUTINE Metropolis(Energy, Alpha, InitStep, StepNum, WalkNum)
  ! The metropolis procedure. StepNum steps are carried out. In one step,
  ! all WalkNum walkers are moved to new positions. The first InitStep
  ! steps are used for equilibration. The energy is evaluated for the 
  ! variational parameter Alpha.  
  use m_global_ho
  implicit none

  real(8) :: Energy, Alpha, X, NewX, TransProb, ELocal, &
             GaussWidth, Rand1, Rand2, RandArr(MaxWalkNum), &
             RealRand, EnerSq

  INTEGER :: K, Step, NAcc, I, Hist(-100:100), NetStep, &
             InitStep, StepNum, WalkNum
      
  Energy = 0.D0
  EnerSq = 0.D0
  NAcc = 0
  GaussWidth = 1.D0
  
  DO I = -100,100
    Hist(I) = 0
  ENDDO
  
  DO Step = 1, StepNum
    DO K=1, WalkNum/2
      CALL ExpRand(Rand1, Rand2)
      RandArr(2*K-1) = Rand1
      RandArr(2*K)   = Rand2
    END DO
    Energy = 0.D0
    EnerSq = 0.D0
    DO K=1, WalkNum
      X    = Walker(K)
      NewX = Walker(K) + RandArr(K)*GaussWidth
      ! Again, the expression for the transition probability is simple as a 
      ! result of the Gaussian form of the trial function
      IF (RealRand() < TransProb(X, NewX, Alpha)) THEN
        NAcc = NAcc+1
        Walker(K) = NewX
        X = NewX
      ENDIF
      !
      IF(Step > InitStep) THEN
        CALL CalcLocal(ELocal, X, Alpha)
        Energy = Energy + ELocal
        EnerSq = EnerSq + ELocal*ELocal
        IF (ABS(X) < 2.D0) THEN
          Hist(NINT(X*5)) = Hist(NINT(X*5))+1
        ENDIF
      ENDIF
    ENDDO
    !
    ! Dynamic adaptation of the width of the Gaussian distribution
    IF( MOD(Step,100) == 0) THEN
      GaussWidth = GaussWidth*NAcc/(50.D0*WalkNum)
      NAcc = 0
    ENDIF
    !
    IF( Step > InitStep) THEN
      Energy = Energy/WalkNum
      !write(*,*) 'Energy in Metropolis() = ', Energy
      WRITE(11, *) Energy
      WRITE(12, *) EnerSq/WalkNum - Energy**2
      !Energy = 0.d0
      !EnerSq = 0.d0
    END IF
  ENDDO
  !
  NetStep = StepNum-InitStep
  Energy = Energy/(NetStep*WalkNum)
  DO I=-10,10
    write(8, '(F10.5 F10.5)') I*0.2D0, Hist(I)/DBLE(WalkNum*NetStep)
  ENDDO
END SUBROUTINE
      

