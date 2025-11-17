PROGRAM Ising
! This program performs a Monte Carlo simulation for the 2D 
! Ising model on a quadratic lattice. 
! The temperature is included in the coupling strength
! J = K/(kT), where K is the "bare" coupling strength, and
! H = B/(kT), where B is the magnetic field.

CALL initialize()
CALL MonteCarlo()

END    


SUBROUTINE initialize()
  USE m_globals
  IMPLICIT NONE
  !
  INTEGER :: I, IX, IY

  PRINT *, 'Give reduced coupling K/(kT): '
  READ *, J
  PRINT *, 'Give reduced magnetic field strength B/(kT): '
  READ *, H
  PRINT *, 'Give linear lattice size: '
  READ *, Size
  PRINT *, 'Give number of MC steps: '
  READ *, MCSteps
  
  DO IX = 0, Size-1
    DO IY = 0, Size-1      
      Spins(IX, IY) = 1
    END DO
  END DO
  Magnetisation = Size**2
  Energy = -2.D0*Size*Size*J
  DO I = 0, 4
    MCWeights(I) = EXP(-2*I*J)
  END DO

  CALL InitRand(12357)

END SUBROUTINE

 
SUBROUTINE MonteCarlo()
  use m_globals
  implicit none
  INTEGER :: MCIter, IX, IY, CostFactor, CurSpin
  real(8) :: RealRand

  OPEN(8, File='Energy.dat')
  OPEN(9, File='Magnet.dat')

  DO MCIter = 1, MCSteps
    DO IX = 0, Size-1
      DO IY = 0, Size-1
        CurSpin = Spins(IX, IY)
        CostFactor = CurSpin*( Spins(MOD(IX+1, Size), IY) + &
                               Spins(MOD(IX-1+Size, Size),IY) + &
                               Spins(IX, MOD(IY+1, Size)) + &
                               Spins(IX, MOD(IY-1+Size,Size)))
        IF (CostFactor > 0) THEN
          IF (RealRand() < MCWeights(CostFactor)) THEN
            Spins(IX, IY) = -Spins(IX, IY)
            Magnetisation = Magnetisation + 2*Spins(IX,IY)
            Energy = Energy+2*J*CostFactor
          ENDIF
        ELSE
          Spins(IX, IY) = -Spins(IX, IY)
          Magnetisation = Magnetisation + 2*Spins(IX,IY)
          Energy = Energy+2*J*CostFactor
        ENDIF
      END DO
    END DO
    !
    IF( MCIter > 300) THEN
      WRITE (8,'(F12.6)') Energy
      WRITE (9,'(I8)') Magnetisation
    END IF

  END DO

END 
