      PROGRAM Ising
C This program performs a Monte Carlo simulation for the 2D 
C Ising model on a quadratic lattice. 
C The temperature is included in the coupling strength
C J = K/(kT), where K is the "bare" coupling strength, and
C H = B/(kT), where B is the magnetic field.

      CALL Initialise
      CALL MonteCarlo
      CALL Finalise
      END    


      SUBROUTINE Initialise
      INCLUDE 'globIsMC'
      INTEGER I, IX, IY

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
      DO I=0, 4
        MCWeights(I) = EXP(-2*I*J)
      END DO
#ifdef Plot
        CALL InitPlot('lightblue', 700,700, 'out.ps', 1)
        CALL Framing(-0.05D0*Size, -0.05D0*Size, 
     .                1.05D0*Size,  1.05D0*Size)
        CALL PutStopButton()
#endif
      CALL InitRand(12357)
      END

 


      SUBROUTINE MonteCarlo
      INCLUDE 'globIsMC'

      INTEGER MCIter, IX, IY, CostFactor, CurSpin
      DOUBLE PRECISION RealRand

      OPEN (8, File='Energy')
      OPEN (9, File='Magnet')

      DO MCIter = 1, MCSteps
        DO IX = 0, Size-1
          DO IY = 0, Size-1
            CurSpin = Spins(IX, IY)
            CostFactor = CurSpin*(Spins(MOD(IX+1, Size), IY)+
     .                            Spins(MOD(IX-1+Size, Size), IY)+
     .                            Spins(IX, MOD(IY+1, Size))+
     .                            Spins(IX, MOD(IY-1+Size,Size)))
            IF (CostFactor.GT.0) THEN
              IF (RealRand().LT.MCWeights(CostFactor)) THEN
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
        IF (MCIter.GT.300) THEN
          WRITE (8,'(F12.6)') Energy
          WRITE (9,'(I8)') Magnetisation
        END IF
#ifdef Plot
        IF (MOD(MCIter,500).EQ.0) CALL Display
#endif
      END DO
      END 


      SUBROUTINE Display
      INCLUDE "globIsMC"

      INTEGER IX, IY
#ifdef Plot
      CALL NextPage()
      DO IX=1, Size
        DO IY=1, Size
          IF (Spins(IX,IY) .EQ. 1) THEN
            CALL FillRectangle(DBLE(IX), DBLE(IY), 
     .                         DBLE(IX+1), DBLE(IY+1))
          END IF
        END DO
      END DO
#endif
      END


      SUBROUTINE Finalise
#ifdef Plot
      CALL EndPlot()
#endif
      END
