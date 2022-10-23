      PROGRAM HydrogenVMC
C*****************************************************************
C* This program estimates the ground state energy of the         *
C* hydrogen atom using a variational quantum Monte Carlo         *
C* procedure.                                                    *
C*                                                               *
C* Program written by Jos Thijssen	                         *
C* Spring 1999                                                   *
C* This program is described in section 12.2.2 of the book       *
C* "Computational Physics" by Jos Thijssen,                      *
C* Cambridge University Press 1999                               *
C*****************************************************************

      DOUBLE PRECISION Alpha, Energy

      INTEGER InitStep, WalkNum, StepNum

      CALL InitWalkers(Alpha, InitStep, StepNum, WalkNum)

      CALL Metropolis(Energy, Alpha, InitStep, StepNum, WalkNum)
      print *, Energy
      
      END



      SUBROUTINE InitWalkers(Alpha, InitStep, StepNum, WalkNum)
C Initialise the parameters of the procedure and initialse the 
C 3D random walkers to random initial positions. 

      include "globHydr"

      INTEGER K, L, InitStep, StepNum, WalkNum

      DOUBLE PRECISION RealRand, Alpha

      print *, 'Give nr of steps and nr of equilibration steps'
      read *, StepNum, InitStep
      print *, 'Give nr of Walkers'
      read *, WalkNum
      IF (WalkNum.GT.MaxWalkNum) THEN
        print *, 'This number should be smaller than', MaxWalkNum
        STOP
      END IF 
      print *, 'Give exponential parameter alpha'
      read *, Alpha
      CALL InitRand(4537)
      DO L=1, 3
        DO K=1, WalkNum
          Walker(K,L) = RealRand()-0.5D0
        END DO
      ENDDO
      OPEN (12, File='energy.dat')
      OPEN (8, File='ground_state.dat')
      END



      DOUBLE PRECISION FUNCTION TransProb(R, NewR, Alpha)
C Calculate the transition probability for a walker located
C at Pos to step to NewPos, with variational parameter Alpha

      include "globHydr"

      DOUBLE PRECISION Alpha, R, NewR

      TransProb = EXP(-2*Alpha*(NewR-R))

      END


      SUBROUTINE CalcLocal (ELocal, R, Alpha)
C Calculate the local energy for an electron at distance R 
C from the nucleus. 

      include "globHydr"

      DOUBLE PRECISION ELocal, Alpha, R, ETemp

      ETemp = -1/R-0.5D0*alpha*(alpha-2.d0/R)
      ELocal = ETemp
      END
      


      SUBROUTINE Metropolis (Energy, Alpha, InitStep, StepNum, WalkNum)
C The metropolis procedure. StepNum steps are carried out. In one step,
C all WalkNum walkers are moved to new positions. The first InitStep
C steps are used for equilibration. The energy is evaluated for the 
C variational parameter Alpha.  


      include "globHydr"

      DOUBLE PRECISION Energy, Alpha, Pos(3), NewPos(3), TransProb, 
     .                 ELocal, MaxStep, RealRand, R, NewR

      INTEGER K, L, Step, NAcc, I, Hist(0:40), NetStep, 
     .        InitStep, StepNum, WalkNum
      
      Energy = 0.D0
      NAcc = 0
      MaxStep = 1.0D0
      DO I=0,40
        Hist(I) = 0
      END DO
      DO Step = 1, StepNum
        DO K=1, WalkNum
          R = 0.D0
          NewR = 0.D0
          DO L=1, 3
            Pos(L)    = Walker(K,L)
            NewPos(L) = Pos(L) + (RealRand()-0.5D0)*MaxStep
            R = R + Pos(L)**2
            NewR = NewR + NewPos(L)**2
          END DO
          R = SQRT(R)
          NewR = SQRT(NewR)
          IF (RealRand().LT.TransProb(R, NewR, Alpha)) THEN
            NAcc = NAcc+1
            DO L=1, 3
              Walker(K, L) = NewPos(L)
              R = NewR
            END DO
          END IF
          IF (Step.GT.InitStep) THEN
            CALL CalcLocal(ELocal, R, Alpha)
            Energy = Energy + ELocal
            IF (MOD(Step, 50).EQ.0) WRITE (12,'(F12.6)') ELocal
            IF (R.LT.4.D0) THEN
              Hist(NINT(R*10.D0)) = Hist(NINT(R*10.D0))+1
            END IF
          END IF
        ENDDO
C Dynamic adaptation of the width of the Gaussian distribution
        IF (MOD(Step,100) .EQ. 0) THEN
          MaxStep = MaxStep*NAcc/(50.D0*WalkNum)
          NAcc = 0
        END IF
      ENDDO
      NetStep = StepNum-InitStep
      Energy = Energy/(NetStep*WalkNum)
      DO I=0,40
        write (8, '(F10.5 F10.5)') I*0.1D0, 
     .                             Hist(I)/DBLE(WalkNum*NetStep)
      ENDDO
      END
      

