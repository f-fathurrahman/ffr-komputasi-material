      PROGRAM MC
C*****************************************************************
C* This program performs a Monte Carlo (MC) simulation           *
C* for a monatomic material.                                     *
C* Periodic boundary conditions are assumed. The potential is    *
C* cut off beyond a distance "Rcutoff". Furthermore, the minimum *
C* image convention is used to calculate forces etc.             *
C*                                                               *
C* Global variables are defined in the file "globmc.f"           *
C*                                                               *
C* Program written by Jos Thijssen	                         *
C* Spring 1999                                                   *
C* This program is described in section 10.3.2 of the book       *
C* "Computational Physics" by Jos Thijssen,                      *
C* Cambridge University Press 1999                               *
C*****************************************************************

C**********************Global Variables***************************
       include "globmc.f"                                       
C*****************************************************************

      LOGICAL EquilPhase

      Call Initialise
      Call InitDraw

      EquilPhase = .TRUE.
      Call Simulation (EquilPhase, InitStep)
      EquilPhase = .FALSE.
      Call Simulation (EquilPhase, SimStep)

      CALL WriteCorr
#ifdef Plot
      Call WithDraw
#endif

      END



      SUBROUTINE InitDraw

C**********************Global Variables***************************
       include "globmc.f"                                       
C*****************************************************************
      DOUBLE PRECISION WinBor

#ifdef Plot
      CALL InitPlot ("lightblue", 750, 750, "mcplot", 1)
      WinBor = VolSize+0.1D0
      print *, winbor
      CALL Framing (-0.1d0, -0.1d0, WinBor, WinBor)
      CALL PutStartButton()
      CALL PutStopButton()
#endif
      END


      SUBROUTINE WithDraw

#ifdef Plot
      CALL EndPlot
#endif

      END
 


      SUBROUTINE WriteCorr

C**********************Global Variables***************************
       include "globmc.f"                                       
C*****************************************************************

      INTEGER I

      CALL CalcPotent
   
      DO I=1, 250
        WRITE (11, '(I4, F12.6)') I, dble(CorrArray(I))/dble(I*I)
      ENDDO
C Close open files
      CLOSE(2)
      CLOSE(11)
      CLOSE(7)
      END 
     


      SUBROUTINE Initialise

      CALL InitParameters
      CALL InitPositions
      CALL OpenFiles

      END


      SUBROUTINE OpenFiles

C**********************Global Variables***************************
       include "globmc.f"                                       
C*****************************************************************

      OPEN (2, file='potential')
      OPEN (7, file='virial')
      OPEN (11, file='correl')
      END

      SUBROUTINE InitParameters

C************  User input ****************************************

C**********************Global Variables***************************
      include  "globmc.f"                                            
C*****************************************************************

      INTEGER I, Num

      OPEN (UNIT=8, FILE='mc.in')

      WRITE (6,*) 'Give temperature in reduced units'
      READ (8,*) Temperature
      Beta = 1.d0/Temperature

      WRITE (6,*) 'Give number of equilibration steps'
      READ (8,*) InitStep

      WRITE (6,*) 'Give number of steps between rescaling'
      WRITE (6,*) 'of displacement volume.'
      READ (8,*) ScaleStep

      WRITE (6,*) 'Give number of steps between updating'
      WRITE (6,*) 'the Verlet pair list.'
      READ (8,*) UpdateStep

      WRITE (6,*) 'Give number of production steps'
      READ (8,*) SimStep

      WRITE (6,*) 'Give the density in reduced units'
      READ (8,*) Density

      WRITE (6,*) 'Give the random number generator seed'
      READ (8,*) num

      WRITE (6,*) 'Give number of steps between output '
      WRITE (6,*) 'of data'
      READ  (8,*) OutSteps

#ifdef Plot
      WRITE (6,*) 'Give number of steps between displays'
      READ (8,*) DispInt
#endif

C****** Density and particle number define the volume of the system **
      Volume = DBLE(PartNum/Density)

      CorrStep = 0.02
      DO I=1, MaxCorr
        CorrArray(I) = 0.D0
      ENDDO

      CALL InitRand(num)

      RInner = 2.5D0
      ROuter = 3.3D0

      END




      SUBROUTINE InitPositions

C *** Positions are stored on a regular fcc lattice. Lattice constant* 
C *** is adjusted to volume size to fill the volume homogeneously ****

C**********************Global Variables***************************
      include   "globmc.f"
C*****************************************************************

      INTEGER LinCell, IX, IY, IZ, Teller
      DOUBLE PRECISION  LattConst, Third


C*******Calculate Volume Size from Volume ************************
C*******LinCell is number of cells along one side ****************
      Third = 1.D0/3.D0
      LinCell = NINT ((DBLE(PartNum)/4)**Third)
      VolSize = Volume**Third
      LattConst = VolSize/LinCell

      Teller = 0
      DO IX = 0, LinCell - 1
        DO IY = 0, LinCell - 1
          DO IZ = 0, LinCell - 1
            Teller = Teller + 1
            Qx(Teller) = (Ix+0.25D0)*LattConst
            Qy(Teller) = (Iy+0.25D0)*LattConst
            Qz(Teller) = (Iz+0.25D0)*LattConst
            Teller = Teller + 1
            Qx(Teller) = (Ix+0.75D0)*LattConst
            Qy(Teller) = (Iy+0.75D0)*LattConst
            Qz(Teller) = (Iz+0.25D0)*LattConst
            Teller = Teller + 1
            Qx(Teller) = (Ix+0.75D0)*LattConst
            Qy(Teller) = (Iy+0.25D0)*LattConst
            Qz(Teller) = (Iz+0.75D0)*LattConst
            Teller = Teller + 1
            Qx(Teller) = (Ix+0.25D0)*LattConst
            Qy(Teller) = (Iy+0.75D0)*LattConst
            Qz(Teller) = (Iz+0.75D0)*LattConst
          END DO
        END DO
      END DO
C Maximum displacement:
      Displace = LattConst*0.15D0
      CALL CalcPotent 
      END







      SUBROUTINE Simulation (EquilPhase, StepNum)

C *******************************************************************
C ****  Subroutine in which the actual simulation is performed. *****
C *******************************************************************

C**********************Global Variables***************************
      include   "globmc.f"
C*****************************************************************

      INTEGER StepNum, Step, AcceptRate
      LOGICAL EquilPhase

      CALL UpdatePairList
      AcceptRate = 0

      DO Step = 1, StepNum
        CALL MCSweep(AcceptRate)
C Graphical display:       
#ifdef Plot     
        IF (MOD(Step, DispInt).EQ.0) THEN
          WRITE (6,*) Step
          CALL DrawConfig (PartNum, VolSize, Qx, Qy)
        END IF
#endif
C Data ouput on disk:
        IF (.NOT.EquilPhase.AND.(MOD(Step, OutSteps).EQ.0)) THEN
          CALL OutputMC(Step, EquilPhase)
        END IF
C Scaling displacement volume:
        IF (EquilPhase.AND.(MOD(Step, ScaleStep).EQ.0)) THEN
          CALL ScaleDisplVol(AcceptRate)
          AcceptRate = 0
        END IF
C Updating pair list:
        IF (MOD(Step,UpdateStep).EQ.0) THEN
          CALL UpdatePairList
        END IF
      END DO

      END


      SUBROUTINE ScaleDisplVol(AcceptRate)

C**********************Global Variables***************************
      include   "globmc.f"
C*****************************************************************

      INTEGER AcceptRate

      DOUBLE PRECISION Rate

      Rate = DBLE(AcceptRate)/ScaleStep/PartNum
      IF (Rate.LT.0.5D0) THEN
        Displace = Displace*0.95D0
      ELSE
        Displace = Displace*1.05
      ENDIF
      END


      SUBROUTINE OutputMC(Step, EquilPhase)

C**********************Global Variables***************************
      include   "globmc.f"
C*****************************************************************

      INTEGER Step

      LOGICAl EquilPhase

      CALL CalcPotent 
      WRITE (2,'(F12.5)') Potential
      WRITE (7, '(F12.5)') Virial
C The following ensures that the output is stored in the files even
C if the user interrupts the program
      ENDFILE 7
      ENDFILE 2
      END


      SUBROUTINE DrawConfig (PartNum, VolSize, Qx, Qy)
C Simple graphical display of the configuration

      INTEGER PartNum, I
      DOUBLE PRECISION VolSize, QX(PartNum), QY(PartNum)

C Only the projection onto the XY plane of the particles
C is shown. VolSize is the linear volume size, Qx, Qy
C are the X and Y coordinates of the particles.
C There are PartNum particles
#ifdef Plot
      DO I=1, PartNum
	CALL SetPoint (Qx(I), Qy(I))
      END DO
#endif
      END
      

      INTEGER FUNCTION TINT (X)
C returns the largest integer smaller than X
      
      DOUBLE PRECISION X

      TINT = NINT(INT(X)-0.5D0*(1-X/ABS(X)))
      END
        




      SUBROUTINE MCSweep (AcceptRate)


C *** A sweep over all particles. For each particle, an attempt
C *** is made to move it, and this move is carried out if according
C *** to the Metropolis acceptance criterion.
C *** AcceptRate is the acceptance rate which is used to determine 
C *** the size of the volume inside which the particle is moved.

C**********************Global Variables***************************
      include   "globmc.f"
C*****************************************************************

      DOUBLE PRECISION PotDiff, RealRand, BoltzFac,
     .       NewQx, NewQy, NewQz

      INTEGER I, AcceptRate, TINT

      DO I=1,PartNum
        NewQx = Qx(I) + Displace*(2*RealRand()-1.d0)
        NewQy = Qy(I) + Displace*(2*RealRand()-1.d0)
        NewQz = Qz(I) + Displace*(2*RealRand()-1.d0)
        NewQx = NewQx-TINT(NewQx/VolSize)*VolSize
        NewQy = NewQy-TINT(NewQy/VolSize)*VolSize
        NewQz = NewQz-TINT(NewQz/VolSize)*VolSize
        CALL CalcPotDiff(NewQx, NewQy, NewQz, I, PotDiff)
        BoltzFac = EXP(-Beta*PotDiff)
        IF (BoltzFac.GT.RealRand()) THEN
          Qx(I) = NewQx
          Qy(I) = NewQy
          Qz(I) = NewQz
          AcceptRate = AcceptRate+1
        END IF
      END DO
      CALL CalcPotent
      END

      SUBROUTINE CalcPotDiff(NewQx, NewQy, NewQz, I, PotDiff)
C *** Calculates the potential energy difference PotDiff between the old
C *** configuration Qx, Qy and Qz, and the new one, NewQx etcetera.
C *** Particle number I is the only one which is moved to
C *** (NewQx, NewQy, NewQz)
C *** VirDiff is the difference in the Virial

C**********************Global Variables***************************
      include   "globmc.f"
C*****************************************************************

      DOUBLE PRECISION NewQx, NewQy, NewQz, PotDiff, R2, R4, R6, R12, 
     .       RMin2, NewPot, OldPot

      INTEGER I, J, NeighCnt

      PotDiff = 0.d0
      DO NeighCnt=1, NeighNum(I) 
        J = PairList(I, NeighCnt)
        CALL CalcDistIJ (I, J, R2)
        IF (R2.LT.RInner*RInner) THEN
          RMin2 = 1.d0/R2
          R4 = RMin2*RMin2
          R6 = R4*RMin2
          R12 = R6*R6
          OldPot = R12-2.d0*R6
        ELSE
          OldPot = 0.d0
        ENDIF

        CALL CalcDist (NewQx, NewQy, NewQz, 
     .                 Qx(J), Qy(J), Qz(J), R2)
        IF (R2.LT.RInner*RInner) THEN
          RMin2 = 1.d0/R2
          R4 = RMin2*RMin2
          R6 = R4*RMin2
          R12 = R6*R6
          NewPot = R12-2.d0*R6
        ELSE
          NewPot = 0.d0
        ENDIF
        PotDiff = PotDiff+ NewPot-OldPot
      ENDDO
      END


      SUBROUTINE CalcDist(X1, Y1, Z1, X2, Y2, Z2, R2)

C Calculates the distance between (X1,Y1,Z1) and (X2,Y2,Z2)
C within the periodic boundary conditions
 
C**********************Global Variables***************************
      include   "globmc.f" 
C*****************************************************************

      DOUBLE PRECISION X1, Y1, Z1, X2, Y2, Z2, R2, Dx, Dy, Dz

c      print *, x1, y1, z1, x2, y2, z2

      Dx = X1 - X2
      IF (Dx .GT. VolSize/2) THEN
        Dx = Dx - VolSize
      ENDIF
      IF (Dx .LT. -VolSize/2) THEN
        Dx = Dx+VolSize
      ENDIF

      Dy = Y1 - Y2
      IF (Dy .GT. VolSize/2) THEN
        Dy = Dy - VolSize
      ENDIF
      IF (Dy .LT. -VolSize/2) THEN
        Dy = Dy+VolSize
      ENDIF

      Dz = Z1 - Z2
      IF (Dz .GT. VolSize/2) THEN
        Dz = Dz - VolSize
      ENDIF
      IF (Dz .LT. -VolSize/2) THEN
        Dz = Dz+VolSize
      ENDIF
      R2 = Dx*Dx+Dy*Dy+Dz*Dz
      END

      SUBROUTINE CalcDistIJ(I, J, R2)

C Calculates the distance between particles I and J modulo PBC,
C and returns the square of this distance in the variable R2
 
C**********************Global Variables***************************
      include   "globmc.f" 
C*****************************************************************

      DOUBLE PRECISION Dx, Dy, Dz, R2
      INTEGER I, J

      Dx = Qx(I) - Qx(J)
      IF (Dx .GT. VolSize/2) THEN
        Dx = Dx - VolSize
      ENDIF
      IF (Dx .LT. -VolSize/2) THEN
        Dx = Dx+VolSize
      ENDIF

      Dy = Qy(I) - Qy(J)
      IF (Dy .GT. VolSize/2) THEN
        Dy = Dy - VolSize
      ENDIF
      IF (Dy .LT. -VolSize/2) THEN
        Dy = Dy+VolSize
      ENDIF

      Dz = Qz(I) - Qz(J)
      IF (Dz .GT. VolSize/2) THEN
        Dz = Dz - VolSize
      ENDIF
      IF (Dz .LT. -VolSize/2) THEN
        Dz = Dz+VolSize
      ENDIF
      R2 = Dx*Dx+Dy*Dy+Dz*Dz
      END



      SUBROUTINE CalcPotent

C *** Calculation of total potential energy. We assume that the*******
C *** total potential energy can be written as a a superposition of **
C *** central-symmetric forces between two particles.  *************** 

C**********************Global Variables***************************
      include   "globmc.f" 
C*****************************************************************

      DOUBLE PRECISION R2, R4, R6, R12, 
     .       RMin2   

      INTEGER NeighCnt, I, J

      Potential = 0.D0
      Virial    = 0.D0
      DO I = 1,PartNum
        DO NeighCnt = 1, NeighNum(I)
          J = PairList(I, NeighCnt)
          CALL CalcDistIJ (I, J, R2)
          IF (R2.LT.RInner*RInner) THEN
            RMin2 = 1/R2
            R4 = RMin2*RMin2
            R6 = R4*RMin2
            R12= R6*R6
            Potential = Potential + (R12-2.d0*R6)
            Virial = Virial + 12.d0*(R12-R6)
          ENDIF
        ENDDO
      ENDDO
C The factor 1/2 arises because of the fact that pairs
C have been counted doubly
      Potential = Potential*0.5D0
      Virial = Virial*0.5D0
      END





      SUBROUTINE UpdatePairList

C *** Update of the Verlet pair list.
C *** Only particle pairs with distances smaller than
C *** Router are stored in this list. 
C *** The list is coded in the array PairList. The first index 
C *** of this array runs over the particles, and the second 
C *** is the coordinate of a neighbour of the first index. 
C *** For a particle I, the number of its neighbours is stored in 
C *** the array NeighNum(I)

C**********************Global Variables***************************
      include   "globmc.f" 
C*****************************************************************

      DOUBLE PRECISION R2, Dx, Dy, Dz

      INTEGER I, J

      DO I=1, PartNum
        NeighNum(I) = 0
      ENDDO

      DO I = 1,PartNum
        DO J= I+1,PartNum
           
           Dx = Qx(I) -  Qx(J)
           IF (Dx .GT. VolSize/2) THEN
               Dx = Dx - VolSize
           ENDIF
           IF (Dx .LT. -VolSize/2) THEN
               Dx = Dx+VolSize
           ENDIF

           Dy = Qy(I) -  Qy(J)
           IF (Dy .GT. VolSize/2) THEN
               Dy = Dy - VolSize
           ENDIF
           IF (Dy .LT. -VolSize/2) THEN
               Dy = Dy+VolSize
           ENDIF
           Dz = Qz(I) -  Qz(J)
           IF (Dz .GT. VolSize/2) THEN
               Dz = Dz - VolSize
           ENDIF
           IF (Dz .LT. -VolSize/2) THEN
               Dz = Dz+VolSize
           ENDIF

           R2 = Dx*Dx+Dy*Dy+Dz*Dz
           IF (R2.LT.ROuter*ROuter) THEN
             NeighNum(I) = NeighNum(I)+1
             NeighNum(J) = NeighNum(J)+1
             PairList(I, NeighNum(I)) = J
             PairList(J, NeighNum(J)) = I
           ENDIF
        ENDDO
      ENDDO
      END




