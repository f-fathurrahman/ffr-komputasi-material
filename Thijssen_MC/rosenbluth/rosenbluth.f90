PROGRAM Rosenbluth
! Rosenbluth algorithm for growing polymers according to the canonical distribution
! See section 10.6.1 of `Computational Physics', by J.M. Thijssen
! Program written by Jos Thijssen, 2006-2007
! The present program does not produce any output, but it is provided as
! an example which should enable the reader to read the program
! recur.f90 more easily...

  IMPLICIT NONE
  INTEGER, PARAMETER :: N=73, StepNum = 30000, ThetaNum = 6
! ThetaNum is the number of different theta-values
  DOUBLE PRECISION :: Polymer(N, 2), PolymerWeight, TotWeight
  DOUBLE PRECISION :: sigma, epsilon, sig6, sig12, End2End, PI
! End2End is the end-to-end length of the polymer.

  CALL Initialise
  CALL Simulate
  print *, End2End/TotWeight

CONTAINS

  SUBROUTINE Initialise
    CALL InitRand(810371231)
    End2End = 0.D0
    PI = 4.d0*DATAN(1.d0)
    Sigma = 0.8D0
    Epsilon = 0.25D0
    Sig6 = Sigma**6
    Sig12 = Sigma**12
    TotWeight = 0.D0
  END SUBROUTINE Initialise



  SUBROUTINE Simulate
! Main simulation routine
    INTEGER :: Step
    Polymer(1,:) = (/0.d0, 0.d0/)
    Polymer(2,:) = (/1.d0, 0.d0/)
    DO Step = 1, StepNum
      CALL GrowPolymer()
      CALL AnalysePolymer
    END DO
  END SUBROUTINE Simulate


  SUBROUTINE AnalysePolymer
! Calculates the end-to-end length of the polymer and this contributes 
! to a statistical average of this length
    DOUBLE PRECISION :: Dist, RelPos(2)

    RelPos = Polymer(1,:) - Polymer(N,:)
    Dist = SUM(RelPos*RelPos)
    End2End = End2End + Dist*PolymerWeight
    TotWeight = TotWeight + PolymerWeight
  END SUBROUTINE AnalysePolymer

  SUBROUTINE GrowPolymer
! The actual Rosenbluth algorithm is carried out
    INTEGER :: J, K
    DOUBLE PRECISION :: ProbArr(ThetaNum), Theta0, Theta, &
                        ProbSum, R, TProb, NewPos(2), RealRand, VraiPoids

    PolymerWeight =1.D0  
    DO J=2, N-1
      Theta0 = RealRand()
      ProbSum = 0.D0
      DO K=1, ThetaNum
        Theta = Theta0 + 2*PI*(K-1)/DBLE(ThetaNum)
        NewPos = Polymer(J,:) + (/COS(Theta), SIN(Theta)/)
        ProbArr(K) = EXP(-Potential(J,Polymer,NewPos))
        ProbSum = ProbSum + ProbArr(K)
      END DO
      ProbArr = ProbArr/ProbSum
      R = RealRand()
      TProb = 0.D0
      DO K=1, ThetaNum
        TProb = TProb + ProbArr(K)
        IF (R<TProb) EXIT
      END DO
      Theta = Theta0 + 2*PI*(K-1)/DBLE(ThetaNum)
      Polymer(J+1,:) = Polymer(J,:)+(/COS(Theta), SIN(Theta)/)
      PolymerWeight =PolymerWeight*ProbSum/(0.5D0*ThetaNum)
!      print *, thetanum, probarr(K)
    END DO
  END SUBROUTINE GrowPolymer
  
  
  DOUBLE PRECISION FUNCTION Potential(J, Polymer, NewPos)
! Lennard-Jones interaction potential for the polymer. Contains 
! only the interaction between the last added segment and the remainder.
    INTEGER :: J, K
    DOUBLE PRECISION :: Polymer(N, 2), NewPos(2), RelPos(2), Dist

    Potential = 0.D0
    DO K=1, J-1
      RelPos = Polymer(K,:)-NewPos(:)
      Dist = 1.d0/SUM(RelPos*RelPos)
      Potential = Potential + 4*Epsilon*(-sig6*Dist**3+Sig12*Dist**6)
    END DO
  END FUNCTION Potential
END PROGRAM Rosenbluth
