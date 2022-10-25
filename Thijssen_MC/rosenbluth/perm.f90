PROGRAM PERM
! Implementation of Grassberger's PERM algorithm, according 
! to section 10.6.1 of `Computational Physics' by J.M. Thijssen, Cambridge 2007 
! The program produces a table containg the size (number of beads) N, 
! the end-to-end length squared, and a statistical error in this length.
! The results should match the form R^2 ~ (N-1)^1.5
! In figure 10.3, the straight line has the form (N-1)^1.5, and not with 
! exponent 0.75 as given in the caption.
! Program written by Jos Thijssen, 2006-2007
 IMPLICIT NONE
  INTEGER, PARAMETER :: N=150, StepNum = 5000, ThetaNum = 9, MaxConf=5*StepNum
  DOUBLE PRECISION :: Polymer(N, 2), PolymerWeight, TotWeight, &
                      sigma, epsilon, sig6, sig12, PI, Weight, &
                      End2End(N, MaxConf), WeightArr(N, MaxConf)
  INTEGER :: StC(N), K, addcount
  
  CALL Initialise
  Polymer(1,:) = (/0.d0, 0.d0/)
  Polymer(2,:) = (/1.d0, 0.d0/)
  DO K=1, StepNum
    Weight = 1.D0
    CALL SimulStep(Polymer, 2, Weight)
  END DO
  CALL PrintResults
  
CONTAINS
  SUBROUTINE Initialise
    
    CALL InitRand(831231)
    End2End = 0.D0
    PI = 4.d0*DATAN(1.d0)
    Sigma = 0.8D0
    Epsilon = 0.25D0
    Sig6 = Sigma**6
    Sig12 = Sigma**12
    TotWeight = 0.D0
    StC = 0
    End2End = 0.D0
    addcount = 0
  END SUBROUTINE Initialise
  
  
  SUBROUTINE PrintResults
! Prints out final results after completing the simulation
  INTEGER,PARAMETER :: BlockNum = 10
  DOUBLE PRECISION :: Average(BlockNum),Weight, AvA, AvSq, ERrA
  INTEGER :: Length, Num, I, Lo, Hi, BlockSize
!  print *, addcount
  DO Length = 3, N
    Num = StC(Length)
    BlockSize = (Num/BlockNum)
    DO I=1, BlockNum-1
      Lo = (I-1)*BlockSize+1
      Hi = I*BlockSize
      Weight = SUM(WeightArr(Length, Lo:Hi)) 
      Average(I) = SUM(End2End(Length, Lo:Hi)*WeightArr(Length, Lo:Hi))/Weight
    END DO
    AvSq = SUM(Average(1:BlockNum-1)*Average(1:BlockNum-1))/(BlockNum-1)
    AvA = SUM(Average(1:BlockNum-1))/(BlockNum-1)
    ErrA = SQRT((AvSq-AvA*AvA)/(BlockNum-1))
    print '(I8, 2F12.4)', Length, AvA, ErrA
  END DO
  END SUBROUTINE PrintResults
 
  RECURSIVE SUBROUTINE SimulStep(Polymer, Length,Weight)
! Main routine, which recursively runs the perm algorithm. Very similar to
! the routine in recur.f90, except for the pruning and enrichment.
    DOUBLE PRECISION, INTENT(INOUT) :: Polymer(N, 2), Weight
    INTEGER, INTENT(IN) :: Length
    INTEGER :: K, CopyNum
    DOUBLE PRECISION :: ProbArr(ThetaNum), Theta0, Theta, &
                        ProbSum, R, TProb, NewPos(2), RealRand, &
                        RelPos(2), Dist, UpLim, LowLim, AvWeight,&
                        ZeroWeight, NewWeight
    
    Theta0 = RealRand()
    ProbSum = 0.D0
    addcount = addcount + 1
    DO K=1, ThetaNum
      Theta = Theta0 + 2*PI*(K-1)/DBLE(ThetaNum)
      NewPos = Polymer(Length,:) + (/COS(Theta), SIN(Theta)/)
      ProbArr(K) = EXP(-Potential(Length,Polymer,NewPos))
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
    Polymer(Length+1,:) = Polymer(Length,:)+(/COS(Theta), SIN(Theta)/)
    Weight = Weight*ProbSum/(0.75*ThetaNum)
    StC(Length+1) = StC(Length+1)+1
    RelPos = Polymer(1,:) - Polymer(Length+1,:)
    Dist = SUM(RelPos*RelPos)
    End2End(Length+1, StC(Length+1)) = Dist
    WeightArr(Length+1, StC(Length+1)) = Weight
    AvWeight = SUM(WeightArr(Length+1,1:StC(Length+1)))!* &
    ZeroWeight = SUM(WeightArr(3,1:StC(3)))
    UpLim = 2.0*AvWeight/ZeroWeight
    LowLim = 1.2*AvWeight/ZeroWeight
    NewWeight = Weight
    IF ((StC(Length+1)>MaxConf-1))THEN !.OR.(StC(Length+1)<=10)) THEN
      IF (Length<N-1) CALL SimulStep(Polymer, Length+1, NewWeight)
    ELSE
      IF (Length<N-1) THEN
        IF (Weight>UpLim) THEN ! Enrich
          NewWeight = Weight*0.5D0
          CALL SimulStep(Polymer, Length+1, NewWeight)
          NewWeight = Weight*0.5D0
          CALL SimulStep(Polymer, Length+1, NewWeight)
        ELSE IF (Weight<LowLim) THEN ! Prune
          NewWeight = 2.d0*Weight
          IF (RealRand()<0.5D0) THEN
             CALL SimulStep(Polymer, Length+1, NewWeight)
          END IF
        ELSE
          CALL SimulStep(Polymer, Length+1, NewWeight)
        END IF
      END IF
    END IF
  END SUBROUTINE SimulStep



  DOUBLE PRECISION FUNCTION Potential(J, Polymer, NewPos)
    INTEGER :: J, K
    DOUBLE PRECISION :: Polymer(N, 2), NewPos(2), RelPos(2), Dist

    Potential = 0.D0
    DO K=1, J-1
      RelPos = Polymer(K,:)-NewPos(:)
      Dist = 1.d0/SUM(RelPos*RelPos)
      Potential = Potential + 4*Epsilon*(-sig6*Dist**3+Sig12*Dist**6)
    END DO
  END FUNCTION Potential
END PROGRAM PERM