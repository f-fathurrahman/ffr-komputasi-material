! Recursive implementation of the Rosenbluth algorithm, according 
! to section 10.6.1 of `Computational Physics' by J.M. Thijssen, Cambridge 2007 
! The program produces a table containg the size (number of beads) N, 
! the end-to-end length squared, and a statistical error in this length.
! The results should match the form R^2 ~ (N-1)^1.5
! In figure 10.3, the straight line has the form (N-1)^1.5, and not with 
! exponent 0.75 as given in the caption.
! Program written by Jos Thijssen, 2006-2007

PROGRAM RecurRosenbluth
  IMPLICIT NONE
  INTEGER, PARAMETER :: N=250, StepNum = 10000, ThetaNum = 9, MaxConf=2*StepNum
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
    CALL InitRand(8101231)
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
  INTEGER,PARAMETER :: BlockNum = 10 ! Nr of blocks usded for data-blocking
  DOUBLE PRECISION :: Average(BlockNum),Weight, AvA, AvSq, ERrA, Alpha
  INTEGER :: Length, Num, I, Lo, Hi, BlockSize, LNum, LC
  Alpha = LOG(DBLE(N))/50.D0
  LNum = LOG(DBLE(N))/Alpha
  DO LC = 3, LNum
    Length = EXP(Alpha*LC) ! Lengths are scaled such as to obtain a more or 
!less equidistant set of points on the x-axis
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
! Recursive routine which implements the rosenbluth algorithm
    DOUBLE PRECISION, INTENT(INOUT) :: Polymer(N, 2), Weight
    INTEGER :: Length,K
    DOUBLE PRECISION :: ProbArr(ThetaNum), Theta0, Theta, &
                        ProbSum, R, TProb, NewPos(2), RealRand, &
                        RelPos(2), Dist

    addcount = addcount +1
    Theta0 = RealRand()
    ProbSum = 0.D0
    DO K=1, ThetaNum
      Theta = Theta0 + 2*PI*(K-1)/DBLE(ThetaNum)
      NewPos = Polymer(Length,:) + (/COS(Theta), SIN(Theta)/)
! Position of the new bead
      ProbArr(K) = EXP(-Potential(Length,Polymer,NewPos))
! Boltzmann weight
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
    Weight = Weight*ProbSum/(0.6D0*ThetaNum)
    
    StC(Length+1) = StC(Length+1)+1
!   Do Statistics    
    RelPos = Polymer(1,:) - Polymer(Length+1,:)
    Dist = SUM(RelPos*RelPos)
    End2End(Length+1, StC(Length+1)) = Dist
    WeightArr(Length+1, StC(Length+1)) = Weight
    IF (Length<N-1) THEN
      CALL SimulStep(Polymer, Length+1, Weight) ! recursive call
    END IF
  END SUBROUTINE SimulStep
  
    
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
END PROGRAM RecurRosenbluth
