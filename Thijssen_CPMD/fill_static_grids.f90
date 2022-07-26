!-----------------------------
SUBROUTINE fill_static_grids()
!-----------------------------
  USE globals, ONLY: GridSize, Gmin2Grid, G2Grid, GridPos, GridIndex, &
                     Gmax, GGrid
  IMPLICIT NONE 
  ! In this subroutine the grids are defined that are independent on the
  ! position of the Ions
  ! G2 stands for G squared, Gmin2 for 1/G Squared
  INTEGER :: I, J, K, Ind
  INTEGER :: II, JJ, KK
  real(8) :: G2

  ! Function
  REAL(8) :: G2_Short

  Ind = 0
    
  GGrid = 0
  GridIndex = 0
  GridPos = 0
  DO I = 0, GridSize-1
    DO J = 0, GridSize-1
      DO K = 0, GridSize-1
        G2 = G2_Short(I,J,K)
        II = I - INT( (2.D0*I)/GridSize )*GridSize
        JJ = J - INT( (2.D0*J)/GridSize )*GridSize
        KK = K - INT( (2.D0*K)/GridSize )*GridSize
        IF( G2 < 4*GMax**2 ) THEN
          GGrid(:,I,J,K) = (/ II, JJ, KK /)
        ELSE
          GGrid(:,I,J,K) = (/0, 0, 0/)
        ENDIF
        IF( G2 < Gmax**2 ) THEN
          Ind = Ind + 1
          GridIndex(I, J, K) = Ind
          GridPos(Ind,:) = (/I,J,K/)
          G2Grid(Ind) = G2
        ENDIF
        IF( G2 /= 0 ) THEN
          Gmin2Grid(I,J,K) = 1.D0/G2
        ELSE
          Gmin2Grid(I,J,K) = 0.D0
        ENDIF
      ENDDO !K
    ENDDO !J
  ENDDO !I
  
END SUBROUTINE

