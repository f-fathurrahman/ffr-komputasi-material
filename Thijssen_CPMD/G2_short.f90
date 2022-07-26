FUNCTION G2_Short(I,J,K) RESULT(res)
  !
  USE globals, ONLY: GridSize
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: I, J, K
  INTEGER :: II, JJ, KK
  REAL(8) :: res

  II = I - INT((2.D0*I)/GridSize)*GridSize
  JJ = J - INT((2.D0*J)/GridSize)*GridSize
  KK = K - INT((2.D0*K)/GridSize)*GridSize

  res = II*II + JJ*JJ + KK*KK
END FUNCTION

