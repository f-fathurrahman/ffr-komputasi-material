real(8) FUNCTION G2_Short(I,J,K)
  use globals
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: I, J, K
  INTEGER :: II, JJ, KK

  II = I-INT((2.D0*I)/GridSize)*GridSize
  JJ = J-INT((2.D0*J)/GridSize)*GridSize
  KK = K-INT((2.D0*K)/GridSize)*GridSize

  G2_Short = II*II+JJ*JJ+KK*KK
END FUNCTION G2_Short