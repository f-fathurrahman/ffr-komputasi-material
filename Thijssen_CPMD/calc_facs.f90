SUBROUTINE calc_facs(I,J,K,N,G2,StructFac)
  use globals
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: I, J, K, N
  INTEGER, INTENT(OUT) :: G2
  complex(8), INTENT(OUT) :: StructFac
  INTEGER :: II, JJ, KK
  INTEGER :: pos(3)

  II = I-INT((2.D0*I)/GridSize)*GridSize
  JJ = J-INT((2.D0*J)/GridSize)*GridSize
  KK = K-INT((2.D0*K)/GridSize)*GridSize

  G2 = II*II+JJ*JJ+KK*KK
  pos = (/II, JJ, KK/)
  StructFac = EXP(-Im*2*PI*DOT_PRODUCT(Ions(N)%R_I,pos)/BoxL)
END SUBROUTINE


