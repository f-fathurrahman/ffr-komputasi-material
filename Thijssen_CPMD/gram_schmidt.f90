SUBROUTINE Gram_Schmidt(Vectors, Number, Dimen)
  use globals
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Number, Dimen
  complex(8), INTENT(INOUT) :: Vectors(Number, Dimen)
  INTEGER Iorb1, Iorb2
  complex(8) :: IP

  DO Iorb1 = 1, Number
    DO Iorb2 = 1, Iorb1-1
      IP = DOT_PRODUCT(Vectors(Iorb2,:),Vectors(Iorb1,:))
      Vectors(Iorb1,:) = Vectors(IOrb1,:)-IP*Vectors(Iorb2,:)
    END DO
    IP = DOT_PRODUCT(Vectors(Iorb1,:),Vectors(Iorb1,:))
    IP = 1/SQRT(IP)
    Vectors(Iorb1,:) = Vectors(Iorb1,:)*IP
  END DO

END SUBROUTINE Gram_Schmidt