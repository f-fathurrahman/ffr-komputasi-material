SUBROUTINE Gram_Schmidt(Vectors, nrow, Dimen)
  USE globals
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nrow, Dimen
  COMPLEX(8), INTENT(INOUT) :: Vectors(nrow, Dimen)
  INTEGER :: Iorb1, Iorb2
  COMPLEX(8) :: IP

  DO Iorb1 = 1, nrow
    DO Iorb2 = 1, Iorb1-1
      IP = DOT_PRODUCT(Vectors(Iorb2,:),Vectors(Iorb1,:))
      Vectors(Iorb1,:) = Vectors(IOrb1,:)-IP*Vectors(Iorb2,:)
    END DO
    IP = DOT_PRODUCT(Vectors(Iorb1,:),Vectors(Iorb1,:))
    IP = 1/SQRT(IP)
    Vectors(Iorb1,:) = Vectors(Iorb1,:)*IP
  ENDDO

END SUBROUTINE

