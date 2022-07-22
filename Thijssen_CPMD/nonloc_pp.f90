!--------------------------------------------------------
complex(8) FUNCTION nonloc_pp(I, J, K, AtomNum, L, M, II)
!--------------------------------------------------------
! Returns the term P_\alpha^I from Eq. (178)
  ! Parameters: I, J, K: G-vector
  !             L, M : Angular Momentum
  !             II : Index of projector (for s-waves)
  use globals
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: I, J, K, L, M, II, AtomNum;
  INTEGER :: G2
  real(8) :: Pl, Xi, RInv
  complex(8) :: Ylm
  INTEGER :: ind
  ! Functions
  real(8) :: G2_Short
  INTEGER :: get_index_pp

  ind = get_index_pp(AtomNum)
  Ylm = CMPLX(1/SQRT(4*PI))

  G2 = G2_Short(I,J,K)
  IF (L==0) THEN
    Xi = PP_Params(ind)%r_s/BoxL
    IF (II==1) THEN
      Pl = Xi**1.5D0*PI**1.25D0*4*SQRT(2.D0)*EXP(-2*PI*PI*G2*Xi*Xi)
    ELSE
      Pl = Xi**1.5D0*PI**1.25D0*8/SQRT(7.5D0)*EXP(-2*PI*PI*G2*Xi*Xi)*&
                (3-4*PI*PI*G2*Xi*Xi)
    END IF
  ELSE
    Xi = PP_Params(ind)%r_p/BoxL
    Pl = Xi**2.5D0*PI**1.25D0*8*2*PI*SQRT(G2/3.D0)*EXP(-2*PI*PI*G2*Xi*Xi)
    IF (G2 /=0) THEN
      RInv = 1/SQRT(DBLE(G2))
      SELECT CASE (M)
        CASE( 1); Ylm = SQRT(0.5D0)*Ylm*CMPLX(-DBLE(I), -DBLE(J))
        CASE( 0); Ylm = Ylm*CMPLX(DBLE(K), 0)
        CASE(-1); Ylm = SQRT(0.5D0)*Ylm*CMPLX(DBLE(I), -DBLE(J))
      END SELECT
      Ylm = SQRT(3.D0)*Ylm*RInv
    ELSE
      Ylm = CMPLX(0.D0)
    END IF   
  END IF
  nonloc_pp = Pl*Ylm
END FUNCTION

