real(8) FUNCTION get_E_ovrl()
  use globals
  IMPLICIT NONE
  ! functions
  integer :: get_index_pp

   real(8) :: RPos, Xi1, Xi2, AvXi, DPos(3), G2
   REAL :: erfc
   INTEGER :: AtomNum1, AtomNum2, Zion1, Zion2, I, J, K, N1, N2, Pos(3)
   INTEGER :: ind1, ind2
   Get_E_ovrl = 0.D0
   DO N1 = 1, N_ion
     AtomNum1 = Ions(N1)%AtomNum
     ind1 = get_index_pp(AtomNum1)
     Zion1 = PP_Params(ind1)%Zion
     Xi1 = PP_Params(ind1)%Xi
     DO N2 = N1, N_ion
       AtomNum2 = Ions(N2)%AtomNum
       ind2 = get_index_pp(AtomNum2)
       Zion2 = PP_Params(ind2)%Zion
       Xi2 = PP_Params(ind2)%Xi
       DPos = Ions(N1)%R_I(:)-Ions(N2)%R_I(:)
       AvXi = SQRT(2.D0*(Xi1**2+Xi2**2))
       DO I=-2, 2
         DO J=-2, 2
           DO K=-2, 2
             IF ((N1 .NE. N2) .OR. (I.NE.0) .OR. (J.NE.0) .OR. (K.NE.0)) THEN
               Pos = (/I*BoxL,J*BoxL,K*BoxL/)
               RPos = SQRT(DOT_PRODUCT(DPos-Pos, DPos-Pos))
               Get_E_ovrl = Get_E_ovrl+Zion1*Zion2/Rpos*DBLE(erfc(REAL(Rpos/AvXi)))
             END IF
           END DO !K
         END DO !J
       END DO !I
     END DO !N2
   END DO !N1
END FUNCTION