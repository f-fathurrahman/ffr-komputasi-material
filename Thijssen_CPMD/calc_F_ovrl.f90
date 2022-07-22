SUBROUTINE calc_F_ovrl(OvrlForce)
  use globals
  ! Overlay part of the ionic force is calculated  
    IMPLICIT NONE
  
    complex(8), INTENT(OUT) :: OvrlForce(N_ion,3)
    
    real(8)    :: RPos(3), AbsRPos, Xi1, Xi2, AvXi, DPos(3), &
                                prefac1, prefac2, prefactot
    REAL                :: erfc
    
    INTEGER :: AtomNum1, AtomNum2, Zion1, Zion2, Pos(3), I, J, K, M, N1, N2
    INTEGER :: ind1, ind2
    ! Function
    integer :: get_index_pp

    OvrlForce = CMPLX(0.D0)

    DO N1 = 1, N_ion
      AtomNum1 = Ions(N1)%AtomNum
      ind1 = get_index_pp(AtomNum1)
      Zion1 = PP_Params(ind1)%Zion
      Xi1 = PP_Params(ind1)%Xi
      DO N2 = N1, N_ion
        AtomNum2 = Ions(N2)%AtomNum
        ind2 = get_index_PP(AtomNum2)
        Zion2 = PP_Params(ind2)%Zion
        Xi2 = PP_Params(ind2)%Xi
        
        DPos = Ions(N1)%R_I - Ions(N2)%R_I
        AvXi = SQRT(2.D0*(Xi1**2+Xi2**2))
        prefac1 = Zion1*Zion2
        prefac2 = prefac1*2/(AvXi*sqrt(PI))
        DO I=-2, 2
          DO J=-2, 2
            DO K=-2, 2
              IF ((N1 .NE. N2) .OR. (I.NE.0) .OR. (J.NE.0) .OR. (K.NE.0)) THEN
                Pos = (/I, J, K/)     
                AbsRPos = 0.D0
                RPos = DPos-Pos*BoxL
                AbsRPos = SQRT(DOT_PRODUCT(RPos,RPos))
                prefactot = prefac1/AbsRPos**3*DBLE(erfc(Real(AbsRpos/AvXi))) + &
                            prefac2/AbsRPos**2*exp(-(AbsRPos/AvXi)**2)
                OvrlForce(N1,:)= OvrlForce(N1,:) + prefactot * RPos 
                OvrlForce(N2,:)= OvrlForce(N2,:) - prefactot * RPos 
              END IF
            END DO
          END DO
        END DO
      END DO
    END DO
END SUBROUTINE