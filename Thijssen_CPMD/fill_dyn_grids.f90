SUBROUTINE fill_dyn_grids()
  use globals
  ! In this subroutine the grids are filled that are dependent on the position
  ! of the Ions
  IMPLICIT NONE
  
  INTEGER :: I, J, K, N, G2, Ind, I1, J1, K1, PP_Index, AtomNum
  complex(8) :: StructFac
  ! Function
  complex(8) :: nonloc_pp
  real(8) :: core_dens, short_pp, local_pp, G2_Short
  INTEGER :: get_index_pp

  totShortLocal = CMPLX(0.D0)
  PseudoGrid = CMPLX(0.D0)
  ShortLocal = CMPLX(0.D0)
  totCoreCharge = CMPLX(0.D0)
  NonLocal = CMPLX(0.D0)

  Ind = 0
  DO I=0, GridSize-1
    DO J=0, GridSize-1
      DO K=0, GridSize-1
        G2 = G2_Short(I,J,K)
        IF (G2<Gmax**2) ind = ind + 1
          DO N = 1, N_ion
            AtomNum = Ions(N)%AtomNum
            PP_Index = get_index_pp(AtomNum)
            I1 = I-INT((2.D0*I)/GridSize)*GridSize
            J1 = J-INT((2.D0*J)/GridSize)*GridSize
            K1 = K-INT((2.D0*K)/GridSize)*GridSize
            CALL calc_facs(I, J, K, N, G2, StructFac)
            PseudoGrid(I,J,K)=PseudoGrid(I,J,K)+ & 
                Local_PP(G2,Ions(N)%AtomNum)*StructFac
            ShortLocal(N,I,J,K) = ShortLocal(N,I,J,K) +&
                              Short_PP(G2,Ions(N)%AtomNum)*StructFac
            totShortLocal(I,J,K) = totShortLocal(I,J,K) + &
                              ShortLocal(N,I,J,K)
            CoreCharge(N,I,J,K) = Core_Dens(G2,Ions(N)%AtomNum)*StructFac
            totCoreCharge(I,J,K) = totCoreCharge(I,J,K) + CoreCharge(N,I,J,K)
            IF (G2<Gmax**2) THEN
              NonLocal(Ind,N,1) = NonLoc_pp(I1,J1,K1,AtomNum,0,0,1)*StructFac
              IF (PP_Params(PP_Index)%MaxL>0) THEN
                NonLocal(Ind,N,2) = NonLoc_pp(I1,J1,K1,AtomNum,0,0,2)*StructFac
                NonLocal(Ind,N,3) = NonLoc_pp(I1,J1,K1,AtomNum,1,1,1)*StructFac
                NonLocal(Ind,N,4) = NonLoc_pp(I1,J1,K1,AtomNum,1,0,1)*StructFac
                NonLocal(Ind,N,5) = NonLoc_pp(I1,J1,K1,AtomNum,1,-1,1)*StructFac
              END IF
            END IF !G2<Gmax**2
        END DO !N
      END DO !K
    END DO !J
  END DO !I
END SUBROUTINE
