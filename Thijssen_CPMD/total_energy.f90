SUBROUTINE Total_Energy(NZCoeffs, E_KS)
  use globals
  use energy_components
  ! Evaluate the energy for a given solution Coeffs_K
  IMPLICIT NONE
  complex(8), INTENT(IN) :: NZCoeffs(N_orbitals, NoOfPW)
  INTEGER :: I1, J1, K1, I2, J2, K2, IIndex, JIndex, G2, N, IT, JT, KT, &
             IElec, M, PP_Index, AtomNum, IOrb
  complex(8), ALLOCATABLE :: TempVec(:,:,:), TempVec2_R(:,:,:), &
                                 TempVec2_K(:,:,:)
  complex(8) :: II, PreFac
  complex(8), INTENT(OUT) :: E_KS

  real(8) :: RPos, h_1s, h_2s, h_1p

  ! Function
  real(8) :: Get_E_Self, Get_E_ovrl
  INTEGER :: get_index_pp
  real(8) :: epsilon_xc

  E_KS = CMPLX(0.D0)
  E_kin = CMPLX(0.D0)
  E_locPP = CMPLX(0.D0)
  E_nonlocPP = CMPLX(0.D0)
  E_locPPsr = CMPLX(0.D0)
  E_xc = CMPLX(0.D0)
  E_ovrl = CMPLX(0.D0)
  E_self = CMPLX(0.D0)
  E_totSelf = CMPLX(0.D0)
  
  II = CONJG(Im)

  ! KINETIC ENERGY
  
  PreFac = 2*PI*PI/BoxL**2
  DO Iorb = 1, N_orbitals
    E_kin = E_kin + FillFac(Iorb)*PreFac*&
        SUM(NZCoeffs(Iorb,:)*CONJG(NZCoeffs(Iorb,:))*G2Grid)
  END DO
     
  !
  ! SHORT RANGE PART OF LOCAL PP
  !
  E_locPPsr = Omega*SUM(totShortLocal*CONJG(Density_K))

  E_locPP = Omega*SUM(PseudoGrid*CONJG(Density_K))

  ! EXCHANGE CORRELATION ENERGY
  DO I1 = 0, GridSize-1
    DO J1 = 0, GridSize-1
      DO K1 = 0, GridSize-1
        ! Evaluate this expression in real space!!!!!!!!!!!!
        E_xc = E_xc+CONJG(Density_R(I1, J1, K1))*&
               CMPLX(epsilon_xc(I1,J1,K1))*Omega/GridSize**3
      END DO
    END DO
  END DO

  ! HARTREE
  PreFac = BoxL**2*Omega/(2*PI)
  E_hartree = PreFac*SUM(Gmin2Grid*Density_K*CONJG(Density_K))

  ! Nonlocal PsP
  E_nonlocPP = CMPLX(0.D0)
  DO IOrb = 1, N_orbitals
    DO N = 1, N_ion
      AtomNum = Ions(N)%AtomNum
      IF (AtomNum > 4) THEN
        PP_Index = get_index_pp(AtomNum)
        h_1s = PP_Params(PP_Index)%h_1s
        h_2s = PP_Params(PP_Index)%h_2s
        h_1p = PP_Params(PP_Index)%h_1p
        PreFac = SUM(NonLocal(:,N,1)*CONJG(NZCoeffs(IOrb,:)))
        E_nonlocPP  = E_nonlocPP + FillFac(IOrb)*h_1s*PreFac*CONJG(PreFac)
        IF (PP_Params(PP_Index)%MaxL>0) THEN
          PreFac = SUM(NonLocal(:,N,2)*CONJG(NZCoeffs(IOrb,:)))
          E_nonlocPP  = E_nonlocPP + FillFac(IOrb)*h_2s*PreFac*CONJG(PreFac)
          PreFac = SUM(NonLocal(:,N,3)*CONJG(NZCoeffs(IOrb,:)))
          E_nonlocPP  = E_nonlocPP + FillFac(IOrb)*h_1p*PreFac*CONJG(PreFac)
          PreFac = SUM(NonLocal(:,N,4)*CONJG(NZCoeffs(IOrb,:)))
          E_nonlocPP  = E_nonlocPP + FillFac(IOrb)*h_1p*PreFac*CONJG(PreFac)
          PreFac = SUM(NonLocal(:,N,5)*CONJG(NZCoeffs(IOrb,:)))
          E_nonlocPP  = E_nonlocPP + FillFac(IOrb)*h_1p*PreFac*CONJG(PreFac)
        END IF
      END IF
    END DO  
  END DO

  ! CORE ENERGY
  PreFac = BoxL**2*Omega/(2*PI)
  E_core = PreFac*SUM(Gmin2Grid*totCoreCharge*CONJG(totCoreCharge))

  E_self = PreFac*SUM(Gmin2Grid*(totCoreCharge+Density_K)*CONJG(totCoreCharge+Density_K))

  E_ovrl = Get_E_ovrl()

  DO N=1, N_ion
    E_totSelf = E_totSelf + Get_E_Self(Ions(N)%AtomNum)
  END DO
  E_KS = E_kin + E_locPPsr + E_xc + E_nonlocPP + E_self + E_ovrl - E_totSelf

END SUBROUTINE