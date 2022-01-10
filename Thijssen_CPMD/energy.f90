MODULE energy
! In this module, the Kohn-Sham matrix is set up and diagonalised.
! The density is then fed into a new KS Hamiltonian etcetera
USE pseudo
USE excorr
USE Globals
USE utilities


CONTAINS


  SUBROUTINE Total_Energy(NZCoeffs, E_KS)
! Evaluate the energy for a given solution Coeffs_K

  IMPLICIT NONE
  DOUBLE COMPLEX, INTENT(IN) :: NZCoeffs(N_orbitals, NoOfPW)
  INTEGER :: I1, J1, K1, I2, J2, K2, IIndex, JIndex, G2, N, IT, JT, KT, &
             IElec, M, PP_Index, AtomNum, IOrb
  DOUBLE COMPLEX, ALLOCATABLE :: TempVec(:,:,:), TempVec2_R(:,:,:), &
                                 TempVec2_K(:,:,:)
  DOUBLE COMPLEX :: II, PreFac
  DOUBLE COMPLEX, INTENT(OUT) :: E_KS
  DOUBLE COMPLEX :: E_kin, E_locPP, E_locPPsr, E_hartree, E_xc, E_nonlocPP, &
                    E_ovrl, E_self, E_totSelf, E_core

  DOUBLE PRECISION :: RPos, h_1s, h_2s, h_1p

  E_KS=CMPLX(0.D0)
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

  print *, 'nr el. ', N_electron
  PreFac = 2*PI*PI/BoxL**2
  DO Iorb = 1, N_orbitals
    E_kin = E_kin + FillFac(Iorb)*PreFac*&
        SUM(NZCoeffs(Iorb,:)*CONJG(NZCoeffs(Iorb,:))*G2Grid)
  END DO
  print '(A23 F15.8)', 'kin:', DBLE(E_kin)
   
! SHORT RANGE PART OF LOCAL PP
  
  E_locPPsr = Omega*SUM(totShortLocal*CONJG(Density_K))
  print '(A23 F15.8)', 'pp_sr:', DBLE(E_locPPsr)

  E_locPP = Omega*SUM(PseudoGrid*CONJG(Density_K))
  print '(A23 F15.8)', 'loc pp:', DBLE(E_locPP)

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
  print '(A23 F15.8)', 'xc:', DBLE(E_xc)

! HARTREE
  PreFac = BoxL**2*Omega/(2*PI)
  E_hartree = PreFac*SUM(Gmin2Grid*Density_K*CONJG(Density_K))
  print '(A23 F15.8)', 'Hartree energy:', DBLE(E_hartree)

! Nonlocal PsP
  E_nonlocPP = CMPLX(0.D0)
  DO IOrb = 1, N_orbitals
    DO N = 1, N_ion
      AtomNum = Ions(N)%AtomNum
      IF (AtomNum > 4) THEN
        PP_Index = GetIndexPP(AtomNum)
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
  print '(A23 F15.8)', 'Nonlocal Psp:', DBLE(E_nonlocPP)

! CORE ENERGY
  PreFac = BoxL**2*Omega/(2*PI)
  E_core = PreFac*SUM(Gmin2Grid*totCoreCharge*CONJG(totCoreCharge))

  print '(A23 F15.8)', 'Local core energy:', DBLE(E_core)

  E_self = PreFac*SUM(Gmin2Grid*(totCoreCharge+Density_K)*CONJG(totCoreCharge+Density_K))

  print '(A23 F15.8)', 'self-energy:', DBLE(E_self)

  E_ovrl = Get_E_ovrl()
  print '(A23 F15.8)', 'ovrl:', DBLE(E_ovrl)
  DO N=1, N_ion
    E_totSelf=E_totSelf+Get_E_Self(Ions(N)%AtomNum)
  END DO
  E_KS = E_kin + E_locPPsr + E_xc + E_nonlocPP + E_self + E_ovrl - E_totSelf
  print '(A23 F15.8)', 'Total energy:', DBLE(E_KS) 
  print *
  END SUBROUTINE Total_Energy



  SUBROUTINE Check_ConstEnergy(NZCoeffsDot, R_ionDot)
  IMPLICIT NONE
  DOUBLE COMPLEX, INTENT(IN) :: NZCoeffsDot(N_orbitals, NoOfPW), &
                                R_ionDot(N_ion, 3)
  DOUBLE COMPLEX :: E
  INTEGER        :: N

  E = CMPLX(0.D0)
  DO N=1, N_ion
    E = E + SUM(R_ionDot(N,:)*CONJG(R_ionDot(N,:)))/(2*Ions(N)%Mass) 
  END DO
  E = E + SUM(NZCoeffsDot*CONJG(NZCoeffsDot))/(2*mu)
  print '(A23 D15.3 )', 'Constant Energy check:', DBLE(E)
  print *
  END SUBROUTINE Check_ConstEnergy



  
  DOUBLE PRECISION FUNCTION Get_E_ovrl()
   IMPLICIT NONE

   DOUBLE PRECISION :: RPos, Xi1, Xi2, AvXi, DPos(3), G2
   REAL :: erfc
   INTEGER :: AtomNum1, AtomNum2, Zion1, Zion2, I, J, K, N1, N2, Pos(3)
   INTEGER :: ind1, ind2
   Get_E_ovrl = 0.D0
   DO N1 = 1, N_ion
     AtomNum1 = Ions(N1)%AtomNum
     ind1 = GetIndexPP(AtomNum1)
     Zion1 = PP_Params(ind1)%Zion
     Xi1 = PP_Params(ind1)%Xi
     DO N2 = N1, N_ion
       AtomNum2 = Ions(N2)%AtomNum
       ind2 = GetIndexPP(AtomNum2)
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
  END FUNCTION Get_E_ovrl


  
  DOUBLE PRECISION FUNCTION Get_E_self(AtomNum)
  ! Returns the long-and short range part of the FT of the local pseudopotential
  ! This is useful for calculating the Kohn-Sham Hamiltonian, but not for the 
  ! total energy, as the long range part should be treated differently in that case
  ! (Ewald sum)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: AtomNum
  DOUBLE PRECISION :: Xi, Zion
  INTEGER :: ind

  ind = GetIndexPP(AtomNum)
  Zion = DBLE(PP_Params(ind)%Zion)
  Xi = PP_Params(ind)%Xi

  Get_E_Self = Zion*Zion/(2*SQRT(PI)*Xi)

END FUNCTION Get_E_Self


  
END MODULE energy

