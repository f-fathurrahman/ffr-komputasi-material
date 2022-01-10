MODULE force
! In this module, the Kohn-Sham matrix is set up and diagonalised.
! The density is then fed into a new KS Hamiltonian etcetera
USE pseudo
USE grids
USE excorr
USE energy
USE Globals
USE utilities


CONTAINS

  SUBROUTINE Calc_OrbForce(NZCoeffs, OrbForce)
! Force vector of the wavegradient is calculated

  IMPLICIT NONE
  DOUBLE COMPLEX, INTENT(OUT) :: OrbForce(N_orbitals, NoOfPW)
  DOUBLE COMPLEX, INTENT(IN) :: NZCoeffs(N_orbitals, NoOfPW)
                              
  INTEGER :: I1, J1, K1, I2, J2, K2, IIndex, JIndex, N, IT, JT, KT, &
             G2, PP_Index, AtomNum, L, Iorb
  DOUBLE COMPLEX, ALLOCATABLE :: TempVec_K(:,:,:), TempVec_R(:,:,:), &
                                 TempForce_K(:,:,:,:), TempForce_R(:,:,:,:), &
                                 Coeffs_R(:,:,:,:)
  DOUBLE COMPLEX :: PreFac, II, F
  DOUBLE PRECISION :: h_1s, h_1p, h_2s, hfac
! ALLOCATE STORAGE
  ALLOCATE (Coeffs_R(N_orbitals, 0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
  ALLOCATE (TempVec_K(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
  ALLOCATE (TempVec_R(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
  ALLOCATE (TempForce_K(N_orbitals, 0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
  ALLOCATE (TempForce_R(N_orbitals, 0:GridSize-1, 0:GridSize-1, 0:GridSize-1))

  CALL CalcDensAndCoeffs_R(NZCoeffs, Coeffs_R)

! KINETIC TERM
  PreFac = CMPLX(2*PI*PI/BoxL**2)

  DO Iorb = 1, N_orbitals
    OrbForce(Iorb,:) = -PreFac*G2Grid*NZCoeffs(Iorb,:)
  END DO
! EXCHANGE CORRELATION 

  DO I1 = 0, GridSize-1
    DO J1 = 0, GridSize-1
      DO K1 = 0, GridSize-1
        TempVec_R(I1, J1, K1) = CMPLX(Vxc(I1, J1, K1))
      END DO
    END DO
  END DO
  CALL BackWard_FFT(GridSize, TempVec_R, TempVec_K)

! LOCAL PSEUDOPOTENTIAL
  TempVec_K = TempVec_K + PseudoGrid

! HARTREE
  PreFac = BoxL**2/PI
  TempVec_K = TempVec_K + PreFac*Density_K*Gmin2Grid
  
! Calculate contrbution from local potential
  CALL Forward_FFT(GridSize, TempVec_K, TempVec_R)
  DO N=1, N_orbitals
    TempForce_R(N,:,:,:) = TempVec_R*Coeffs_R(N,:,:,:)*SQRT(Omega)
    CALL Backward_FFT(GridSize, TempForce_R(N,:,:,:), TempForce_K(N,:,:,:))
  END DO
  DO IIndex = 1, NoOfPW
    I1 = GridPos(IIndex,1)
    J1 = GridPos(IIndex,2)
    K1 = GridPos(IIndex,3)
    DO N=1, N_orbitals
      OrbForce(N, IIndex) = OrbForce(N, IIndex) - TempForce_K(N, I1, J1, K1)
    END DO
  END DO

! NONLOCAL PSEUDOPOTENTIAL
  DO Iorb = 1, N_orbitals
    DO N=1, N_ion
      AtomNum = Ions(N)%AtomNum
      IF (AtomNum > 4) THEN
        PP_Index = GetIndexPP(AtomNum)
        h_1s = PP_Params(PP_Index)%h_1s
        IF (PP_Params(PP_Index)%MaxL>0) THEN
          h_2s = PP_Params(PP_Index)%h_2s
          h_1p = PP_Params(PP_Index)%h_1p
        END IF
        F = 0.D0
        DO IIndex=1, NoOfPW
          F = F+NonLocal(IIndex,N,1)*CONJG(NZCoeffs(Iorb,IIndex))
        END DO
        DO IIndex=1, NoOfPW
          OrbForce (Iorb, IIndex) = OrbForce (Iorb, IIndex) - &
             CONJG(F)*h_1s*NonLocal(IIndex,N, 1)
        END DO
        IF (PP_Params(PP_Index)%MaxL>0) THEN
          DO L=2, 5
            IF (L==2) THEN
              hfac = h_2s
            ELSE
              hfac = h_1p
              END IF
            F = 0.D0
            DO IIndex=1, NoOfPW
              F=F+NonLocal(IIndex,N,L)*CONJG(NZCoeffs(Iorb,IIndex))
            END DO
            DO IIndex=1, NoOfPW
              OrbForce (Iorb, IIndex) = OrbForce (Iorb, IIndex) - &
                 CONJG(F)*hfac*NonLocal(IIndex,N,L)
            END DO
          END DO
        END IF
      END IF
    END DO
  END DO
  
  DO N=1, N_orbitals
    OrbForce(N,:) = OrbForce(N,:)*FillFac(N)
  END DO
  END SUBROUTINE Calc_OrbForce


  
  SUBROUTINE Calc_IonForce(NZCoeffs, IonForce)
! The force on the Ions are calculated 
    IMPLICIT NONE
    DOUBLE COMPLEX, INTENT(IN)  :: NZCoeffs(N_orbitals, NoOfPW)
    DOUBLE COMPLEX, INTENT(OUT) :: IonForce(N_ion,3)
    DOUBLE COMPLEX, ALLOCATABLE :: ESForce(:,:), LocalForce(:,:), &
                                   NonLocalForce(:,:)
    INTEGER                     :: I
! Calculation of the ElectroStatic part
    ALLOCATE(ESForce(N_ion,3))
    CALL Calc_F_ES(ESForce)
 
! Calculation of the PP local part
    ALLOCATE(LocalForce(N_ion,3))
    CALL Calc_F_Local(LocalForce)
    
! Calculation of the PP nonlocal part
    ALLOCATE(NonLocalForce(N_ion,3))
    CALL Calc_F_nonlocal(NZCoeffs,NonLocalForce)
 
    IonForce = EsForce + LocalForce + NonLocalForce
    IF (PrintOut) THEN
      DO I=1, N_Ion 
         print '(A23 I3, F15.8, F15.8, F15.8)', 'ES Force ion:',I, DBLE(EsForce(I,:)) 
      END DO
      DO I=1, N_Ion 
         print '(A23 I3, F15.8, F15.8, F15.8)', 'Local Force ion:',I, DBLE(LocalForce(I,:)) 
      END DO
      DO I=1, N_Ion 
         print '(A23 I3, F15.8, F15.8, F15.8)', 'NonLocal  Force ion:',I, DBLE(NonLocalForce(I,:)) 
      END DO
      DO I=1, N_Ion 
         print '(A23 I3, F15.8, F15.8, F15.8)', 'Total Ionic Force ion:',I, DBLE(IonForce(I,:)) 
      END DO
      print *
    END IF
  END SUBROUTINE Calc_IonForce

  
   

  
  SUBROUTINE Calc_F_ES(ESForce)
    IMPLICIT NONE

    DOUBLE COMPLEX, INTENT(OUT) :: ESForce(N_ion,3)
    DOUBLE COMPLEX, ALLOCATABLE :: OvrlForce(:,:)
    DOUBLE COMPLEX              :: PreFac
    INTEGER                     :: N, M
    
    ALLOCATE(OvrlForce( N_ion, 3))
    
    ESForce = CMPLX(0.D0) 
    
    PreFac = 2*Im*BoxL**4
    DO N = 1, N_ion
      DO M = 1, 3
        ESForce(N,M) = ESForce(N,M) + PreFac*SUM(Gmin2Grid*GGrid(M,:,:,:)*CoreCharge(N,:,:,:)*CONJG(Density_K+totCoreCharge))
      END DO
    END DO
    CALL Calc_F_ovrl( OvrlForce )
    ESForce = ESForce + OvrlForce
    
  END SUBROUTINE Calc_F_ES

  
  SUBROUTINE  Calc_F_ovrl(OvrlForce)
! Overlay part of the ionic force is calculated  
    IMPLICIT NONE
  
    DOUBLE COMPLEX, INTENT(OUT) :: OvrlForce(N_ion,3)
    
    DOUBLE PRECISION    :: RPos(3), AbsRPos, Xi1, Xi2, AvXi, DPos(3), &
                                prefac1, prefac2, prefactot
    REAL                :: erfc
    
    INTEGER :: AtomNum1, AtomNum2, Zion1, Zion2, Pos(3), I, J, K, M, N1, N2
    INTEGER :: ind1, ind2
    
    OvrlForce = CMPLX(0.D0)

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
  END SUBROUTINE Calc_F_ovrl



  
  SUBROUTINE Calc_F_local(LocalForce)
     IMPLICIT NONE
     DOUBLE COMPLEX              :: PreFac
     INTEGER                     :: N, M
     
     DOUBLE COMPLEX, INTENT(OUT) :: LocalForce(N_ion,3)
     PreFac = 2*PI*Im*BoxL**2
     LocalForce = CMPLX(0.D0)
     DO N=1, N_ion
       DO M=1, 3
            LocalForce(N,M) = LocalForce(N,M) + &
                        preFac*SUM(GGrid(M,:,:,:)*ShortLocal(N,:,:,:)*CONJG(Density_K))
       END DO
     END DO
  END SUBROUTINE



  
  SUBROUTINE Calc_F_nonlocal(NZCoeffs,NonLocalForce)
     IMPLICIT NONE
     DOUBLE COMPLEX, INTENT(OUT) :: NonLocalForce(N_ion,3)
     DOUBLE COMPLEX, INTENT(IN)  :: NZCoeffs(N_orbitals, NoOfPW)
     DOUBLE COMPLEX              :: F(N_ion), dFdRI(N_ion,3)
     DOUBLE PRECISION            :: h(5), prefac
     INTEGER                     :: N, M, K, Iorb, MaxM, AtomNum, &
                                    PP_Index, GCnt, G

     prefac = 2*PI/BoxL
     NonLocalForce = CMPLX(0.D0)
     DO N = 1, N_ion
       DO Iorb = 1, N_orbitals  
         AtomNum = Ions(N)%AtomNum
         IF (AtomNum > 4) THEN
           PP_Index = GetIndexPP(AtomNum)
           h = (/PP_Params(PP_Index)%h_1s, & 
                 PP_Params(PP_Index)%h_2s, &
                 PP_Params(PP_Index)%h_1p, &
                 PP_Params(PP_Index)%h_1p, &
                 PP_Params(PP_Index)%h_1p/)
           IF (PP_Params(PP_Index)%MaxL>0) THEN
               MaxM = 5
           ELSE
               MaxM = 1
           END IF
           DO M = 1, MaxM
              F(N) = SUM( NonLocal(:,N,M) * CONJG(NZCoeffs(IOrb,:)) )
              DO K = 1, 3
                dFdRI(N,K) = CMPLX(0.D0)
                DO GCnt = 1, NoOfPW
                  G = GridPos(GCnt,K)
                  G = G-INT((2.D0*G)/GridSize)*GridSize
                  dFdRI(N,K) = dFdRI(N,K)-Im*prefac*G*NonLocal(GCnt,N,M)*CONJG(NZCoeffs(IOrb,GCnt))
                END DO
              END DO !K
              NonLocalForce(N,:) = NonLocalForce(N,:) - &
                     (CONJG(dFdRI(N,:))*h(M)*F(N)+dFdRI(N,:)*h(M)*CONJG(F(N)))* FillFac(Iorb)
           END DO !M
         END IF !AtomNum>4
       END DO !Iorb
     END DO !N
  END SUBROUTINE Calc_F_nonlocal
  
END MODULE force

