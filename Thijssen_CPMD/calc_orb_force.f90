SUBROUTINE calc_orb_force(NZCoeffs, OrbForce)
  ! Force vector of the wavegradient is calculated
  USE globals
  IMPLICIT NONE
  COMPLEX(8), INTENT(OUT) :: OrbForce(N_orbitals, NoOfPW)
  COMPLEX(8), INTENT(IN) :: NZCoeffs(N_orbitals, NoOfPW)
                              
  INTEGER :: I1, J1, K1, IIndex, N, &
             PP_Index, AtomNum, L, Iorb
  COMPLEX(8), ALLOCATABLE :: TempVec_K(:,:,:), TempVec_R(:,:,:), &
                             TempForce_K(:,:,:,:), TempForce_R(:,:,:,:), &
                             Coeffs_R(:,:,:,:)
  COMPLEX(8) :: PreFac, F
  REAL(8) :: h_1s, h_1p, h_2s, hfac
  ! Functions
  REAL(8) :: Vxc
  INTEGER :: Get_index_pp

  ! ALLOCATE STORAGE
  ALLOCATE( Coeffs_R(N_orbitals, 0:GridSize-1, 0:GridSize-1, 0:GridSize-1) )
  ALLOCATE( TempVec_K(0:GridSize-1, 0:GridSize-1, 0:GridSize-1) )
  ALLOCATE( TempVec_R(0:GridSize-1, 0:GridSize-1, 0:GridSize-1) )
  ALLOCATE( TempForce_K(N_orbitals, 0:GridSize-1, 0:GridSize-1, 0:GridSize-1) )
  ALLOCATE( TempForce_R(N_orbitals, 0:GridSize-1, 0:GridSize-1, 0:GridSize-1) )

  CALL calc_dens_and_coeffs_r(NZCoeffs, Coeffs_R)

  ! KINETIC TERM
  PreFac = CMPLX(2*PI*PI/BoxL**2, kind=8)

  DO Iorb = 1, N_orbitals
    OrbForce(Iorb,:) = -PreFac*G2Grid*NZCoeffs(Iorb,:)
  ENDDO

  ! EXCHANGE CORRELATION 
  DO I1 = 0, GridSize-1
    DO J1 = 0, GridSize-1
      DO K1 = 0, GridSize-1
        TempVec_R(I1, J1, K1) = CMPLX(Vxc(I1, J1, K1), kind=8)
      ENDDO
    ENDDO
  ENDDO
  CALL BackWard_FFT(GridSize, TempVec_R, TempVec_K)

  ! LOCAL PSEUDOPOTENTIAL
  TempVec_K = TempVec_K + PseudoGrid

  ! HARTREE
  PreFac = BoxL**2/PI
  TempVec_K = TempVec_K + PreFac*Density_K*Gmin2Grid
  
  ! Calculate contrbution from local potential
  CALL Forward_FFT(GridSize, TempVec_K, TempVec_R)
  DO N = 1, N_orbitals
    TempForce_R(N,:,:,:) = TempVec_R*Coeffs_R(N,:,:,:)*SQRT(Omega)
    CALL Backward_FFT(GridSize, TempForce_R(N,:,:,:), TempForce_K(N,:,:,:))
  ENDDO
  DO IIndex = 1, NoOfPW
    I1 = GridPos(IIndex,1)
    J1 = GridPos(IIndex,2)
    K1 = GridPos(IIndex,3)
    DO N=1, N_orbitals
      OrbForce(N, IIndex) = OrbForce(N, IIndex) - TempForce_K(N, I1, J1, K1)
    ENDDO
  ENDDO

  ! NONLOCAL PSEUDOPOTENTIAL
  DO Iorb = 1, N_orbitals
    DO N=1, N_ion
      AtomNum = Ions(N)%AtomNum
      IF( AtomNum > 4 ) THEN
        PP_Index = Get_index_pp(AtomNum)
        h_1s = PP_Params(PP_Index)%h_1s
        IF( PP_Params(PP_Index)%MaxL > 0 ) THEN
          h_2s = PP_Params(PP_Index)%h_2s
          h_1p = PP_Params(PP_Index)%h_1p
        ENDIF
        F = 0.D0
        DO IIndex=1, NoOfPW
          F = F + NonLocal(IIndex,N,1)*CONJG(NZCoeffs(Iorb,IIndex))
        ENDDO
        DO IIndex=1, NoOfPW
          OrbForce (Iorb, IIndex) = OrbForce (Iorb, IIndex) - &
             CONJG(F)*h_1s*NonLocal(IIndex,N, 1)
        ENDDO
        IF( PP_Params(PP_Index)%MaxL > 0 ) THEN
          DO L = 2,5
            IF( L == 2 ) THEN
              hfac = h_2s
            ELSE
              hfac = h_1p
            ENDIF
            F = 0.D0
            DO IIndex = 1, NoOfPW
              F = F + NonLocal(IIndex,N,L)*CONJG(NZCoeffs(Iorb,IIndex))
            ENDDO
            DO IIndex = 1, NoOfPW
              OrbForce(Iorb, IIndex) = OrbForce(Iorb, IIndex) - &
                 CONJG(F)*hfac*NonLocal(IIndex,N,L)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  
  DO N = 1,N_orbitals
    OrbForce(N,:) = OrbForce(N,:)*FillFac(N)
  ENDDO

END SUBROUTINE


