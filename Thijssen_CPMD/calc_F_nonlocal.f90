SUBROUTINE calc_F_nonlocal(NZCoeffs,NonLocalForce)
  use globals
     IMPLICIT NONE
     complex(8), INTENT(OUT) :: NonLocalForce(N_ion,3)
     complex(8), INTENT(IN)  :: NZCoeffs(N_orbitals, NoOfPW)
     complex(8)              :: F(N_ion), dFdRI(N_ion,3)
     real(8)            :: h(5), prefac
     INTEGER                     :: N, M, K, Iorb, MaxM, AtomNum, &
                                    PP_Index, GCnt, G
  ! Functions
  integer :: get_index_pp


     prefac = 2*PI/BoxL
     NonLocalForce = CMPLX(0.D0)
     DO N = 1, N_ion
       DO Iorb = 1, N_orbitals  
         AtomNum = Ions(N)%AtomNum
         IF (AtomNum > 4) THEN
           PP_Index = get_index_pp(AtomNum)
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
END SUBROUTINE