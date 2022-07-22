SUBROUTINE calc_dens_and_coeffs_r(NZCoeffs, Coeffs_R)
  use globals
  ! Calculates the density_K from Coeffs_K, and normalize all
  IMPLICIT NONE
  complex(8), INTENT(IN) :: NZCoeffs(N_orbitals, NoOfPW)
  complex(8), INTENT(INOUT) :: Coeffs_R(N_orbitals,0:GridSize-1,0:GridSize-1,0:GridSize-1)
  complex(8), ALLOCATABLE :: Coeffs_K(:,:,:,:)
  
  INTEGER :: ElecCnt, IIndex, N
  real(8) :: Norm
  
  ALLOCATE(Coeffs_K(N_orbitals,0:GridSize-1,0:GridSize-1,0:GridSize-1)) 
  Coeffs_K = CMPLX(0.D0)
  DO IIndex = 1, NoOfPW
    Coeffs_K(:, GridPos(IIndex,1), GridPos(IIndex,2), GridPos(IIndex,3)) = &
         NZCoeffs(:,IIndex)
  END DO 
  DO N=1, N_orbitals
    CALL Forward_FFT(GridSize, Coeffs_K(N,:,:,:),Coeffs_R(N,:,:,:))
  END DO
  Coeffs_R = Coeffs_R/SQRT(Omega)
  Density_R = CMPLX(0.D0)
  DO N = 1, N_orbitals
    Density_R = Density_R + FillFac(N)*Coeffs_R(N,:,:,:)*CONJG(Coeffs_R(N,:,:,:))
  END DO
  CALL Backward_FFT(GridSize, Density_R, Density_K)

END SUBROUTINE
