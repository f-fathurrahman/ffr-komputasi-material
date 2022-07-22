SUBROUTINE calc_F_ES(ESForce)
  use globals
    IMPLICIT NONE

    complex(8), INTENT(OUT) :: ESForce(N_ion,3)
    complex(8), ALLOCATABLE :: OvrlForce(:,:)
    complex(8)              :: PreFac
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
    
END SUBROUTINE