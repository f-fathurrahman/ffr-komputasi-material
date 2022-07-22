SUBROUTINE calc_F_local(LocalForce)
  use globals
  IMPLICIT NONE
  complex(8)              :: PreFac
  INTEGER                     :: N, M
  
  complex(8), INTENT(OUT) :: LocalForce(N_ion,3)
  PreFac = 2*PI*Im*BoxL**2
  LocalForce = CMPLX(0.D0)
  DO N=1, N_ion
    DO M=1, 3
         LocalForce(N,M) = LocalForce(N,M) + &
                     preFac*SUM(GGrid(M,:,:,:)*ShortLocal(N,:,:,:)*CONJG(Density_K))
    END DO
  END DO
END SUBROUTINE
