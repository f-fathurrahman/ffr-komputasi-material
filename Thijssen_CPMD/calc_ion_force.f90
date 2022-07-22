SUBROUTINE calc_ion_force(NZCoeffs, IonForce)
  use globals
! The force on the Ions are calculated 
    IMPLICIT NONE
    complex(8), INTENT(IN)  :: NZCoeffs(N_orbitals, NoOfPW)
    complex(8), INTENT(OUT) :: IonForce(N_ion,3)
    complex(8), ALLOCATABLE :: ESForce(:,:), LocalForce(:,:), &
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
END SUBROUTINE