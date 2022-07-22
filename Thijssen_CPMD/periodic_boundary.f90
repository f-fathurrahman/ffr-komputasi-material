SUBROUTINE periodic_boundary()
  use globals
  IMPLICIT NONE
  INTEGER       :: N, I 
     DO N=1, N_Ion
        DO I=1, 3
          Ions(N)%R_I(I) =  Ions(N)%R_I(I) + BoxL/2 - & 
              BoxL*DNINT(Ions(N)%R_I(I)/BoxL) 
        END DO
     END DO
END SUBROUTINE
