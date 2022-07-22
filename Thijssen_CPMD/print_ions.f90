SUBROUTINE print_ions()
  use globals
  IMPLICIT NONE
  INTEGER I
  OPEN (UNIT=11, FILE='Ions.dat', POSITION='APPEND')
  DO I= 1, N_Ion
        print '(A23 I3)', 'Ion Nr:', I
        print '(A23 F15.8 F15.8 F15.8)', 'Ion R_I:', Ions(I)%R_I
        print *
        WRITE (11,'(3F15.8 $)') Ions(I)%R_I
  END DO
  WRITE (11,*)
  CLOSE (11) 
END SUBROUTINE