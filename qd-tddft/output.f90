SUBROUTINE output(f, filename)
  USE mesh
  IMPLICIT NONE 
  REAL(8), INTENT(in) :: f(N, N)    ! function to output
  CHARACTER(len=*), INTENT(in) :: filename ! file to output
  !
  INTEGER, parameter :: iunit = 99
  INTEGER :: ix, iy

  OPEN(iunit, file=trim(filename), status='unknown')
  DO ix = 1, N
    DO iy = 1, N
      write(iunit, *) x(ix, iy), y(ix, iy), f(ix, iy)
    ENDDO 
    write(iunit, '(1x)')
  ENDDO 

  CLOSE(iunit)
END SUBROUTINE 

