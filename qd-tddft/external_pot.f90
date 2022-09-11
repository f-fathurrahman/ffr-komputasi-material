!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE EXTERNAL_POT
! =======================
!
! OUTPUT:
!   v [real(8), dimension(n, n)] : The variable where the external potential is put.
!
! This subroutine places in the "v" argument the external potential that
! defines the quantum dot.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine external_pot(v)
  use mesh
  implicit none
  real(8), intent(out) :: v(N, N)    ! the external potential
  integer :: ix, iy
  real(8) :: r2, omega

  v = 0.0
  ! This defines a harmonic potential:
  omega = 0.22d0
  do ix = 1, N
  do iy = 1, N
     r2 = x(ix, iy)**2 + y(ix, iy)**2
     v(ix, iy) = 0.5d0 * omega**2 * r2 ! - 2.d0 ! offset to get bound states
  enddo
  enddo

  write(*,*) 'Initialized external pot, harmonic with omega = ', omega

end subroutine external_pot
