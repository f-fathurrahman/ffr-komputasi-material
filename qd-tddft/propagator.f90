subroutine propagator(v, dt)
  use states
  use mesh
  use expo
  implicit none

  real(8), intent(in) :: v(N, N)
  real(8), intent(in) :: dt

  integer :: j
  real(8) :: ex, ec
  real(8), allocatable :: vext(:, :), vx(:, :), vc(:, :), rho(:, :), vh(:, :), vp(:, :)
  complex(8), allocatable :: zwfsp(:, :, :)


  allocate(zwfsp(N, N, N_wf), vp(N, N), vext(N, N), vc(N, N), vh(N, N), vx(N, N), rho(N, N))

  ! This is the first, auxiliary propagation, to get v(t+dt)
  zwfsp = zwfs
  do j = 1, N_occ
     call exponential(v, zwfsp(:, :, j), dt)
  enddo
  call external_pot(vext)
  call zbuild_rho(zwfsp, rho)
  call interaction_pot(rho, vh, vx, vc, ex, ec)
  vp = vext + vh + vx + vc

  do j = 1, N_occ
     call exponential(v, zwfs(:, :, j), dt/2.0_8)
     call exponential(vp, zwfs(:, :, j), dt/2.0_8)
  enddo

  deallocate(zwfsp, vp, vext, vc, rho)
end subroutine propagator
