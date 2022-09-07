subroutine propagate()
  use states
  use mesh
  implicit none

  real(8), external :: energy
  complex(8), allocatable :: hwf(:, :)
  real(8), allocatable :: rho(:, :)
  real(8), allocatable :: vext(:,:), vh(:,:), vx(:,:), vc(:,:), vtot(:, :) ! the potentials
  real(8), allocatable :: eigenval(:)

  integer :: niter, i, j, ix, iy, dipole_unit
  real(8) :: dt, t, ex, ec, etot, dipole(1:2), prop_time

  allocate(rho(N, N), vext(N, N), vh(N, N), vx(N, N), vc(N, N), vtot(N, N), hwf(N, N))
  allocate(eigenval(N_wf))

  prop_time = 1000.d0
  dt = 1.0d0
  !dt = 0.5d0
  niter = nint(prop_time/dt)

  dipole_unit = 12
  open(unit = dipole_unit, file="dipole.dat", action ="write", status="replace", form ="formatted")
  close(unit = dipole_unit)
  
  write(*,'(/,a,/)') 'Time dependent propagation follows.'
  do i = 0, niter, 1

     call external_pot(vext)
     call zbuild_rho(zwfs, rho)
     call interaction_pot(rho, vh, vx, vc, ex, ec)
     vtot = vext + vh + vx + vc
     do j = 1, N_wf
        call zhpsi(vtot, zwfs(:, :, j), hwf(:, :))
        eigenval(j) = zdotproduct(zwfs(:, :, j), hwf)
     enddo
        
     etot =  energy(eigenval, rho, vh, vc, vx, ec, ex)

     dipole = 0.0_8
     do ix = 1, n
        do iy = 1, n
           dipole(1) = dipole(1) + x(ix, iy)*rho(ix, iy)
           dipole(2) = dipole(2) + y(ix, iy)*rho(ix, iy)
        enddo
     enddo

     t = i*dt

     open(unit = dipole_unit, file="dipole.dat", position="append")
     write(dipole_unit, '(3e18.8)') t, dipole(1), dipole(2)
     close(unit = dipole_unit)

     write(*, *) i, t, etot

     call propagator(vtot, dt)

  enddo

  close(unit=dipole_unit)
  
  deallocate(rho, vext, vh, vx, vc, vtot, hwf, eigenval)

end subroutine propagate
