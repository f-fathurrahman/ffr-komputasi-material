!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE SCF
! ==============
!
! scf performs the self-consistent cycle that solves the Kohn-Sham system of
! equation.
!
! It takes no arguments, but acts on the wfs variable defined in states module.
! The Hamiltonian must be built by making use of the external_pot and
! interaction_pot subroutines; the density may be built from the wfs variables
! by calling the build_rho subroutines.
!
! There are a couple of parameters that must be fixed: tol, which fixes
! the threshold below which the problem is considered to be solved, and
! alpha which is the mixing parameter for the linear mixing scheme.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine scf()
  use states
  use mesh
  implicit none

  real(8), external :: energy
  real(8) :: ex, ec, etot, diff
  real(8), parameter :: alpha = 0.2_8, &
                        tol = 1.0e-7_8
  real(8), allocatable :: eigenval(:), residues(:)
  real(8), allocatable :: vext(:,:), vh(:,:), vx(:,:), vc(:,:), vtot(:, :) ! the potentials
  real(8), allocatable :: rho(:,:), rhoprev(:, :), phi(:, :)                             ! the density
  
  integer :: i, j
  integer, parameter :: max_iter = 100
  
  allocate(eigenval(N_wf), residues(N_wf))
  allocate(vext(N, N), vh(N, N), vx(N, N), vc(N, N), vtot(N, N))
  allocate(rho(N, N), rhoprev(N, N), phi(N, N))

  ! The first iteration will use the independent particle approximation to the problem
  call external_pot(vext)
  call output(vext, 'vext.dat')
  vh = 0.0_8; vc = 0.0_8; vx = 0.0_8; ex = 0.0_8; ec = 0.0_8
  vtot = vext + vh + vx + vc

  do i = 1, max_iter

     call conjugate_gradients(N_wf, wfs, vtot, eigenval, residues)

     call external_pot(vext)
     rhoprev = rho
     call build_rho(wfs, rho)

     diff = sqrt(dotproduct(rho-rhoprev,rho-rhoprev))
     if(diff < tol) exit

     rho = alpha*rho + (1.0 - alpha)*rhoprev
     call interaction_pot(rho, vh, vx, vc, ex, ec)
     vtot = vext + vh + vx + vc

     etot =  energy(eigenval, rho, vh, vc, vx, ec, ex)

     write(*, '(/,a,i6)') 'SCF CYCLE ITER # ',i
     write(*, '(5x,a,es18.4)') 'diff = ',diff
     do j = 1, N_wf
        write(*, '(5x,i4, 2es20.8)') j, eigenval(j), residues(j)
     enddo
     write(*,'(a)')

  enddo

  ! Gather all the final numbers...
  call external_pot(vext)
  call build_rho(wfs, rho)
  call interaction_pot(rho, vh, vx, vc, ex, ec)
  vtot = vext + vh + vx + vc
  do j = 1, N_wf
     call hpsi(vtot, wfs(:, :, j), phi)
     eigenval(j) = dotproduct(wfs(:, :, j), phi)
  enddo
  etot =  energy(eigenval, rho, vh, vc, vx, ec, ex)


  ! Write down all the information.
  write(*, '(/,a,/)') 'SCF CYCLE ENDED '
  write(*, '(5x,a,es18.4)') 'diff = ',diff
  write(*, '(5x,a,es18.4)') 'Etot = ', etot
  do j = 1, N_wf
     write(*, '(5x,i4, 2es20.8)') j, eigenval(j), residues(j)
  enddo
  write(*,'(a)')

  ! Output the final functions
  call output(rho, 'rho.dat')
  call output(vh, 'vh.dat')
  call output(vc, 'vc.dat')
  call output(vx, 'vx.dat')

  deallocate(vext, vh, vx, vc, rho, vtot)

end subroutine scf


