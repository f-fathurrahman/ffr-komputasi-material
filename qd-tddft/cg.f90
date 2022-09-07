!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE CONJUGATE_GRADIENTS
! ==============================
!
! INPUT:
!   nst [integer] : number of eigenfunctions/eigenvalues that are seeked.
!   wf [real(8), dimension(n, n, nst)] : the initial guess for the eigenfunctions
!   v [real(8), dimension(n, n)] : the potential that defines the Kohn-Sham
!     Hamiltonian that defines the problem, and which is internally passed to
!     hpsi subroutine to perform the matrix-vector multiplications.
! OUTPUT:
!   wf [real(8), dimension(n, n, nst)] : the resulting approximated eigenvectors.
!   eigenval [real(8), dimension(n)] : the eigenvalues.
!   residues [real(8), dimension(n)] : the residues of the approximations:
!     residues(i) = || H|wf(i)> - <wf(i)|H|wf(i)>|wf(i)> ||
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine conjugate_gradients(nst, wf, v, eigenval, residues)
  use mesh
  implicit none
  integer, intent(in) :: nst
  real(8), intent(inout) :: wf(n, n, nst)
  real(8), intent(in)    :: v(n, n)
  real(8), intent(out)   :: eigenval(nst), residues(nst)

  integer :: matvec, ist, k, niter, i, j
  real(8) :: dump, stheta, stheta2, ctheta, ctheta2, gamma, alpha, beta, theta, theta2, mu, res, tol, lambda, sol(2)
  real(8), allocatable :: psi(:, :), phi(:, :), hcgp(:, :), cg(:, :), cgp(:, :), sd(:, :)
  real(8), parameter :: pi = 3.141592653589793_8

  allocate(psi(n, n), phi(n, n), hcgp(n, n), cg(n, n), cgp(n, n), sd(n, n))

  tol = 1.0e-5
  niter = 5000
  matvec = 0

  states:    do ist = 1, nst

    ! 1. Orthonormalize current vector to previously calculated ones.
    psi(:, :) = wf(:, :, ist)
    do k = 1, ist - 1
       dump = dotproduct(wf(:, :, k), psi(:, :))
       psi(:, :) = psi(:, :) - dump*wf(:, :, k)
    enddo
    dump = sqrt(dotproduct(psi,psi))
    psi(:, :) = psi(:, :)/dump

    ! 2. Initial gradient.
    call hpsi(v, psi, phi)

    ! 3. 
    ctheta = 1.0_8
    stheta = 0.0_8
    hcgp = 0.0_8
    cg   = 0.0_8
    mu   = 1.0_8

    ! 4. Iterations for this band, up to niter iterations maximum.
    do i = 1, niter

       ! 4.1. Get H|phi>, and <phi|H|phi>
       phi = ctheta*phi + stheta*hcgp
       lambda = dotproduct(psi,phi)

       ! 4.2. Check convergence
       residues(ist) = residue(phi, lambda, psi)
       if(residues(ist) < tol) exit

       ! 4.3. Get steepest descent vector.
       sd = lambda*psi - phi
       do k = 1, ist - 1
          dump = dotproduct(wf(:, :, k),sd)
          sd(:, :) = sd(:, :) - dump*wf(:, :, k)
       enddo

       ! 4.4. Get conjugate-gradient vector
       gamma = dotproduct(sd,sd)/mu
       mu = dotproduct(sd,sd)
       cg = sd + gamma*cg

       ! 4.5.
          dump = dotproduct(psi,cg)
          cgp = cg - dump*psi
          dump = sqrt(dotproduct(cgp,cgp))
       cgp = cgp/dump
     
       ! 4.6.
       call hpsi(v, cgp, hcgp)

       ! 4.7.
       alpha = - lambda + dotproduct(cgp,hcgp)
       beta = 2.0_8*dotproduct(cgp,phi)
       theta = 0.5_8*atan(-beta/alpha)
       ctheta = cos(theta)
       stheta = sin(theta)

       ! I am not sure wether this is necessary or not.
       theta2 = theta + pi/2.0_8
       ctheta2 = cos(theta2)
       stheta2 = sin(theta2)
       sol(1) = ctheta**2*lambda + stheta**2*dotproduct(cgp,hcgp) + 2.0_8*stheta*ctheta*dotproduct(cgp,phi)
       sol(2) = ctheta2**2*lambda + stheta2**2*dotproduct(cgp,hcgp) + 2.0_8*stheta2*ctheta2*dotproduct(cgp,phi)

       if(sol(2) < sol(1)) then
          theta = theta2
          stheta = stheta2
          ctheta = ctheta2
       endif

       ! 4.8.
       psi = ctheta*psi + stheta*cgp

    enddo

    ! Put the generated result into the original variables...
    wf(:, :, ist) = psi(:, :)
    eigenval(ist) = lambda

  enddo states

  niter = matvec
 
  deallocate(psi, phi, hcgp, cg, cgp, sd)

  contains

  real(8) function residue(hf, e, f) result (r)
       real(8), intent(in) :: f(:, :), e, hf(:, :)
       real(8), allocatable :: res(:, :)
       allocate(res(n, n))
       res = hf - e*f
       r = sqrt(dotproduct(res,res))
       deallocate(res)
  end function residue

end subroutine conjugate_gradients
