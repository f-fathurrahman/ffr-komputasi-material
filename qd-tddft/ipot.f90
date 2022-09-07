!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE INTERACTION_POT
! ==========================
!
! INPUT:
!   rho [real(8), dimension(n, n)] : electronic density.
! ---------
! OUTPUT:
!   vx [real(8), dimension(n, n)] : Hartree potential.
!   vx [real(8), dimension(n, n)] : exchange-energy potential.
!   vc [real(8), dimension(n, n)] : correlation-energy potential.
!   ex [real(8)] : exchange energy.
!   ec [real(8)] : correlation energy.
!
! Given an input density, this subroutine must provide the Hartree, exchange
! and correlation potential, as well as the exchange energy and potential.
! For that purpose, it must call "poisson_solve" and "vxc_lda"
!
! It is useful if the code is built in such a way that is easy to choose
! between full use of LDA, use of LDA for exchange only (null correlation),
! use only of the Hartree term (exchange and correlation are null), and
! the independent particle approximation (everything is null).
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interaction_pot(rho, vh, vx, vc, ex, ec)
  use mesh
  use poisson
  implicit none
  real(8), intent(in) :: rho(n, n)
  real(8), intent(out) :: vh(n, n), vx(n, n), vc(n, n), ex, ec
  integer :: mode
  integer, parameter :: INDEPENDENT_PARTICLES = 0, &
                        HARTREE               = 1, &
                        HARTREE_X             = 2, &
                        HARTREE_XC            = 3

  vh = 0.0; vx = 0.0; vc = 0.0; ex = 0.0; ec = 0.0

  !mode = HARTREE
  !mode = INDEPENDENT_PARTICLES
  mode = HARTREE_XC

  select case(mode)
    case(INDEPENDENT_PARTICLES)
      vh = 0.0_8; vx = 0.0_8; vc = 0.0_8; ex = 0.0_8; ec = 0.0_8
    case(HARTREE)
      vx = 0.0_8; vc = 0.0_8; ex = 0.0_8; ec = 0.0_8
      call poisson_solve(rho, vh)
    case(HARTREE_X)
      call vxc_lda(rho, vx, vc, Ex, Ec)
      vc = 0.0_8; ec = 0.0_8
      call poisson_solve(rho, vh)
    case(HARTREE_XC)
      call vxc_lda(rho, vx, vc, Ex, Ec)
      call poisson_solve(rho, vh)
  end select

end subroutine interaction_pot





