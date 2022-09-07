!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE VXC_LDA
! ==================
!
! INPUT:
!   rho [real(8), dimension(n, n)] : electronic density.
! ---------
! OUTPUT:
!   vx [real(8), dimension(n, n)] : exchange-energy potential.
!   vc [real(8), dimension(n, n)] : correlation-energy potential.
!   ex [real(8)] : exchange energy.
!   ec [real(8)] : correlation energy.
!
! Given a electronic density rho, it provides the exchange and correlation
! energies and potentials by calling the "exchange" and correlation"
! subroutines (see below) at each point in the 2D mesh.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vxc_lda(rho, vx, vc, Ex, Ec)
  use mesh

  implicit none

  real(8), intent(in)  :: rho(N, N)  ! the density
  real(8), intent(out) :: vx(N, N)   ! the exchange potential
  real(8), intent(out) :: vc(N, N)   ! the correlation potential
  real(8), intent(out) :: Ex, Ec     ! the exchange and corelation energies

  integer :: i, j
  real(8) :: de

  Ex = 0.0d0
  Ec = 0.0d0

  do i = 1, N
    do j = 1, N
      call exchange(rho(i, j), de, vx(i, j))
      Ex = Ex + de*rho(i, j)

      call correlation(rho(i, j), de, vc(i, j))
      Ec = Ec + de*rho(i, j)
    end do
  end do

  Ex = Ex * delta**2
  Ec = Ec * delta**2

end subroutine vxc_lda




!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE EXCHANGE
! ===================
!
! INPUT:
!   ds [real(8)] : electronic density at a given point.
! ---------
! OUTPUT:
!   vx [real(8)] : exchange potential at the given point
!   ex [real(8)] : exchange energy per volume unit. (surface, in 2D)
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine exchange(ds, ex, vx)
  implicit none

  real(8), intent(in)   :: ds
  real(8), intent(out)  :: ex, vx
  
  real(8), parameter :: a_x = -1.06384608107049d0
  real(8), parameter :: MINDEN = 1.0d-15

  real(8) :: dens
  
  ! Sanity check
  dens = max(0.0d0, ds)

  ! If the density is too small, return zero  
  if(dens < MINDEN) then
     ex = 0.0d0; vx = 0.0d0
    return
  endif

  ex = a_x*sqrt(ds)
  vx = 1.5d0*ex

end subroutine exchange



!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE CORRELATION
! ===================
!
! INPUT:
!   ds [real(8)] : electronic density at a given point.
! ---------
! OUTPUT:
!   vx [real(8)] : exchange potential at the given point
!   ex [real(8)] : exchange energy per volume unit. (surface, in 2D)
!
! Correlation energy per particle and potentials for a homogeneous electron
! gas in 2D, as parametrized by Attacalite et al.
! Refs: [1] C. Attacalite et al, Phys. Rev. Lett. 88, 256601 (2002) for 2D.
!       [2] C. Attacalite, PhD thesis. 
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine correlation(ds, ec, vc)
  implicit none

  real(8), intent(in)   :: ds
  real(8), intent(out)  :: ec, vc
  
  real(8), parameter :: &
     a = -0.1925d0, b = 0.0863136d0, c = 0.0572384d0, &
     e = 1.0022d0, f = -0.02069d0, g = 0.33997d0, h = 1.747d-2, &
     d = -a*h
  real(8), parameter :: pi = 3.141592653589793d0

  real(8), parameter :: MINDEN = 1.0d-15
  real(8) :: rs, dens

  ! Sanity check
  dens = max(0.0d0, ds)

  ! If the density is too small, return zero  
  if(dens < MINDEN) then
    ec = 0.0d0; vc = 0.0d0
    return
  endif

  ! Wigner radius
  rs = sqrt(1.0d0 / (pi*dens))

  ! In unpolarized cases, expressions are fairly simple.
  ec = alpha()
  vc = ec - 0.5d0*rs*dalphadrs()

  return

contains

  real(8) function alpha() ! Eq.[1]4

    alpha = a + (b*rs + c*rs**2 + d*rs**3) * &
       log(1.0d0 + 1.0d0/(e*rs + f*sqrt(rs)**3 + g*rs**2 + h*rs**3) )

  end function alpha

  real(8) function dalphadrs() ! Eq.[2]C3
    real(8) :: efe, efep, lg, x

    efe  = e*rs + f*sqrt(rs)**3 + g*rs**2 + h*rs**3 ! Eq.[2]C5
    efep = e + 1.5d0*f*sqrt(rs) + 2.0d0*g*rs + 3.0d0*h*rs**2 ! Eq. [2]C6
    lg = log(1.0d0 + 1.0d0/efe)
    x  = ((b*rs + c*rs**2 + d*rs**3)*efep)/(efe**2+efe)
    dalphadrs = (b + 2.0d0*c*rs + 3.0d0*d*rs**2)*lg - x ! Eq.[2]C3
  end function dalphadrs

end subroutine correlation




!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE FXC_LDA
! ==================
!
! INPUT:
! n [real(8)] : electronic density at a given point.
! ---------
! OUTPUT:
! fxc [real(8)] : exchange-and correlation term, eq. [whatever]
!
! NOTES:
! ------
! Correlation term is missing!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/!
subroutine fxc_lda(n, fxc)
  real(8), intent(in)  :: n
  real(8), intent(out) :: fxc

  real(8), parameter :: a_x = -1.06384608107049d0
  real(8), parameter :: MINDEN = 1.0d-15
  real(8) :: dens

  ! Sanity check
  dens = max(0.0d0, n)

  ! If the density is too small, return zero  
  if(dens < MINDEN) then
    fxc = 0.0d0
    return
  endif

  fxc = 0.75d0*a_x/sqrt(dens)

end subroutine fxc_lda


