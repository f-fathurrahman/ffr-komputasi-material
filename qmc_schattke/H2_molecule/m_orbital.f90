!---------------
module m_orbital
!---------------

!   A p-wave extension has been introduced replacing the  
!   linear 2s-wave term in the single particle product ansatz

use m_highlevel
use m_midlevel

implicit none

public :: ORBWAV,ORBDER,ORBLCAOWAV,ORBLCAODER,ORBPRODWAV
public :: ORBPRODDER,ORBPRODPWAV,ORBPRODPDER
integer, parameter, public :: NORBM=1 ! max number of orbital
integer, public :: NORB
real(dp), public :: CKPOINT
real(dp), public :: WAVEC,ALPHA
character(40), public :: ORBNAME
       

  contains !----------------------------------



!-----------------------
subroutine ORBWAV(r, psi)
!-----------------------
  implicit none
  real(dp), intent(in), dimension(3) :: r
  real(dp), intent(out), dimension(NORBM,NatomsMax) :: psi
  !
  select case(ORBNAME) 
  !
  case("orbital composition from LCAO           ")
    call ORBLCAOWAV(r,psi)
  !
  case("orbital composition from product        ")
    call ORBPRODWAV(r,psi)
  !
  case("orbital composition from product_p      ")
    call ORBPRODPWAV(r,psi)
  !
  case default
    call ORBLCAOWAV(r,psi)
  !
  end select
end subroutine ORBWAV


!-------------------------------
subroutine ORBDER(r,psi,pgr,pla)
!-------------------------------
  implicit none
  real(dp),intent(in),dimension(3) :: r
  real(dp),intent(out),dimension(NORBM,NatomsMax)   :: psi,pla
  real(dp),intent(out),dimension(3,NORBM,NatomsMax) :: pgr
  ! local variables
  select case (ORBNAME) 
    !
    case("orbital composition from LCAO           ")
      call ORBLCAODER(r,psi,pgr,pla)
    !
    case("orbital composition from product        ")
      call ORBPRODDER(r, psi, pgr, pla)
    !
    case("orbital composition from product_p      ")
      call ORBPRODPDER(r, psi, pgr, pla)
    !
    case default
      call ORBLCAODER(r, psi, pgr, pla)
    !
  end select
end subroutine ORBDER


!---------------------------
subroutine ORBLCAOWAV(r,psi)
!---------------------------
  implicit none
! Calculate orbital phi=(1+c*x)exp(-alpha*r) and wavefunction psi as 
! linear combination of phi.
! Up spin runs from IE=1,NelectronsPerSpin, down spin from IE=NelectronsPerSpin+1,Nelectrons
! Should be tabulated later when we use determinants!
  real(dp), intent(in), dimension(3) :: r
  real(dp), intent(out), dimension(NORB,Natoms) :: psi
  integer :: k
  real(dp) :: s,sx
  real(dp), dimension(3) :: rad
  real(dp), dimension(Natoms) :: phi
  
  do k=1,Natoms
    rad(1:3) = r(1:3) - RK(1:3,k)
    s = sqrt(rad(1)**2 + rad(2)**2 + rad(3)**2)
    sx = rad(1)
    ! Single orbitals look like (note: orbitals differ only by shifting)
    phi(k) = (1.0d0+WAVEC*sx)*dexp(-ALPHA*s)
  enddo
  ! Associate single particle wavefunction with electron, here 1 to 1
  ! psi(1:NORB,1:Natoms)=phi(1:Natoms) !single center orbital
  psi(1:NORB,1) = phi(1) + CKPOINT*phi(2) ! LCAO
  psi(1:NORB,2) = CKPOINT*phi(1) + phi(2) ! LCAO
end subroutine ORBLCAOWAV



!-----------------------------------
subroutine ORBLCAODER(r,psi,pgr,pla)
!-----------------------------------
  implicit none
  ! Gradient, Laplacian, and wave function
  real(dp), intent(in), dimension(3) :: r
  real(dp), intent(out), dimension(NORB,Natoms) :: psi,pla
  real(dp), intent(out), dimension(3,NORB,Natoms) :: pgr
  integer :: i,k
  real(dp) :: s,sx
  real(dp),dimension(3) :: rad,rx
  real(dp) :: w
  real(dp),dimension(Natoms) :: phi,lap
  real(dp),dimension(3,Natoms) :: grp
  
  rx(1:3) = 0.d0
  rx(1) = 1.d0
  do k = 1,Natoms
    rad(1:3) = r(1:3) - RK(1:3,k)
    s = sqrt(sum(rad(1:3)**2))
    sx = rad(1)
    ! Single orbitals look like (note: orbitals differ only by shifting)
    w = 1.0_dp+WAVEC*sx
    phi(k) = w*dexp(-ALPHA*s)
    grp(1:3,k) = rx(1:3)*WAVEC/w - rad(1:3)/s*ALPHA
    lap(k) = ALPHA**2 - 2.0_dp/s*ALPHA*(1+WAVEC*sx/w)
  enddo
  
  ! Associate single particle wavefunction with electron, here 1 to 1
  psi(1:NORB,1) = phi(1) + CKPOINT*phi(2) ! LCAO
  psi(1:NORB,2) = CKPOINT*phi(1) + phi(2) ! LCAO
  do i = 1,3
    pgr(i,1:NORB,1) = (grp(i,1)*phi(1) + CKPOINT*grp(i,2)*phi(2))/psi(1:NORB,1)
    pgr(i,1:NORB,2) = (CKPOINT*grp(i,1)*phi(1) + grp(i,2)*phi(2))/psi(1:NORB,2)
  enddo
  
  pla(1:NORB,1) = (lap(1)*phi(1) + CKPOINT*lap(2)*phi(2))/psi(1:NORB,1)
  pla(1:NORB,2) = (CKPOINT*lap(1)*phi(1) + lap(2)*phi(2))/psi(1:NORB,2)
end subroutine ORBLCAODER


!---------------------------
subroutine ORBPRODWAV(r,psi)
!---------------------------
  implicit none
! Calculate orbital phi=exp(-alpha*r) and wavefunction psi as 
! product of phi.Only one orbital.
! Up spin runs from IE=1,NelectronsPerSpin, down spin from IE=NelectronsPerSpin+1,Nelectrons
! Should be tabulated later when we use determinants!
  real(dp),intent(in), dimension(3) :: r
  real(dp),intent(out), dimension(NORB,Natoms) :: psi
  integer :: k
  real(dp) :: s
  real(dp),dimension(3) :: rad
  real(dp),dimension(NORB,Natoms) :: phi
  do k=1,Natoms
    rad(1:3)=r(1:3)-RK(1:3,k)
    s = dsqrt (sum (rad(1:3)**2)) 
    ! Single orbitals look like (note: orbitals differ only by shifting)
    ! The extension by WAVEC .NEQ. 0.0 seems useless as covered by change 
    ! of ALPHA, for small values at least
    phi(1:NORB,k)=(1.0_dp+WAVEC*s)*dexp(-ALPHA*s)
  enddo
  ! Associate single particle  wavefunction with electron, here 1 to 1
  psi(1:NORB,1) = phi(1:NORB,1)*phi(1:NORB,2) ! product ansatz
  psi(1:NORB,2) = phi(1:NORB,1)*phi(1:NORB,2) ! product ansatz
end subroutine ORBPRODWAV



!-----------------------------------
subroutine ORBPRODDER(r,psi,pgr,pla)
!-----------------------------------
  implicit none
  ! Gradient, Laplacian, and wave function
  ! NORB=1, only 1st component is affected
  real(dp),intent(in),dimension(3) :: r
  real(dp),intent(out),dimension(NORB,Natoms) :: psi,pla
  real(dp),intent(out),dimension(3,NORB,Natoms) :: pgr
  integer :: i,k
  real(dp) :: sk,ww
  real(dp),dimension(3)         :: rad
  real(dp),dimension(3,Natoms)      :: w
  real(dp),dimension(NORB)      :: wwave
  real(dp),dimension(NORB,Natoms)   :: phi,lap
  real(dp),dimension(3,NORB,Natoms) :: grp
  
  do k=1,Natoms
    rad(1:3) = r(1:3)-RK(1:3,k)
    sk = sqrt( sum(rad(1:3)**2) )
    ! Single orbitals look like (note: orbitals differ only by shifting)
    if( sk < EMACH ) then
      w(1:3,k) = 0.d0
      sk = 2.d0*EMACH/3.d0
    else
      w(1:3,k)=rad(1:3)/sk
    endif
    wwave(1:NORB) = (1.0_dp + WAVEC*sk)
    phi(1:NORB,k) = wwave(1:NORB)*dexp(-ALPHA*sk)
    do i=1,NORB
      grp(1:3,i,k) = w(1:3,k)*(-ALPHA+WAVEC/wwave(i))
    enddo
    lap(1:NORB,k) = ALPHA**2 - 2.0_dp/sk*ALPHA + 2.0_dp*WAVEC/wwave(1:NORB)*(-ALPHA + 1.0_dp/sk)
  enddo
  ww = dot_product(grp(1:3,1,1),grp(1:3,1,2))
  ! Associate single particle wavefunction with electron
  psi(1,1) = phi(1,1)*phi(1,2) ! product ansatz
  psi(1,2) = phi(1,1)*phi(1,2) ! product ansatz
  pgr(1:3,1,1) = grp(1:3,1,1) + grp(1:3,1,2)
  pgr(1:3,1,2) = grp(1:3,1,1) + grp(1:3,1,2)
  pla(1,1) = lap(1,1) + lap(1,2) + 2.d0*ww
  pla(1,2) = lap(1,1) + lap(1,2) + 2.d0*ww
end subroutine ORBPRODDER


!----------------------------
subroutine ORBPRODPWAV(r,psi)
!----------------------------
  implicit none
  ! Calculate orbital phi=exp(-alpha*r) and wavefunction psi as 
  ! product of phi. 2p_x-wave instead of 2s-wave 
  ! Up spin runs from IE=1,NelectronsPerSpin, down spin from IE=NelectronsPerSpin+1,Nelectrons
  ! Should be tabulated later when we use determinants!
  real(dp), intent(in), dimension(3) :: r
  real(dp), intent(out), dimension(NORB,Natoms) :: psi
  integer :: k
  real(dp) :: s,sx
  real(dp), dimension(3) :: rad
  real(dp), dimension(NORBM,NatomsMax) :: phi
  !
  do k=1,Natoms
    rad(1:3) = r(1:3) - RK(1:3,k)
    s = sqrt(rad(1)**2 + rad(2)**2 + rad(3)**2) 
    sx = rad(1)
    ! Single orbitals look like (note: orbitals differ only by shifting)
    ! The extension by WAVEC .NEQ. 0.0 seems useless as covered by change 
    ! of ALPHA, for small values at least
    phi(1:NORB,k)=(1.0_dp+WAVEC*sx)*dexp(-ALPHA*s)
  enddo
  ! Associate single particle  wavefunction with electron, here 1 to 1
  psi(1,1)=phi(1,1)*phi(1,2) ! product ansatz
  psi(1,2)=phi(1,1)*phi(1,2) ! product ansatz
end subroutine ORBPRODPWAV



!------------------------------------
subroutine ORBPRODPDER(r,psi,pgr,pla)
!------------------------------------
  implicit none
  ! Gradient, Laplacian, and wave function
  real(dp), intent(in), dimension(3) :: r
  real(dp), intent(out), dimension(NORB,Natoms) :: psi,pla
  real(dp), intent(out), dimension(3,NORB,Natoms) :: pgr
  integer :: i,k
  real(dp) :: sk
  real(dp), dimension(3) :: rad,sx
  real(dp), dimension(3,Natoms) :: w
  real(dp), dimension(NORB) :: ww,wwave
  real(dp), dimension(NORB,Natoms) :: phi,lap
  real(dp), dimension(3,NORB,Natoms) :: grp
  
  do k=1,Natoms
    sx(1:3) = 0.d0
    rad(1:3) = r(1:3) - RK(1:3,k)
    sk = sqrt(sum(rad(1:3)**2))
    sx(1) = rad(1)
    ! Single orbitals look like (note: orbitals differ only by shifting)
    if( sk .lt. EMACH ) then
      w(1:3,k) = 0.d0
      sk = 2.d0*EMACH/3.d0
    else
      w(1:3,k) = rad(1:3)/sk
    endif
    !
    wwave(1:NORB) = (1.d0 + WAVEC*sx(1))
    phi(1:NORB,k) = wwave(1:NORB)*exp(-ALPHA*sk)
    grp(1,1:NORB,k) = ( -w(1,k)*ALPHA + WAVEC/wwave(1:NORB) )
    do i=2,3
      grp(i,1:NORB,k)=-w(i,k)*ALPHA
    enddo
    lap(1:NORB,k) = ALPHA**2 - 2.d0/sk*ALPHA - 2.d0*WAVEC/wwave(1:NORB)*ALPHA*sx(1)/sk
  enddo
  !
  ww(1:NORB) = grp(1,1:NORB,1)*grp(1,1:NORB,2) + &
             & grp(2,1:NORB,1)*grp(2,1:NORB,2) + &
             & grp(3,1:NORB,1)*grp(3,1:NORB,2)
  ! Associate single particle wavefunction with electron
  psi(1:NORB,1) = phi(1:NORB,1)*phi(1:NORB,2) ! product ansatz
  psi(1:NORB,2) = phi(1:NORB,1)*phi(1:NORB,2) ! product ansatz
  do i=1,3
    pgr(i,1:NORB,1) = grp(i,1:NORB,1) + grp(i,1:NORB,2)
    pgr(i,1:NORB,2) = grp(i,1:NORB,1) + grp(i,1:NORB,2)
  enddo
  pla(1:NORB,1) = lap(1:NORB,1) + lap(1:NORB,2) + 2.d0*ww
  pla(1:NORB,2) = lap(1:NORB,1) + lap(1:NORB,2) + 2.d0*ww
end subroutine ORBPRODPDER


end module

