module m_highlevel

! Lithium version
! 25.9.2009: Uses new module strategy
! Here: quantities as parameters, variables on highest program level
implicit none

! Double Precision
integer, parameter,public :: dp=selected_real_kind(2*precision(1.0))
real(dp), parameter       :: EMACH=1.0d-6,PI=3.1415927

! Physical constants
real(dp),parameter,public :: HARTREE=27.21168_dp, BOHR=0.52917706_dp

! NE number electrons
! NES(2) number electrons per spin, 1=up-spin, 2=down-spin
! NK number of nuclei
integer,parameter,public :: NEMAX=3,NKMAX=1,NE=3,NK=1
integer,dimension(2),public  :: NES=(/2,1/)

end module


!-----------------------------------------------------------------------
module m_midlevel
!-----------------------------------------------------------------------

use m_highlevel

! Here: quantities may change during run
implicit none

SAVE
! IE=1,NE index of electrons
! IES=1,2 index of spin: 1=up, 2=down
! IEES=1,NES(IES) index of electrons of same spin
! IENS=NES(IES) no. of electrons with specified spin
! IK=1,NK index of nuclei
integer,public :: IE,IES,IK,IEES,IENS,IMC,IMCR,IMCA,IBLOCKA, MCOUNT
logical,public :: MCSTEP,MCRUN
!  RE position array of electron
!  RNEW position of electron after move
!  RK position array of nuclei
!  DIST actual distances between electrons,4th component for modulus
!  DISTNEW updated distances from a moved electron
!  JASEMACH lowest distance used in Jastrow factor for finiteness
  real(dp),parameter,public :: JASEMACH=2.D0/3.D0*EMACH
  real(dp),public                :: STEPMAX
  real(dp),dimension(3,NE),public :: RE,RENEW
  real(dp),dimension(3,NK),public :: RK
  real(dp),dimension(4,NE,NE),public :: DIST,DISTNEW
  real(kind=dp),public,dimension(NE,NE,2)  :: PSIMAT
  public :: RDIST

contains

subroutine RDIST(r,rn)
       real(kind=dp),intent(in),dimension(3,NE) :: r,rn
        integer :: k,n
        real(dp) :: woo,won
       ielek:do k=1,NE
         if (k .eq. IE) then
          cycle ielek
         end if
         woo=0.0_dp
         won=0.0_dp
         do n=1,3
          DIST(n,IE,k)=r(n,IE)-r(n,k)
!          DIST(n,k,IE)= -DIST(n,IE,k)
          DISTNEW(n,IE,k)=rn(n,IE)-r(n,k)
!          DISTNEW(n,k,IE)=-DISTNEW(n,IE,k)
          woo=woo+DIST(n,IE,k)**2
          won=won+DISTNEW(n,IE,k)**2
         end do
         woo = dsqrt(woo)
         won = dsqrt(won)
!  Cut-off at small distance for Jastrow factor
         if (woo .lt. EMACH)  woo = JASEMACH
         if (won .lt. EMACH)  won = JASEMACH
         DIST(4,IE,k) = woo
!         DIST(4,k,IE) = woo
         DISTNEW(4,IE,k) = won
!         DISTNEW(4,k,IE) = won
       end do ielek
end subroutine RDIST


end module
