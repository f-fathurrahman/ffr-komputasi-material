C-----------------------------------------------------------------------
      module highlevel
C   Lithium version
C   25.9.2009: Uses new module strategy
C   Here: quantities as parameters, variables on highest program level
       implicit none
C   Double Precision
       integer, parameter,public ::
     & dp=selected_real_kind(2*precision(1.0))
       real(dp), parameter       :: EMACH=1.0d-6,PI=3.1415927_dp
C   Physical constants
       real(dp),parameter,public :: HARTREE=27.21168_dp,
     &   BOHR=0.52917706_dp
C  NE number electrons
C  NES(2) number electrons per spin, 1=up-spin, 2=down-spin
C  NK number of nuclei
       integer,parameter,public :: NEMAX=3,NKMAX=1,NE=3,NK=1
       integer,dimension(2),public  :: NES=(/2,1/)
      end module highlevel
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      module midlevel
       use highlevel
C  Here: quantities may change during run
       implicit none
       SAVE
C  IE=1,NE index of electrons
C  IES=1,2 index of spin: 1=up, 2=down
C  IEES=1,NES(IES) index of electrons of same spin
C  IENS=NES(IES) no. of electrons with specified spin
C  IK=1,NK index of nuclei
       integer,public :: IE,IES,IK,IEES,IENS,IMC,IMCR,IMCA,IBLOCKA,
     &                   MCOUNT,IMCOPT
       logical,public :: MCSTEP,MCRUN
C  RE position array of electron
C  RNEW position of electron after move
C  RK position array of nuclei
C  DIST actual distances between electrons,4th component for modulus
C  DISTNEW updated distances from a moved electron
C  JASEMACH lowest distance used in Jastrow factor for finiteness
       real(dp),parameter,public :: JASEMACH=2.0_dp/3.0_dp*EMACH
       real(dp),public       :: STEPMAX
       real(dp),public       :: MINOPT,EGUESS,OPTVAR,WPSI,SWEIGHT
       real(dp),dimension(3,NE),public :: RE,RENEW
       real(dp),dimension(3,NK),public :: RK
       real(dp),dimension(4,NE,NE),public :: DIST,DISTNEW
       real(kind=dp),public,dimension(NE,NE,2)  :: PSIMAT,PSIMATOPT
       public :: RDIST
      contains
C-----------------------------------------------------------------------
      subroutine RDIST(r,rn)
       real(kind=dp),intent(in),dimension(3,NE) :: r,rn
        integer :: i,k,n
        real(dp) :: as,woo,won
       ielek:do k=1,NE
         if (k .eq. IE) then
          cycle ielek
         end if
         woo=0.0_dp
         won=0.0_dp
         do n=1,3
          DIST(n,IE,k)=r(n,IE)-r(n,k)
C          DIST(n,k,IE)= -DIST(n,IE,k)
          DISTNEW(n,IE,k)=rn(n,IE)-r(n,k)
C          DISTNEW(n,k,IE)=-DISTNEW(n,IE,k)
          woo=woo+DIST(n,IE,k)**2
          won=won+DISTNEW(n,IE,k)**2
         end do
         woo = dsqrt(woo)
         won = dsqrt(won)
C  Cut-off at small distance for Jastrow factor
         if (woo .lt. EMACH)  woo = JASEMACH
         if (won .lt. EMACH)  won = JASEMACH
         DIST(4,IE,k) = woo
C         DIST(4,k,IE) = woo
         DISTNEW(4,IE,k) = won
C         DISTNEW(4,k,IE) = won
       end do ielek
      end subroutine RDIST
C
C-----------------------------------------------------------------------
      end module midlevel
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
