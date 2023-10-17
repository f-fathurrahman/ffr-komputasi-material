module m_orbital_Li_CR
  ! Lithium version
  ! 25.9.2009: Uses new module strategy
  use m_highlevel
  use m_midlevel
  implicit none

!   Introduce constants for Hartree orbitals [Clementi,Roetti]
!   Calculate single electron wavefunctions and its derivatives
!   Slater parameter SLAP to compress wavefunction (1s,2s only)
!   according to exp(zr)-->exp(zr/SLAP)
! Eckstein: als NDETPAR die Eigenwerte der H-Orbitale 1s und 2s zu Z=3 ?
!   MPHI=number of Slater functions
!   NORB=number of orbitals, 3 only for testing
!   NELORB(NEMAX)=maps the electron index onto the set of orbitals
!   PSINEW(NORB) orbital wavefunction of actual electron
  integer,parameter,public          :: MPHI=6,NORB=3
  integer,dimension(NE),public      :: NELORB
  real(kind=dp),public,dimension(NORB) : SLAP,PSINEW

  ! Internal quantities
  integer,private,dimension(MPHI)            :: phinp
  real(kind=dp),private,dimension(MPHI)      :: phizeta
  real(kind=dp),private,dimension(MPHI,NORB):: phiaa,phiah
  data phinp /2*0,4*1/
  data phizeta / 2.47673_dp,4.69873_dp,0.38350_dp,0.66055_dp,
     &                1.07000_dp,1.63200_dp/
       data phiaa/ 1.97447, 0.639673,
     &            -2.3736E-6, 0.0001294,-0.0008332, 0.0097975,
     &            -0.321705,-0.087116,
     &             0.0001119, 0.11326, 0.042322,-0.12215,
     &    6*1.0_dp/
!   Specify variables
       public :: INITORB,ORBWAV,ORBDER



contains


subroutine INITORB()
!  1. Maps the electron index onto the set of orbitals through the
!     index array i=NELORB(k), i=electron index, i=orbital index
!     For Li: up-spin electron IE=1, IEES=1 occupies 1st orbital
!             up-spin electron IE=2, IEES=2 occupies 2nd orbital
!           down-spin electron IE=3, IEES=1 occupies 1st orbital
!  2. Sets normalized orbital parameters phiah.
!  3. Calculates initial wavefunction matrix PSIMAT(i,k,s)
!     i=orbital index, k= electron index, s= spin index, via
!     PSIMAT(i,k) depending on orbital i for electron k
!
        integer                   :: j,jj,ies,hi
        integer,dimension(NE)     :: hnelorb
        real(dp)                  :: hh,r
        real(dp),dimension(NORB)  :: psi
        data hnelorb / 1, 2, 1 /
!  1. Index array to associate each electron with a wavefunction
        NELORB(1:NE)=hnelorb(1:NE)
!  2. Normalize anew because of SLAP parameter. Original Slater
!     waves are supposed to be normalized. May test by dropping
!     this part that this is not necessary, because Monte-Carlo
!     run is self-normalising.
        do j=1,NORB
         hh=SLAP(j)*dsqrt(SLAP(j))
         do jj=1,MPHI
          phiah(jj,j) = phiaa(jj,j)/(SLAP(j)**phinp(jj)*hh)
         end do
        end do
!  3. The matrix psi_i(k) of the wavefunction is determined.
!     In the columns the index runs over the orbitals i for a
!     fixed electron index k.
!     But one has to discriminate between both spin directions.
!     One begins with spin-up electrons counting them
!     up to the number of electrons NES(1) with this spin.
!     Subsequently follow the spin-down electrons
!     in a new matrix up to NES(2). Take care that the orbital
!     index is correctly associated with the electron index,
!     i.e. for Li we put (El.,Orb.)=(1,1),(2,2) for first matrix and
!     (3,1) for second.
         do j=1,NE
          ies=1
          if (j > NES(1)) ies=2
          hi=(ies-1)*2
          call ORBWAV(RE(1:3,j),psi)
          do jj=1,NES(ies)
           PSIMAT(jj,j-hi,ies)=psi(NELORB(hi+jj))
          end do
         end do
       end subroutine INITORB


subroutine ORBWAV(r,psi)
        real(kind=dp),intent(in),dimension(3)      :: r
        real(kind=dp),intent(out),dimension(NORB) :: psi
!   Calculates the Hartree atomic wavefunction of Lithium
!   according to Clementi & Roetti (1s and 2s orbitals,spherical symmetry)
!   Valid only for phinp=0 or 1!
        integer                              :: j,jj
        real(kind=dp)                        :: s
        real(kind=dp),dimension(MPHI)        :: rhonp
        s=dsqrt(r(1)**2+r(2)**2+r(3)**2)
        do j=1,MPHI
         rhonp(j) = 1.0_dp
         if (phinp(j) .eq. 1) rhonp(j)=s
        end do
        psi = 0.0_dp
        do j=1,MPHI
         psi(1:NORB-1)=psi(1:NORB-1)+phiah(j,1:NORB-1)*rhonp(j)*
     &             dexp(-phizeta(j)*s/SLAP(1:NORB-1))
        end do
        psi(3)=dsin(3.0*r(1))
       end subroutine ORBWAV
C-----------------------------------------------------------------------
       subroutine ORBDER(r,psi,pgr,pla)
C   Gradient, laplacian, and wave function from Hartree atomic wavefct.
C   via Clementi & Roetti (1s and 2s orbitals, spherical symmetry)
C   Valid only for phinp=0 or 1!
        real(dp),intent(in),dimension(3) :: r
        real(dp),intent(out),dimension(NORB)    :: psi,pla
        real(dp),intent(out),dimension(3,NORB)  :: pgr
C
        integer                         :: j
        real(dp),dimension(NORB)       :: psid1,psid2
        real(kind=dp)                   :: s
        real(kind=dp),dimension(MPHI)   :: rhonp
        s=dsqrt(r(1)**2+r(2)**2+r(3)**2)
        do j=1,MPHI
         rhonp(j) = 1.0_dp
         if (phinp(j) .eq. 1) rhonp(j)=s
        end do
        psi = 0.0_dp
        do j=1,MPHI
         psi(1:NORB-1)=psi(1:NORB-1)+phiah(j,1:NORB-1)*rhonp(j)*
     &             dexp(-phizeta(j)*s/SLAP(1:NORB-1))
        end do
        psi(3)=dsin(3.0*r(1))
        psid1 = 0.0_dp
        do j=1,MPHI
         psid1(1:NORB-1)=psid1(1:NORB-1)+
     &                 phiah(j,1:NORB-1)*(dble(phinp(j))
     &                -phizeta(j)*rhonp(j)/SLAP(1:NORB-1))
     &                *dexp(-phizeta(j)*s/SLAP(1:NORB-1))
        end do
        psid1(3)=0.0_dp
        psid2 = 0.0_dp
        do j=1,MPHI
         psid2(1:NORB-1)=psid2(1:NORB-1)
     &    +phiah(j,1:NORB-1)*phizeta(j)/SLAP(1:NORB-1)
     &    *(-2.0_dp*dble(phinp(j))+phizeta(j)*rhonp(j)/SLAP(1:NORB-1))
     &    *dexp(-phizeta(j)*s/SLAP(1:NORB-1))
        end do
        psid2(3)=-9.0*dsin(3.0*r(1))
	    if (s .lt. EMACH) s=EMACH
        do j=1,3
         pgr(j,1:NORB-1)=r(j)/s*psid1(1:NORB-1)
        end do
        pgr(1:3,3)=0.0_dp
        pgr(1,3)=3.0_dp*dcos(3.0*r(1))
        pla = psid2 + 2.0_dp/s*psid1
       end subroutine ORBDER
C
      end module orbital
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
