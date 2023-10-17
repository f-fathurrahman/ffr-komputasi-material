C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      module orbital
C   Lithium version with optimization, uses reference SLAPOPT and
C   varying SLAP. Therefore additional routines INITORB1 and ORBWAV1
C   are introduced which contain SLAPOPT instead of SLAP
       use highlevel
       use midlevel
       implicit none
C   Introduces constants for HF effective potential orbitals
C   from [Sarsa,Galvez,Buendia].
C   Calculates single electron wavefunctions and its derivatives.
C   Slater parameter SLAP compress wavefunction (1s,2s only)
C   according to exp(zr)-->exp(zr/SLAP).
C
C   MPHI=number of Slater functions
C   NORB=number of orbitals, 3 only for testing
C   NELORB(NEMAX)=maps the electron index onto the set of orbitals
       integer,parameter,public          :: MPHI=7,NORB=2
       integer,dimension(NE),public      :: NELORB
       real(dp),public                   :: SLAP0,SLAPD
       real(kind=dp),public,dimension(NORB)  :: SLAP,PSINEW
       real(kind=dp),public,dimension(NORB)  :: SLAPOPT,PSINEWOPT
C   Internal quantities
       integer,private,dimension(MPHI)            :: phinp,fak
       real(kind=dp),private,dimension(MPHI)      :: phizeta,norm
       real(kind=dp),private,dimension(MPHI,NORB) :: phiaa
       data phinp /2*1,4*2,3/
       data phizeta / 0.72089388_dp, 2.61691643_dp, 0.69257443_dp,
     &  1.37137558_dp, 3.97864549_dp, 13.52900016_dp, 19.30801440_dp/
       data phiaa(1:MPHI,1)/ -0.12220686_dp, 1.11273225_dp,
     &  0.04125378_dp, 0.09306499_dp, -0.10260021_dp, -0.00034191_dp,
     &  0.00021963_dp /
       data phiaa(1:MPHI,2)/ 0.47750469_dp, 0.11140449_dp,
     &  -1.25954273_dp, -0.18475003_dp, -0.02736293_dp, -0.00025064_dp,
     &  0.00057962_dp /
C   Specify variables
       public :: INITORB,ORBWAV,ORBDER,INITORB1,ORBWAV1
      contains
C-----------------------------------------------------------------------
       subroutine INITORB
C  1. Maps the electron index onto the set of orbitals through the
C     index array i=NELORB(k), k=electron index, i=orbital index via
C     PSINEW(i) depending only on the orbital index i.
C     For Li: up-spin electron IE=1, IEES=1 occupies 1st orbital
C             up-spin electron IE=2, IEES=2 occupies 2nd orbital
C           down-spin electron IE=3, IEES=1 occupies 1st orbital
C  2. Calculates the factorial, i.e. (2*phinp)! and normalization
C     coefficients.
C  3. Calculates initial wavefunction matrix PSIMAT(i,k,s) with
C     i=orbital index, k= electron index, s= spin index
C
        integer                   :: j,jj,ies,hi
        integer,dimension(NE)     :: hnelorb
        real(dp)                  :: hh,r
        real(dp),dimension(NORB)  :: psi
        data hnelorb / 1, 2, 1 /
C  1. Index array to associate each electron with a wavefunction
        NELORB(1:NE)=hnelorb(1:NE)
C  2. The factorial (2*phinp)! and normalization
        fak=1
        do j=1,MPHI
         do jj=1,2*phinp(j)-1
          fak(j)=fak(j)*(jj+1)
         end do
        end do
        norm=(2.0_dp*phizeta)**(phinp+0.5_dp)/dsqrt(dble(fak))
C  3. The matrix psi_i(k) of the wavefunction is determined.
C     In the columns the index runs over the orbitals i for a
C     fixed electron index k.
C     But one has to discriminate between both spin directions.
C     One begins with spin-up electrons counting them
C     up to the number of electrons NES(1) with this spin.
C     Subsequently follow the spin-down electrons
C     in a new matrix up to NES(2). Take care that the orbital
C     index is correctly associated with the electron index,
C     i.e. for Li we put (El.,Orb.)=(1,1),(2,2) for first matrix and
C     (3,1) for second.
         do j=1,NE
          ies=1
          if (j > NES(1)) ies=2
          hi=(ies-1)*NES(1)
          call ORBWAV(RE(1:3,j),psi)
          do jj=1,NES(ies)
           PSIMAT(jj,j-hi,ies)=psi(NELORB(hi+jj))
          end do
         end do
       end subroutine INITORB
C-----------------------------------------------------------------------
       subroutine ORBWAV(r,psi)
        real(kind=dp),intent(in),dimension(3)      :: r
        real(kind=dp),intent(out),dimension(NORB) :: psi
C   Calculates the Hartree atomic wavefunction of Lithium
C   according to [Sarsa,Galvez,Buendia].
C   Gives only 1s and 2s orbitals (spherical symmetry).
        integer                              :: j
        real(kind=dp)                        :: s
        s=dsqrt(r(1)**2+r(2)**2+r(3)**2)
        psi = 0.0_dp
        do j=1,MPHI
         psi(1:NORB)=psi(1:NORB)+phiaa(j,1:NORB)*s**(phinp(j)-1)*
     &                   norm(j)*dexp(-phizeta(j)*s/SLAP(1:NORB))
        end do
       end subroutine ORBWAV
C-----------------------------------------------------------------------
       subroutine ORBDER(r,psi,pgr,pla)
C   Gradient, laplacian, and wave function from HF effective potential
C   atomic wavefunction via [Sarsa,Galvez,Buendia].
C   (1s and 2s orbitals, spherical symmetry)
        real(dp),intent(in),dimension(3) :: r
        real(dp),intent(out),dimension(NORB)    :: psi,pla
        real(dp),intent(out),dimension(3,NORB)  :: pgr
C
        integer                        :: j
        real(dp),dimension(NORB)       :: psid1,psid2
        real(dp),dimension(MPHI,NORB)  :: pexd1,pexd2,pexd3
        real(kind=dp)                  :: s
        real(dp),dimension(MPHI)       :: c1,c2
        s=dsqrt(r(1)**2+r(2)**2+r(3)**2)
        if (s < EMACH) s=EMACH
        psi = 0.0_dp
        do j=1,MPHI
         pexd1(j,1:NORB)=phiaa(j,1:NORB)*norm(j)*
     &    dexp(-phizeta(j)*s/SLAP(1:NORB))
         psi(1:NORB)=psi(1:NORB)+pexd1(j,1:NORB)*s**(phinp(j)-1)
        end do

        psid1 = 0.0_dp
        do j=1,MPHI
         if (phinp(j) == 1) then
          c1(j) = 0.0_dp
          c2(j) = 0.0_dp
         else
          c1(j) = dble(phinp(j)-1)*s**(phinp(j)-2)
          if (phinp(j) == 2) then
           c2(j)=0.0_dp
          else
           c2(j) = dble((phinp(j)-2)*(phinp(j)-1))*s**(phinp(j)-3)
          end if
         end if
         pexd2(j,1:NORB)=pexd1(j,1:NORB)*phizeta(j)/SLAP(1:NORB)*
     &                   s**(phinp(j)-1)
         pexd3(j,1:NORB)=pexd1(j,1:NORB)*c1(j)
         psid1(1:NORB)=psid1(1:NORB)+pexd3(j,1:NORB)-pexd2(j,1:NORB)
        end do
        psid2 = 0.0_dp
        do j=1,MPHI
         psid2(1:NORB)=psid2(1:NORB)+pexd2(j,1:NORB)*
     &        phizeta(j)/SLAP(1:NORB)
         psid2(1:NORB)=psid2(1:NORB)-
     &        2.0_dp*pexd3(j,1:NORB)*phizeta(j)/SLAP(1:NORB)
         psid2(1:NORB)=psid2(1:NORB) + pexd1(j,1:NORB)*c2(j)
        end do
        do j=1,3
         pgr(j,1:NORB)=r(j)/s*psid1(1:NORB)
        end do
        pla = psid2 + 2.0_dp/s*psid1
       end subroutine ORBDER
C
C-----------------------------------------------------------------------
       subroutine INITORB1
C  1. Maps the electron index onto the set of orbitals through the
C     index array i=NELORB(k), k=electron index, i=orbital index via
C     PSINEW(i) depending only on the orbital index i.
C     For Li: up-spin electron IE=1, IEES=1 occupies 1st orbital
C             up-spin electron IE=2, IEES=2 occupies 2nd orbital
C           down-spin electron IE=3, IEES=1 occupies 1st orbital
C  2. Calculates the factorial, i.e. (2*phinp)! and normalization
C     coefficients.
C  3. Calculates initial wavefunction matrix PSIMAT(i,k,s) with
C     i=orbital index, k= electron index, s= spin index
C
        integer                   :: j,jj,ies,hi
        integer,dimension(NE)     :: hnelorb
        real(dp)                  :: hh,r
        real(dp),dimension(NORB)  :: psi
C
C     The matrix psi_i(k) of the wavefunction is determined.
C     In the columns the index runs over the orbitals i for a
C     fixed electron index k.
C     But one has to discriminate between both spin directions.
C     One begins with spin-up electrons counting them
C     up to the number of electrons NES(1) with this spin.
C     Subsequently follow the spin-down electrons
C     in a new matrix up to NES(2). Take care that the orbital
C     index is correctly associated with the electron index,
C     i.e. for Li we put (El.,Orb.)=(1,1),(2,2) for first matrix and
C     (3,1) for second.
         do j=1,NE
          ies=1
          if (j > NES(1)) ies=2
          hi=(ies-1)*NES(1)
          call ORBWAV1(RE(1:3,j),psi)
          do jj=1,NES(ies)
           PSIMATOPT(jj,j-hi,ies)=psi(NELORB(hi+jj))
          end do
         end do
       end subroutine INITORB1
C-----------------------------------------------------------------------
       subroutine ORBWAV1(r,psi)
        real(kind=dp),intent(in),dimension(3)      :: r
        real(kind=dp),intent(out),dimension(NORB) :: psi
C   Calculates the Hartree atomic wavefunction of Lithium
C   according to [Sarsa,Galvez,Buendia].
C   Gives only 1s and 2s orbitals (spherical symmetry).
        integer                              :: j
        real(kind=dp)                        :: s
        s=dsqrt(r(1)**2+r(2)**2+r(3)**2)
        psi = 0.0_dp
        do j=1,MPHI
         psi(1:NORB)=psi(1:NORB)+phiaa(j,1:NORB)*s**(phinp(j)-1)*
     &                   norm(j)*dexp(-phizeta(j)*s/SLAPOPT(1:NORB))
        end do
       end subroutine ORBWAV1
C
      end module orbital
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
