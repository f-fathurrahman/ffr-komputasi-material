module m_orbital
!
!   Cubic cluster version wihtin infinite walls
!
  use m_highlevel
  use m_midlevel
  implicit none
!
!   Calculates single electron wavefunctions and its derivatives.
!
!   ALPHA= is the exponent coefficient in the Gaussian adjusting of the 
!        wavefunction, usually negative for compression as opposed to
!        the 2-paricle Jastrow factor
!   NELWAVE(NE)= maps the electron index into the set of orbitals
!   KWAVE(3,NE)= gives the momenta of states in the order of the electron 
!              index in the same wayas NELWAVE 
!   PSINEW= is used here and in the main program could be shifted to
!           midlevel variables
       real(dp)                  :: ALPHA=-0.5
       integer,dimension(NE),public         :: NELWAVE
       real(dp),dimension(3,NE),public      :: KWAVE
       real(kind=dp),public,dimension(NE)   :: PSINEW
!
!   Specify variables
       public :: INITORB,ORBWAV,ORBDER,CUBICQD_WAVE
      contains
!-----------------------------------------------------------------------
       subroutine INITORB
!  1. Maps the electron index onto the set of orbitals through the
!     index array NELWAVE(ie), ie=electron index
!     For a jellium cluster of 8 electrons:
!             up-spin electron IE=1, IEES=1 occupies 1st orbital
!             up-spin electron IE=2, IEES=2 occupies 2nd orbital
!             up-spin electron IE=3, IEES=3 occupies 3rd orbital
!             up-spin electron IE=4, IEES=4 occupies 4th orbital
!             down-spin electron IE=5, IEES=1 occupies 1st orbital
!             down-spin electron IE=6, IEES=2 occupies 2nd orbital
!             down-spin electron IE=7, IEES=3 occupies 3rd orbital
!             down-spin electron IE=8, IEES=4 occupies 4th orbital
!
!  2. Initializes the arrays with the coefficients necessary
!     to obtain the single-particle wave functions.
!
!  3. Calculates initial wavefunction matrix PSIMAT(indw,ie,s) with
!     indw=wave index, ie= electron index, s= spin index, and
!     hnelwave(NE) maps electron index ie on wavefunction enumerated
!     by a sequence of single particle states.
!     hkwave(3,NE) maps electron index ie on wavenumber,ik(1:3)=hkwave(k)
        integer                       :: j,jj,ies,hi
        integer,dimension(NE)         :: hnelwave
        real(dp),dimension(NE,3)      :: hkwave
        real(dp),dimension(NE)        :: psi
! data hnelwave /1, 2, 3, 5, 1, 2, 3, 5/
! data hkwave(1:NE,1) /1, 1, 2, 1, 1, 1, 2, 1/
! data hkwave(1:NE,2) /1, 2, 1, 1, 1, 2, 1, 1/
! data hkwave(1:NE,3) /1, 1, 1, 2, 1, 1, 1, 2/
        data hnelwave /1,1/
        data hkwave(1:NE,1) /1,1/
        data hkwave(1:NE,2) /1,1/
        data hkwave(1:NE,3) /1,1/
        NELWAVE(1:NE)=hnelwave(1:NE)
        do j=1,3
         KWAVE(j,1:NE)=PI*hkwave(1:NE,j)/LCLUSTER
        end do
!
!  wavefunctions are:
!                   cos(kx) * cos(ky) * cos(kz) (hnelwave=1)
!                   cos(kx) * sin(ky) * cos(kz) (hnelwave=2)
!                   sin(kx) * cos(ky) * cos(kz) (hnelwave=3)
!                   sin(kx) * sin(ky) * cos(kz) (hnelwave=4)
!                   cos(kx) * cos(ky) * sin(kz) (hnelwave=5)
!                   cos(kx) * sin(ky) * sin(kz) (hnelwave=6)
!                   sin(kx) * cos(ky) * sin(kz) (hnelwave=7)
!                   sin(kx) * sin(ky) * sin(kz) (hnelwave=8)
!  with (kx,ky,kz)=
!                  (KWAVE(1,IE)*x,KWAVE(2,IE)*y,KWAVE(3,IE)*z)
!
!     The matrix psi_{indw}(ie) = PSIMAT(indw,ie,s) of the wavefunction
!     is determined.
!     In the columns the index runs over the wavefunctions for a
!     fixed electron index ie.
!     One has to discriminate between both spin directions.
!     Begin with spin-up electrons counting them
!     up to the number of electrons NES(1) with this spin.
!     Subsequently follow the spin-down electrons
!     in a new matrix up to NES(2).
         do j=1,NE
          ies=1
          if (j > NES(1)) ies=2
          hi=(ies-1)*NES(1)
          call ORBWAV(RE(1:3,j),psi)
          do jj=1,NES(ies)
           PSIMAT(jj,j-hi,ies)=psi(hi+jj)
          end do
         end do
       end subroutine INITORB
!-----------------------------------------------------------------------
       subroutine ORBWAV(r,phi)
        real(kind=dp),intent(in),dimension(3)     :: r
        real(kind=dp),intent(out),dimension(NE)   :: phi
!
!   Calculates single particle standing waves
!	with screening exponential. Maximum number would be NE, if all
!   electrons are associated with different orbitals. When some
!   orbitals are occupied with both spins they are calculated twice
!   for generality only without furhter necessity.
!	Wavefunction stored in phi(NE), wphi1(3), and wphi2 for itself,
!   gradient, and Laplacian, resp..
!   Screening potential vsce=exp(ALPHA*(r/LCLUSTER)^2)
!   internal variables:
        integer                  :: sd,i,j
        real(dp)                 :: hr,vsce
        real(dp),dimension(NE)   :: w2
        real(dp),dimension(3)    :: k,kx
        real(dp),dimension(3,NE) :: w1
!
        sd = 0      ! to jump over the calculation of the derivatives
        hr = sum (r(1:3)**2)/LCLUSTER**2
        vsce = dexp(ALPHA*hr)
        do i=1,NE
         do j=1,3
          k(j) = KWAVE(j,i)
          kx(j) = k(j)*r(j)
         end do
         call CUBICQD_WAVE(sd,NELWAVE(i),k,kx,phi(i),w1(1:3,i),w2(i))
!   include screening exponential
         phi(i) = phi(i)*vsce
        end do
       end subroutine ORBWAV
!-----------------------------------------------------------------------
       subroutine ORBDER(r,wphi1,wphi2)
!   Gradient and Laplacian forthe cubic cluster.
        real(dp),intent(in),dimension(3)       :: r
        real(dp),intent(out),dimension(NE)     :: wphi2
        real(dp),intent(out),dimension(3,NE)   :: wphi1
!   Screening potential vsce=exp(ALPHA*(r/LCLUSTER)^2)
!   Internal variables:
        integer                  :: sd,i,j
        real(dp)                 :: hr,vsce
        real(dp),dimension(3)    :: k,kx
        real(dp),dimension(NE)   :: w,w2,wphi
        real(dp),dimension(3,NE) :: w1
!
        sd = 1     ! for the derivatives
        hr = sum (r(1:3)**2)/LCLUSTER**2
        vsce = dexp(ALPHA*hr)
        do i=1,NE
         do j=1,3
          k(j) = KWAVE(j,i)
          kx(j) = k(j)*r(j)
         end do
         call CUBICQD_WAVE(sd,NELWAVE(i),k,kx,w(i),w1(1:3,i),w2(i))
         wphi(i) = w(i)*vsce
         wphi1(1:3,i) = w1(1:3,i)*vsce + wphi(i)*ALPHA*2._dp*r(1:3)/LCLUSTER**2
         wphi2(i) = w2(i)*vsce + &
     &    dot_product (w1(1:3,i),r(1:3))*vsce*ALPHA*4._dp/LCLUSTER**2 + &
     &    wphi(i)*ALPHA**2*4._dp*hr/LCLUSTER**2 + &
     &    6._dp*wphi(i)*ALPHA/LCLUSTER**2
         end do
       end subroutine ORBDER


!-----------------------------------------------------------------------
       subroutine CUBICQD_WAVE(sd,no,k,kxx,wphi,wphi1,wphi2)
        integer,intent(in)                 :: sd,no
        real(dp),dimension(3),intent(in)   :: k,kxx
        real(dp),intent(out)               :: wphi,wphi2
        real(dp),dimension(3),intent(out)  :: wphi1
!  Provides for the free wavefunctions of a cubic cluster with infinite
!  walls
!  sd = 0/1 switch, =1 for calculating also gradient and Laplacian
!  no = 1,2,3,4,5,6,7,8 index of wavefunction
!  kx(3) = cartesian coordinates times wave vector
!  k(3) = wave vector
!  wphi = wavefunction
!  wphi1(3) = gradient of wavefunction
!  wphi2 = Laplacian of wavefunction
!  internal variables:
        real(dp)    :: a,kx,ky,kz
!
        a = sum ( k(1:3)**2)
        kx = kxx(1)
        ky = kxx(2)
        kz = kxx(3)
        if (no == 1) then
         wphi = cos(kx)*cos(ky)*cos(kz)
         if (sd == 1) then
           wphi1(1) = -k(1)*sin(kx)*cos(ky)*cos(kz)
           wphi1(2) = -k(2)*cos(kx)*sin(ky)*cos(kz)
           wphi1(3) = -k(3)*cos(kx)*cos(ky)*sin(kz)
           wphi2 = -a*wphi
         end if
        end if
        if (no == 2) then
         wphi = cos(kx)*sin(ky)*cos(kz)
         if (sd == 1) then
           wphi1(1) = -k(1)*sin(kx)*sin(ky)*cos(kz)
           wphi1(2) =  k(2)*cos(kx)*cos(ky)*cos(kz)
           wphi1(3) = -k(3)*cos(kx)*sin(ky)*sin(kz)
           wphi2 = -a*wphi
         end if
        end if
        if (no == 3) then
         wphi = sin(kx)*cos(ky)*cos(kz)
         if (sd == 1) then
          wphi1(1) =  k(1)*cos(kx)*cos(ky)*cos(kz)
          wphi1(2) = -k(2)*sin(kx)*sin(ky)*cos(kz)
          wphi1(3) = -k(3)*sin(kx)*cos(ky)*sin(kz)
          wphi2 = -a*wphi
         end if
        end if
        if (no == 4) then
         wphi = sin(kx)*sin(ky)*cos(kz)
         if (sd == 1) then
           wphi1(1) = +k(1)*cos(kx)*sin(ky)*cos(kz)
           wphi1(2) = +k(2)*sin(kx)*cos(ky)*cos(kz)
           wphi1(3) = -k(3)*sin(kx)*sin(ky)*sin(kz)
           wphi2 = -a*wphi
         end if
        end if
        if (no == 5) then
         wphi = cos(kx)*cos(ky)*sin(kz)
         if (sd == 1) then
           wphi1(1) = -k(1)*sin(kx)*cos(ky)*sin(kz)
           wphi1(2) = -k(2)*cos(kx)*sin(ky)*sin(kz)
           wphi1(3) = +k(3)*cos(kx)*cos(ky)*cos(kz)
           wphi2 = -a*wphi
         end if
        end if
        if (no == 6) then
         wphi = cos(kx)*sin(ky)*sin(kz)
         if (sd == 1) then
           wphi1(1) = -k(1)*sin(kx)*sin(ky)*sin(kz)
           wphi1(2) = +k(2)*cos(kx)*cos(ky)*sin(kz)
           wphi1(3) = +k(3)*cos(kx)*sin(ky)*cos(kz)
           wphi2 = -a*wphi
         end if
        end if
        if (no == 7) then
         wphi = sin(kx)*cos(ky)*sin(kz)
         if (sd == 1) then
           wphi1(1) = +k(1)*cos(kx)*cos(ky)*sin(kz)
           wphi1(2) = -k(2)*sin(kx)*sin(ky)*sin(kz)
           wphi1(3) = +k(3)*sin(kx)*cos(ky)*cos(kz)
           wphi2 = -a*wphi
         end if
        end if
        if (no == 8) then
         wphi = sin(kx)*sin(ky)*sin(kz)
         if (sd == 1) then
           wphi1(1) = +k(1)*cos(kx)*sin(ky)*sin(kz)
           wphi1(2) = +k(2)*sin(kx)*cos(ky)*sin(kz)
           wphi1(3) = +k(3)*sin(kx)*sin(ky)*cos(kz)
           wphi2 = -a*wphi
         end if
        end if
       end subroutine CUBICQD_WAVE

end module