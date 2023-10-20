!-----------------------------------------------------------------------
      module observables
       use highlevel
       use midlevel
       use random
       use orbital
       use determinant
       use jastrow
!  Comprises subroutines for potential, interaction and kinetic energy
!  INTEN  interaction energy, a specific electron is considered
!  POTEN  potential energy (POTENOLD for the old value before the move,
!         in the following OLD appended as here)
!  KINEN  kinetic energy; consider that this is a two particle energy
!         through the Jastrow factor
!  TOTELEN sum of kinetic,potential and interaction energy specific el.
!  AVELEN  average over full QMC run of TOTELEN
!  VARELEN variance of TOTELEN
!  AVALLEL average of TOTELEN over NE electrons, i.e.total
!         energy of all electrons per electron
!  VARALLEL variance of total energy per electron
!  AVBLEL average of total energy as above, intermediately used only
!         for input and control of blocking
!  VARBLEL variance of total energy per electron
!  AVKINEL average of single electron laplace kinetic energy
!  VARKINEL variance of AVKINEL
!  AVVELEL average of single electron velocity kinetic energy
!  VARVELEL variance of AVVELEL
!  AVPOTEL average of single electron potential energy
!  VARPOTEL variance of AVPOTEL
!  AVINTEL average of single electron part of interaction energy
!  VARINTEL variance of AVINTEL
!  AVBLOCKEL block average of single electron total energy
!  VARBLOCKEL block variance of AVBLOCKEL
       implicit none
       public :: INITOBS,OBSERV,OBSERVSTATEL,OBSERVSTATALL
       public :: INTENERGY,POTENERGY,KINENERGY,DENSITY,PAIRCORR
       integer,parameter,public        :: NRHO=1000
       real(dp),public   :: AVTOTALL,AVALLEL,VARTOTALL,DRHO
       real(dp),public,dimension(NE,2) :: TESTPOT
       real(dp),public,dimension(NE)   :: INTEN,POTEN,KINEN,
     &       POTENOLD,KINENOLD,VELEN,VELENOLD,TOTELEN,
     &       AVELEN,VARELEN,AVBLEL,VARBLEL,AVKINEL,VARKINEL,AVVELEL,
     &       VARVELEL,AVPOTEL,VARPOTEL,AVINTEL,VARINTEL,
     &       AVBLOCKEL,VARBLOCKEL
       real(dp),public,dimension(NRHO,NE)  :: RHORAD,RHORADOLD,AVRHORAD
       real(dp),public,dimension(NRHO,2,2) :: PAIR,AVPAIR
      contains
!-----------------------------------------------------------------------
      subroutine INITOBS
       MCRUN=.false.
       INTEN=EMACH
       POTEN=EMACH
       KINEN=EMACH
       VELEN=EMACH
       TOTELEN=0.0_dp
       DRHO=LCLUSTER/dble(NRHO)
       RHORAD=0.0_dp
       PAIR=0.0_dp
      end subroutine INITOBS
!-----------------------------------------------------------------------
      subroutine INITRANOBS
       MCRUN = .true.
       MCOUNT = 0
       IBLOCKA = 1
       IMCA = 1
       AVELEN=0.0_dp
       AVALLEL=0.0_dp
       VARELEN=0.0_dp
       AVTOTALL=0.0_dp
       VARTOTALL=0.0_dp
       AVKINEL=0.0_dp
       VARKINEL=0.0_dp
       AVBLEL=0.0_dp
       VARBLEL=0.0_dp
       AVBLOCKEL=0.0_dp
       VARBLOCKEL=0.0_dp
       AVRHORAD=0.0_dp
       PAIR=0.0_dp
       AVPAIR=0.0_dp
      end subroutine INITRANOBS
!-----------------------------------------------------------------------
      subroutine OBSERV
!  The actual random values of the observables are called if the step
!  is accepted, otherwise the old values are taken except the
!  interaction energy which could have changed because positions
!  of other electrons have changed after the last access to this
!  electron.
       call INTENERGY(INTEN)
       call KINENERGY(KINEN,VELEN)
       call PAIRCORR
       if(MCSTEP) then
        call POTENERGY(POTEN)
        call DENSITY
       else
        POTEN(IE)=POTENOLD(IE)
        RHORAD(1:NRHO,IE)=RHORADOLD(1:NRHO,IE)
       end if
       TOTELEN(IE)=KINEN(IE)+POTEN(IE)+INTEN(IE)
      end subroutine OBSERV
!-----------------------------------------------------------------------
      subroutine OBSERVSTATEL
!  The statistics for the observables is calculated for each electron
!  separately. The new values are named as old ones for the start of a
!  new IMC step.
         POTENOLD(IE)=POTEN(IE)
         RHORADOLD(1:NRHO,IE)=RHORAD(1:NRHO,IE)
         call AVVAR(IMC,TOTELEN(IE),AVELEN(IE),VARELEN(IE))
         call AVVAR(IMC,KINEN(IE),AVKINEL(IE),VARKINEL(IE))
         call AVVAR(IMC,VELEN(IE),AVVELEL(IE),VARVELEL(IE))
         call AVVAR(IMC,POTEN(IE),AVPOTEL(IE),VARPOTEL(IE))
         call AVVAR(IMC,INTEN(IE),AVINTEL(IE),VARINTEL(IE))
         AVRHORAD(1:NRHO,IE)=AVRHORAD(1:NRHO,IE)*dble(IMC)/dble(IMC+1)+
     &              RHORAD(1:NRHO,IE)/dble(IMC+1)
         call BLOCKING(TOTELEN,AVBLEL,VARBLEL,AVBLOCKEL,VARBLOCKEL)
      end subroutine OBSERVSTATEL
!-----------------------------------------------------------------------
      subroutine OBSERVSTATALL
!  Calculates the overall per electron statistics, i.e. for the
!  average of NE electrons.
       integer  ::  ii
       AVALLEL = sum (TOTELEN(1:NE))/dble(NE)
       call AVVAR(IMC,AVALLEL,AVTOTALL,VARTOTALL)
       AVPAIR(1:NRHO,1:2,1:2) = AVPAIR(1:NRHO,1:2,1:2)*
     &   dble(IMC)/dble(IMC+1) + PAIR(1:NRHO,1:2,1:2)/dble((IMC+1)*NE)
       PAIR = 0.0_dp
      end subroutine OBSERVSTATALL
!-----------------------------------------------------------------------
      subroutine INTENERGY(inten)
!  Calculates interaction energy. The factor 0.5 at the Coulomb double
!  sum is taken into account here, inten(IE)=0.5*\sum_k 1/|r_IE - r_k|.
!  inten(IE)= contribution to the interaction energy of IE-th electron
       real(dp),intent(out),dimension(NE) :: inten
!
       integer :: i,k
       real(dp) :: wo,hw
       hw=0.0_dp
       ielekpot:do k=1,NE
        if (k == IE) then
         cycle ielekpot
        end if
        if (MCSTEP) then
         wo = DISTNEW(4,IE,k)
        else
         wo = DIST(4,IE,k)
        end if
        hw=hw+1.0_dp/wo
       end do ielekpot
       inten(IE)=0.5_dp*hw
      end subroutine INTENERGY
!-----------------------------------------------------------------------
      subroutine POTENERGY(poten)
!   Obsolete for infinite wall cluster
!C  Calculates one-particle potential energy.
!
       real(dp),intent(out),dimension(NE)  :: poten
!C
!       real(dp)  ::  r
!C       r = dsqrt(sum(RENEW(1:3,IE)**2))
!C       r = max (r,JASEMACH)
!C       poten(IE) = -3.0_dp/r
       poten = 0._dp
      end subroutine POTENERGY
!-----------------------------------------------------------------------
      subroutine KINENERGY(kinen,velen)
!  Calculates kinetic energy kinen(IE) of electron IE.
       real(dp),intent(out),dimension(NE)    :: kinen,velen
!
       real(dp)                              :: hkin,qsel
       real(dp)                              :: hlad,hlaj
       real(dp),dimension(3)                 :: hgrd,hgrj
       real(dp),dimension(3,NE)            :: pgrnew
       real(dp),dimension(NE)              :: planew
       if (MCSTEP) then
        qsel = QD(IES)
        call ORBDER(RENEW(1:3,IE),pgrnew,planew)
        hgrj(1:3)=GRJAS(1:3,IEES,IES)
        hlaj=LAPJAS(IEES,IES)
       else
        qsel = 1._dp
        call ORBDER(RE(1:3,IE),pgrnew,planew)
        hgrj(1:3)=GRJASOLD(1:3,IEES,IES)
        hlaj=LAPJASOLD(IEES,IES)
       end if
       call SLAKIN(AOLD(1:IENS,1:IENS,IES),
     &        pgrnew,planew,GRADDET(1:3,1:IENS,IES),
     &                          LAPLDET(1:IENS,IES))
       hgrd(1:3)=GRADDET(1:3,IEES,IES)/qsel
       hlad=LAPLDET(IEES,IES)/qsel
       hkin = hlad+hlaj+2.0_dp*dot_product(hgrd,hgrj)
       velen(IE) = 0.5_dp*sum ((hgrd(1:3)+hgrj(1:3))**2)
       kinen(IE) = -0.5_dp*hkin
      end subroutine KINENERGY
!-----------------------------------------------------------------------
      subroutine DENSITY
!   Radial density array RHORAD is a function of distance s from nucleus and
!   normalized to 1=sum_s 4*PI*s**2 ds RHORAD(s). It is discretized
!   in units of DRHO with NRHO sections. Values below DRHO are added
!   to first unit and those above NRHO*DRHO added to last unit. Remember
!   that box extends over [-LCLUSTER/2,+LCLUSTER/2] with maximum
!   electron distance from center at zero being sqrt(3)*LCLUSTER/2.
       integer    :: j
       real(dp)   :: s,h
       RHORAD = 0.0_dp
       h = 4.0_dp*PI*DRHO
       s=max (sqrt (sum (RENEW(1:3,IE)**2)),DRHO+EMACH)
       s=min (s,NRHO*DRHO)
       j=int (s/DRHO)
       RHORAD(j,IE) = 1._dp/h/s**2
      end subroutine DENSITY
!-----------------------------------------------------------------------
      subroutine PAIRCORR
!  The pair-correlation function is described by PAIR(j,i1,i2) with
!  i1 = IES and i2 = 1 or 2 for both spins of both electrons. Maximum
!  distance between two electrons is sqrt(3)*LCLUSTER, i.e. the mesh
!  width is twice as large as that of the density support, because the 
!  same number of points NRHO is used but the extent is doubled.
       integer    :: ii,j,i1,i2
       real(dp)   :: s,h,dpa
       dpa = 2._dp*DRHO
       h = dpa*(NE-1)
       i1 = IES
       do ii = 1,NE
        if ( ii == IE ) cycle
        i2 = 1
        if ( ii > NES(1) ) i2=2
        s=max (DISTNEW(4,IE,ii),dpa+EMACH)
        s=min (s,NRHO*dpa)
        j=int (s/dpa)
        PAIR(j,i1,i2) = PAIR(j,i1,i2) + 1._dp/h
       end do
      end subroutine PAIRCORR
!-----------------------------------------------------------------------
      end module observables
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
