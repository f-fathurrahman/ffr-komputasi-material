module m_observables
  use m_highlevel
  use m_midlevel
  use m_random
  use m_orbital
  use m_determinant
  use m_jastrow
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
       public :: INTENERGY,POTENERGY,KINENERGY,DENSITY
       integer,parameter,public        :: NRHO=250
       real(dp),public   :: AVTOTALL,AVALLEL,VARTOTALL,DRHO
       real(dp),public,dimension(NE,2) :: TESTPOT
       real(dp),public,dimension(NE)   :: INTEN,POTEN,KINEN, &
     &       POTENOLD,KINENOLD,VELEN,VELENOLD,TOTELEN, &
     &       AVELEN,VARELEN,AVBLEL,VARBLEL,AVKINEL,VARKINEL,AVVELEL, &
     &       VARVELEL,AVPOTEL,VARPOTEL,AVINTEL,VARINTEL, &
     &       AVBLOCKEL,VARBLOCKEL
       real(dp),public,dimension(NRHO,NE) :: RHORAD,RHORADOLD,AVRHORAD
      

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
         AVRHORAD(1:NRHO,IE)=AVRHORAD(1:NRHO,IE)*dble(IMC)/dble(IMC+1)+ &
     &              RHORAD(1:NRHO,IE)/dble(IMC+1)
         call BLOCKING(TOTELEN,AVBLEL,VARBLEL,AVBLOCKEL,VARBLOCKEL)
      end subroutine OBSERVSTATEL
!-----------------------------------------------------------------------


subroutine OBSERVSTATALL
!  Calculates the overall per electron statistics, i.e. for the
!  average of NE electrons.
       AVALLEL = sum (TOTELEN(1:NE))/dble(NE)
       call AVVAR(IMC,AVALLEL,AVTOTALL,VARTOTALL)
      end subroutine OBSERVSTATALL


!-----------------------------------------------------------------------
      subroutine INTENERGY(inten)
!  Calculates interaction energy. The factor 0.5 at the Coulomb double
!  sum is taken into account here, inten(IE)=0.5*\sum_k 1/|r_IE - r_k|.
!  inten(IE)= contribution to the interaction energy of IE-th electron
       real(dp),intent(out),dimension(NE) :: inten
!
       integer :: k
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
       call SLAKIN(AOLD(1:IENS,1:IENS,IES), &
     &        pgrnew,planew,GRADDET(1:3,1:IENS,IES), &
     &                          LAPLDET(1:IENS,IES))
       hgrd(1:3)=GRADDET(1:3,IEES,IES)/qsel
       hlad=LAPLDET(IEES,IES)/qsel
       hkin = hlad+hlaj+2.0_dp*dot_product(hgrd,hgrj)
       velen(IE) = 0.5_dp*sum ((hgrd(1:3)+hgrj(1:3))**2)
       kinen(IE) = -0.5_dp*hkin
      end subroutine KINENERGY
!-----------------------------------------------------------------------
      subroutine DENSITY
!   Radial density RHORAD is a function of distance s from nucleus and
!   normalized to 1=sum_s 4*PI*s**2 ds RHORAD(s). It is discretized
!   in units of DRHO with NRHO sections. Values below DRHO are added
!   to first unit and those above NRHO*DRHO added to last unit.
       integer    :: j
       real(dp)   :: s,h
       RHORAD=0.0_dp
       h=4.0_dp*PI*DRHO
       s=max (sqrt(sum (RENEW(1:3,IE)**2)),DRHO+EMACH)
       s=min (s,NRHO*DRHO)
       j=int (s/DRHO)
       RHORAD(j,IE)=1/h/s**2
      end subroutine DENSITY


end module
