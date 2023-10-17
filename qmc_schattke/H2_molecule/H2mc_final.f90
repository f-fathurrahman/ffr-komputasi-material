PROGRAM H2MOL
! Calculates groundstate of H2 molecule
      use m_highlevel
      use m_midlevel
      use m_random
      use m_jastrow
      use m_orbital

      implicit none

      integer,parameter       :: NALPHA=10,NWAVEC=0
      real(kind=dp) :: X,Y,Z,STEPMAX
      real(kind=dp) :: DALPHA,ALPHA0,ALPHA1,DWAVEC,WAVEC0,WAVEC1
      real(kind=dp) :: q,qj,qd,rannumb,qd1,qd2,qd3,qd4
      real(kind=dp) :: LOCEN,LOKIN,ERWEN,VAREN,ERWKIN,VARKIN,ERWPOT
      real(dp) :: VARPOT,LOCPOT,LOP1,LOP2,LOK,LOP,LOCWW,ERWWW,VARWW,TOTALEN
  real(dp) :: LOKOLD,LOCLKD,LKDETAIL,LOCVIR,ERWLKD,VIR,ERWVIR,FORCE
  real(dp) :: ERWVELEL,VARVELEL,LOCVEL,VEL,LOCENVEL,ERWENVEL,VARENVEL,TOTENVEL
      real(kind=dp),dimension(NEMAX*NKMAX+1)        :: AVCHA
      real(kind=dp),dimension(NDIVMX,NDIVMX,NDIVMX) :: AVRHO
  
  character(10) :: PRONAME
! Local variables
      integer           :: n1,n2,n3,n4,i,k,n,nx,ny,nz,n5,n6,s
      real(kind=dp)     :: w,rd
      real(kind=dp)     :: RAD1,RAD2,RADNEU,EK,EP,EPD
      real(kind=dp)     :: minenarr,maxenarr,eminalpha
!    erwenx,varenx
      real(kind=dp),dimension(NORBM,NE)    :: hpsi

      PRONAME="H2MOL"
!      RANDNAME="random generator from REC_PJN           "
       RANDNAME="random generator from TAO               "
!      RANDNAME="random generator from G95               "
!      RANDNAME="random generator from F90/95            "
!     ORBNAME ="orbital composition from LCAO           "
!     ORBNAME ="orbital composition from product_p      "
      ORBNAME ="orbital composition from product        "
      write(*,'(1x,2a)') 'Program ',PRONAME
      write(*,'(1x,a)') RANDNAME
      write(*,'(1x,a)') ORBNAME

!  Number electrons NE and nuclei NK
      NORB = 1
      if ((NE.gt.NEMAX) .or. (NK.gt.NKMAX)) then
        write (*,*)'NE or NK= ',NE,NK,' larger than= ',NEMAX,NKMAX
        stop
      end if
      LENGTH = 10.0_dp ! size of arb. box to display positions

!  Number MC steps
      MCPRE = 100000
      MCMAX = 2000000
      NDIV = 21        ! NDIV must be odd

! Start data:
! One nucleus at origin, other nucleus on x-axis
      SWIRHO = .false. ! true if density to be sampled
      SWICHA = .true.  ! true if Madelung charge to be sampled
      CKPOINT = 1.0_dp ! LCAO phase, KONTUZ: double occupation

! LCAO for biatomic molecule is a bistable system for large
! atom separation; set CKPOINT to zero to enforce atomic limit.
      ldkx:do n5=10,10,1
      DKX=1.40_dp+(n5-10)*0.1_dp       ! distance of both H2 nuclei
      RK(1:3,1:2) = 0.0_dp
      RK(1,2) = DKX

! Maximum step width, KONTUZ: Always check with acceptance ratio!
      STEPMAX = 1.00_dp
      write(*,'(1x,a,2f12.3)')'DKX,STEPMAX = ',DKX,STEPMAX

! Gaussian localization at nuclei
      BETA1=0.01_dp
      BETA2=0.02_dp
      GAM=0.0001
      write(*,'(1x,a,3f12.4)')'BETA1,BETA2,GAM = ',BETA1,BETA2,GAM


! The central Jastrow parameter
      ljas:do n6=10,10,1
      CJAS = 0.00001d-0+(n6-10)*1.0_dp
      write(*,'(1x,a,f12.4/)')'CJAS = ',CJAS
      eminalpha = 0._dp

  ! Parameter scan ALPHA and WAVEC
  ! ALPHA0=0.678571428571_dp ! KR
  ALPHA0=0.680_dp
  ! ALPHA0=0.535714285714_dp ! JC
  
  ALPHA1=0.700_dp
  WAVEC0=+0.0_dp
  WAVEC1=+5.0_dp
  DALPHA=(ALPHA1-ALPHA0)/dble(NALPHA)
  DWAVEC=(WAVEC1-WAVEC0)/dble(NWAVEC+1)

  lalpha1:do n1=1,NALPHA+1
  lwavec1:do n3=1,NWAVEC+1


  ALPHA = ALPHA0+(n1-1)*DALPHA
  WAVEC = WAVEC0+(n3-1)*DWAVEC
  
  ! Maximum step width, KONTUZ: Always check with acceptance ratio!
  STEPMAX = 1.00_dp
  
  ! starting always the same sequence of random numbers
  CALL INITRAN()
  
  ! Random initial electron positions
  do k=1,NE
    s=1
    if (k > NES1) s=2
    do i=1,3
      call GENRAN(rannumb)
      rd = (rannumb-0.5)
      RE(i,k) = RK(i,k)+rd
      RNEU(i,k) = RE(i,k)
      call ORBWAV(RE(1:3,k),hpsi)
      DOLD(s)=hpsi(1,k)
    end do
  end do


! Compute initial distances
      VJAS(1:NE) = 0._dp
      VJDI(1:NE) = 0._dp
      V2POT(1:NE) = 0._dp
      V2PDI(1:NE) = 0._dp
      do i=1,NE
      DISTNEU(1:4,i,i)=0._dp
      DIST(1:4,i,i)=0._dp
       lothers:do k=1,NE
        if (k == i) cycle lothers
        w = 0.0_dp
        DISTNEU(1:3,i,k) = RNEU(1:3,i)-RNEU(1:3,k)
        DIST(1:3,i,k) = DISTNEU(1:3,i,k)
        DISTNEU(4,i,k) = dsqrt (sum ((RNEU(1:3,i)-RNEU(1:3,k))**2))
        DIST(4,i,k) = DISTNEU(4,i,k)
!  take care below for differing Jastrow factor definitions
        VJAS(i) = VJAS(i) + 1._dp/DISTNEU(4,i,k)*(1.0_dp-dexp(-DISTNEU(4,i,k)/CJAS))
        V2POT(i) = V2POT(i) + 1._dp/DISTNEU(4,i,k)
       end do lothers
      end do
! Counts the acceptance number
      MCOUNT = 0
! Observables
      RHO(1:NDIV,1:NDIV,1:NDIV) = 0._dp
      AVRHO(1:NDIV,1:NDIV,1:NDIV) = 0._dp
      CHA(1:NKMAX*NEMAX+1) = 0._dp
      AVCHA(1:NKMAX*NEMAX+1) = 0._dp
      LOCEN = 0._dp
      LOKIN = 0._dp
      LOCPOT = 0._dp
      LOCWW = 0._dp
      LOCLKD = 0._dp
      LOCVIR = 0._dp
      LOCVEL = 0._dp
      LOCENVEL = 0._dp
      ERWEN = 0._dp
      VAREN = 0._dp
      ERWKIN = 0._dp
      VARKIN = 0._dp
      ERWPOT = 0._dp
      VARPOT = 0._dp
      ERWWW = 0._dp
      VARWW = 0._dp
      ERWLKD = 0._dp
      ERWVIR = 0._dp
      minenarr = 0._dp
      maxenarr = 0._dp

! MC loop: prerun for thermalizing
      lprerun:do IMC=1,MCPRE
       lelpre:do IE=1,NE
        IES=1
        if (IE > NES1) IES=2
        do i=1,3
! Shift position at random within +-STEPMAX/2
          call GENRAN(rannumb)
          rd = (rannumb-0.5)*STEPMAX
          RNEU(i,IE) = RE(i,IE)+rd
        end do
! Jastrow factor exponent -0.5*sum_k u_ik without term k=i
        call JEXP(VJAS,VJDI,V2POT,V2PDI)
        qj = dexp(-VJDI(IE))
! Calculate single particle wavefunction part
        call DETUPD(DNEW(IES),DOLD(IES))
        qd =  DNEW(IES)/DOLD(IES)
! Test on acceptance
        q = (qd*qj)**2
        if (q < 1.0_dp) then
         call GENRAN(rannumb)
         MCSCHRITT = (dble(rannumb) < q)
        else
         MCSCHRITT = .true.
        end if
        if (MCSCHRITT) then
         RE(1:3,IE) = RNEU(1:3,IE)
         DOLD(IES) = DNEW(IES)
         MCOUNT = MCOUNT + 1
        else
         RNEU(1:3,IE) = RE(1:3,IE) ! for DETUPD
        end if
       end do lelpre
      end do lprerun



  MCOUNT = 0

  ! MC loop: main run after thermalizing
  lmainrun: do IMC=1,MCMAX
  lelmai: do IE=1,NE
    
    !if( mod(imc, 10000) == 0 ) then
    !  write(*,*) 'imc = ', imc
    !endif

    IES=1
    
    if( IE > NES1 ) IES=2
    
    do i=1,3
      ! Shift position at random within +-STEPMAX/2
          call GENRAN(rannumb)
          rd = (rannumb-0.5)*STEPMAX
          RNEU(i,IE) = RE(i,IE)+rd
        end do
! Calculate with u_12=CJAS**2/r_12*(1-exp(-r_12/CJAS))*(1.0,0.5)
! for (equal,opposite) spin with general
! Jastrow factor exponent -0.5*sum_k u_ik without term k=i
        call JEXP(VJAS,VJDI,V2POT,V2PDI)
        qj = dexp(-VJDI(IE))
! Calculate single particle wavefunction part
        call DETUPD(DNEW(IES),DOLD(IES))
        qd =  DNEW(IES)/DOLD(IES)
! Test on acceptance
        q = (qd*qj)**2
        if (q < 1.0_dp) then
         call GENRAN(rannumb)
         MCSCHRITT = (dble(rannumb) < q)
        else
          MCSCHRITT = .true.
        end if
! Update of observables
         if (MCSCHRITT) then
           RE(1:3,IE) = RNEU(1:3,IE)
           LOP2 = 0.5_dp*V2POT(IE)
           call ERGLOC(LOK,LOP1,LOKOLD)
           call VIRIAL(RE(1:3,IE),VIR)
           DOLD(IES) = DNEW(IES)
           VEL = VELEN(IE)
           MCOUNT = MCOUNT + 1
         else
           RNEU(1:3,IE) =  RE(1:3,IE) ! necessary for DETUPD
           LOP2 = 0.5_dp*(V2POT(IE) - V2PDI(IE))
           call ERGLOC(LOK,LOP1,LOKOLD)
           call VIRIAL(RNEU(1:3,IE),VIR)
           LOK=LOKOLD
           VEL = VELENOLD(IE)
         end if
         LOP=LOP1+LOP2
!        write(*,*)'LOK= ',LOK,'LOP= ',LOP
! Factor 0.5 is correct, LOCPOT=0.5 sum_ik v_ik, sum i appears as
! loop over electrons IE with contributions that are summed
! and divided by NE, thus energy per electron is calculated
        LOCEN = LOCEN + LOK + LOP
        LOKIN = LOKIN + LOK
        LOCPOT = LOCPOT + LOP1
        LOCWW = LOCWW + LOP2
        LOCLKD = LOCLKD + LKDETAIL
        LOCVIR = LOCVIR + VIR
        LOCVEL = LOCVEL + VEL
        LOCENVEL = LOCENVEL+VEL+LOP
       end do lelmai
! energy per particle
      LOCEN = LOCEN/DBLE(NE)
      LOKIN = LOKIN/DBLE(NE)
      LOCPOT = LOCPOT/DBLE(NE)
      LOCWW = LOCWW/DBLE(NE)
      LOCLKD = LOCLKD/DBLE(NE)
      LOCVIR = LOCVIR/DBLE(NE)
      LOCVEL = LOCVEL/DBLE(NE)
      LOCENVEL = LOCENVEL/DBLE(NE)
      ERWEN = DBLE(IMC-1)/DBLE(IMC)*ERWEN+LOCEN/DBLE(IMC)
      ERWKIN = DBLE(IMC-1)/DBLE(IMC)*ERWKIN+LOKIN/DBLE(IMC)
      ERWPOT = DBLE(IMC-1)/DBLE(IMC)*ERWPOT+LOCPOT/DBLE(IMC)
      ERWWW = DBLE(IMC-1)/DBLE(IMC)*ERWWW+LOCWW/DBLE(IMC)
      ERWLKD = DBLE(IMC-1)/DBLE(IMC)*ERWLKD+LOCLKD/DBLE(IMC)
      ERWVIR = DBLE(IMC-1)/DBLE(IMC)*ERWVIR+LOCVIR/DBLE(IMC)
      ERWVELEL = DBLE(IMC-1)/DBLE(IMC)*ERWVELEL+LOCVEL/DBLE(IMC)
      ERWENVEL = DBLE(IMC-1)/DBLE(IMC)*ERWENVEL+LOCENVEL/DBLE(IMC)
      maxenarr = max (maxenarr,LOCEN)
      minenarr = min (minenarr,LOCEN)
      if (IMC.gt.1) then
       VAREN = DBLE(IMC-1)/DBLE(IMC)*VAREN + 1/DBLE(IMC-1)*(ERWEN-LOCEN)**2
       VARKIN = DBLE(IMC-1)/DBLE(IMC)*VARKIN + 1/DBLE(IMC-1)*(ERWKIN-LOKIN)**2
       VARPOT = DBLE(IMC-1)/DBLE(IMC)*VARPOT + 1/DBLE(IMC-1)*(ERWPOT-LOCPOT)**2
       VARWW = DBLE(IMC-1)/DBLE(IMC)*VARWW + 1/DBLE(IMC-1)*(ERWWW-LOCWW)**2
       VARVELEL = DBLE(IMC-1)/DBLE(IMC)*VARVELEL + 1/DBLE(IMC-1)*(ERWVELEL-LOCVEL)**2
       VARENVEL = DBLE(IMC-1)/DBLE(IMC)*VARENVEL + 1/DBLE(IMC-1)*(ERWENVEL-LOCENVEL)**2
      end if
      LOCEN = 0.D0
      LOKIN = 0.D0
      LOCPOT = 0.D0
      LOCWW = 0.D0
      LOCLKD = 0.D0
      LOCVIR = 0.D0
      LOCVEL = 0.D0
      LOCENVEL = 0.D0
  ! Density
  if(SWIRHO) then
    call DENSITY(RHO)
    AVRHO(1:NDIV,1:NDIV,1:NDIV) = DBLE(IMC-1)/DBLE(IMC) * AVRHO(1:NDIV,1:NDIV,1:NDIV)+RHO(1:NDIV,1:NDIV,1:NDIV)/DBLE(IMC)
  endif

! Madelung charge counting
        if (SWICHA) then
         call CHARGE(DKX/2.0_dp,CHA)
         AVCHA(1:NK*NE+1) = DBLE(IMC-1)/DBLE(IMC)*AVCHA(1:NK*NE+1) + CHA(1:NK*NE+1)/DBLE(IMC)
        end if
        CHA=0.0_dp
      end do lmainrun

  write(*,*) 'End of lmainrun'

  FORCE = - ERWVIR/DKX - 0.5_dp/DKX**2
  ! end MC loop

  write(*,'(1x,a,2e12.3)')'minenarr,maxenarr=',minenarr,maxenarr
  write(*,'(1x,a,2e14.5)')'ERWEN,VAREN= ',ERWEN,VAREN
  write(*,'(1x,a,2e14.5)')'ERWKIN,VARKIN= ',ERWKIN,VARKIN
  write(*,'(1x,a,2e14.5)')'ERWPOT,VARPOT= ',ERWPOT,VARPOT
  write(*,'(1x,a,2e14.5)')'ERWWW,VARWW= ',ERWWW,VARWW
  write(*,'(1x,a,2e14.5)')'ERWVELEL,VARVELEL= ',ERWVELEL,VARVELEL
  write(*,'(1x,a,2e14.5)')'ERWENVEL,VARENVEL= ',ERWENVEL,VARENVEL

      write(*,'(1x,2a,3e14.5)')'2*ERWEN/(ERWPOT+ERWWW), ERWVIR/2*(ERWPOT+ERWWW), FORCE=', &
     &  2._dp*ERWEN/(ERWPOT+ERWWW),ERWVIR/(ERWPOT+ERWWW),FORCE
      
      write(*,'(1x,a,i9,a14,f7.3,a4)')'main run: MCMAX= ',MCMAX &
     &     ,', acc. ratio = ',100.*DBLE(MCOUNT)/DBLE(NE*MCMAX),' %  '

! Output density on file
      if (SWIRHO) then
        open(unit=36,file=PRONAME//"_DENSITY.dat",position="append", status="unknown")
        do nz=1,NDIV
         do ny=1,NDIV
          do nx=1,NDIV
            write(36,'(t3,3i3,e12.3)') nx,ny,nz,AVRHO(nx,ny,nz)
          end do
         end do
        end do
        close(36)
      end if
      end do lwavec1
      TOTALEN=(2.0_dp*ERWEN+1.0_dp/DKX+1.0_dp)*HARTREE
      TOTENVEL=(2._dp*ERWENVEL+1.0_dp/DKX+1.0_dp)*HARTREE
      write(*,'(1x,a,2e14.3)')'ALPHA,CJAS= ',ALPHA,CJAS
      write(*,'(1x,a,3e14.4)')'DKX(A),TOTALEN(eV),SIGMA(eV)= ',DKX*BOHR, TOTALEN,2._dp*dsqrt(VAREN/MCMAX)*HARTREE
      write(*,'(1x,a,2e14.4)')'TOTENVEL(eV),SIGMA(eV)=',TOTENVEL, 2._dp*dsqrt(VARENVEL/MCMAX)*HARTREE
      write(*,'(1x,2a,5f7.3)') 'AVCHA ', '(1,2),(2,1),(1v2,0),(0,1v2),(0,0)= ',AVCHA(1:NK*NE+1)
      eminalpha = min(eminalpha,TOTALEN)
      write(*,*)
      end do lalpha1
      write(*,'(1x,a,2e14.4)')'E-min_ALPHA (eV) = ',eminalpha

      end do ljas
      end do ldkx

end
