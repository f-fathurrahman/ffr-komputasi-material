PROGRAM debug_H2MOL
! Calculates groundstate of H2 molecule
  use m_highlevel
  use m_midlevel
  use m_random
  use m_jastrow
  use m_orbital

  implicit none
  integer, parameter :: NALPHA=1,NWAVEC=0
  real(dp) :: STEPMAX
  real(dp) :: q,qj,qd,rannumb
  real(dp) :: LOCEN,LOKIN,ERWEN,VAREN,ERWKIN,VARKIN,ERWPOT
  real(dp) :: VARPOT,LOCPOT,LOP1,LOP2,LOK,LOP,LOCWW,ERWWW,VARWW,TOTALEN
  real(dp) :: LOKOLD,LOCLKD,LKDETAIL,LOCVIR,ERWLKD,VIR,ERWVIR,FORCE
  real(dp) :: ERWVELEL,VARVELEL,LOCVEL,VEL,LOCENVEL,ERWENVEL,VARENVEL,TOTENVEL
  real(kind=dp),dimension(NelectronsMax*NatomsMax+1)        :: AVCHA
  real(kind=dp),dimension(NDIVMX,NDIVMX,NDIVMX) :: AVRHO
  
  character(10) :: PRONAME
  ! Local variables
  integer :: i,k,nx,ny,nz,n5,n6, s
  real(dp) :: w,rd
  real(dp) :: minenarr,maxenarr,eminalpha
  real(kind=dp), dimension(NORBM,Nelectrons) :: hpsi

  PRONAME = "H2MOL"
  RANDNAME = "random generator from TAO       "
  ORBNAME = "orbital composition from product        "
  
  write(*,'(1x,2a)') 'Program ', PRONAME
  write(*,'(1x,a)') RANDNAME
  write(*,'(1x,a)') ORBNAME

  ! Number electrons Nelectrons and nuclei Natoms
  NORB = 1
  if( (Nelectrons > NelectronsMax) .or. (Natoms > NatomsMax)) then
    write(*,*) 'Nelectrons or Natoms is too large'
    write(*,*) 'Nelectrons = ', Nelectrons, ' NelectronsMax = ', NelectronsMax
    write(*,*) 'Natoms     = ', Natoms, ' NatomsMax = ', NatomsMax
    stop
  endif

  LENGTH = 10.0d0 ! size of arb. box to display positions

  ! Original values
  !MCPRE = 100000
  !MCMAX = 2000000

  ! For debugging purpose
  MCPRE = 1000  ! 100000
  MCMAX = 20000  ! 2000000

  NDIV = 21        ! NDIV must be odd

  ! Start data:
  ! One nucleus at origin, other nucleus on x-axis
  SWIRHO = .false. ! true if density to be sampled
  SWICHA = .true.  ! true if Madelung charge to be sampled
  CKPOINT = 1.0d0 ! LCAO phase, KONTUZ: double occupation

  ! LCAO for biatomic molecule is a bistable system for large
  ! atom separation; set CKPOINT to zero to enforce atomic limit.
  n5 = 10
  ! Set up atomic coordinates
  DKX = 1.40_dp + (n5 - 10)*0.1_dp       ! distance of both H2 nuclei
  RK(1:3,1:2) = 0.0_dp
  RK(1,2) = DKX

  ! Maximum step width, KONTUZ: Always check with acceptance ratio!
  STEPMAX = 1.0d0
  write(*,'(1x,a,2f12.3)') 'DKX, STEPMAX = ', DKX, STEPMAX

  ! Gaussian localization at nuclei
  BETA1 = 0.01d0
  BETA2 = 0.02d0
  GAM = 0.0001d0
  write(*,'(1x,a,3f12.4)') 'BETA1, BETA2, GAM = ', BETA1, BETA2, GAM

  ! The central Jastrow parameter
  n6 = 10
  CJAS = 0.00001d0 + (n6-10)*1.d0
  write(*,*) 'CJAS = ', CJAS
  
  eminalpha = 0.0d0

  ! Parameter scan ALPHA and WAVEC
  ALPHA = 0.679d0
  WAVEC = 0.1d0 
  write(*,*) 'alpha = ', alpha
  write(*,*) 'wavec = ', wavec

  !
  ! Maximum step width, KONTUZ: Always check with acceptance ratio!
  STEPMAX = 1.0d0
  !
  ! starting always the same sequence of random numbers
  CALL INITRAN()
  !
  ! Random initial electron positions
  do k = 1,Nelectrons
    s = 1 ! what's this?
    if( k > NelectronsPerSpin ) s = 2
    write(*,*)
    write(*,*) 'Before: k = ', k, ' rneu = ', rneu(1:3,k)
    do i = 1,3
      call GENRAN(rannumb)
      rd = (rannumb - 0.5)
      write(*,*) 'k, i = ', k, i
      RE(i,k) = RK(i,k) + rd ! relative from nucleus
      RNEU(i,k) = RE(i,k) ! also assign RNEU <= RE
    enddo
    call ORBWAV( RE(1:3,k), hpsi ) ! evaluate psi at this position
    DOLD(s) = hpsi(1,k) ! save to DOLD
     ! 1= index of orbital
    write(*,*) 'After: k = ', k, ' rneu = ', rneu(1:3,k)
  enddo
  write(*,*) 'DOLD = ', DOLD

  !
  ! Compute initial distances
  VJAS(1:Nelectrons) = 0.0d0
  VJDI(1:Nelectrons) = 0.0d0
  V2POT(1:Nelectrons) = 0.0d0
  V2PDI(1:Nelectrons) = 0.0d0
  !
  do i = 1,Nelectrons
    !
    DISTNEU(1:4,i,i) = 0.0d0 ! idx 1:3 = x,y,z components, idx=4 is the magnitude
    DIST(1:4,i,i) = 0.0d0
    !
    lothers: do k = 1,Nelectrons
      !
      if( k == i ) cycle lothers ! skip itself
      !
      w = 0.d0 ! ????
      !
      DISTNEU(1:3,i,k) = RNEU(1:3,i) - RNEU(1:3,k)
      DIST(1:3,i,k) = DISTNEU(1:3,i,k)
      !
      DISTNEU(4,i,k) = sqrt( sum( (RNEU(1:3,i) - RNEU(1:3,k))**2) )
      DIST(4,i,k) = DISTNEU(4,i,k)
      ! take care below for differing Jastrow factor definitions
      VJAS(i) = VJAS(i) + 1.0d0/DISTNEU(4,i,k)*(1.0d0 - exp(-DISTNEU(4,i,k)/CJAS))
      V2POT(i) = V2POT(i) + 1.0d0/DISTNEU(4,i,k)

      write(*,'(1x,A,2I5,A,4F18.10)') '(', i, k, ') = ', DIST(1:3,i,k), DIST(4,i,k)

    enddo lothers

  enddo ! loop over Nelectrons

  !stop 'DEBUG STOP HERE ...'


  !
  ! Counts the acceptance number
  MCOUNT = 0
  !
  ! Observables
  RHO(1:NDIV,1:NDIV,1:NDIV) = 0._dp
  AVRHO(1:NDIV,1:NDIV,1:NDIV) = 0._dp
  CHA(1:NatomsMax*NelectronsMax+1) = 0._dp
  AVCHA(1:NatomsMax*NelectronsMax+1) = 0._dp
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
  !
  ! MC loop: prerun for thermalizing
  write(*,*)
  write(*,*) 'Begin lpreprun'
  lprerun: do IMC = 1,MCPRE
    !    
    lelpre: do IE = 1,Nelectrons
      ! 
      IES = 1
      if( IE > NelectronsPerSpin ) IES = 2
      !
      do i=1,3
        ! Shift position at random within +-STEPMAX/2
        call GENRAN(rannumb)
        rd = (rannumb - 0.5)*STEPMAX
        RNEU(i,IE) = RE(i,IE) + rd
      enddo
      !
      ! Jastrow factor exponent -0.5*sum_k u_ik without term k=i
      call JEXP(VJAS, VJDI, V2POT, V2PDI)
      qj = exp(-VJDI(IE))
      !
      ! Calculate single particle wavefunction part
      call DETUPD(DNEW(IES), DOLD(IES))
      qd = DNEW(IES)/DOLD(IES)
      ! Test on acceptance
      q = (qd*qj)**2
      if( q < 1.0_dp ) then
        call GENRAN(rannumb)
        MCSCHRITT = (dble(rannumb) < q)
      else
        MCSCHRITT = .true.
      endif
      !
      if (MCSCHRITT) then
        RE(1:3,IE) = RNEU(1:3,IE)
        DOLD(IES) = DNEW(IES)
        MCOUNT = MCOUNT + 1
      else
        RNEU(1:3,IE) = RE(1:3,IE) ! for DETUPD
      endif
    !
    enddo lelpre
  !
  enddo lprerun
  write(*,*)
  write(*,*) 'End of lprerun'
  !
  !
  MCOUNT = 0
  !
  ! MC loop: main run after thermalizing
  write(*,*)
  write(*,*) 'Begin lmainrun'
  write(*,*)
  lmainrun: do IMC=1,MCMAX
    !
    lelmai: do IE=1,Nelectrons
      !
      !if( mod(imc, 10000) == 0 ) then
      !  write(*,*) 'imc = ', imc
      !endif
      !
      IES = 1
      if( IE > NelectronsPerSpin ) IES=2
      !
      do i=1,3
        ! Shift position at random within +-STEPMAX/2
        call GENRAN(rannumb)
        rd = (rannumb-0.5)*STEPMAX
        RNEU(i,IE) = RE(i,IE) + rd
      enddo
      !
      ! Calculate with u_12=CJAS**2/r_12*(1-exp(-r_12/CJAS))*(1.0,0.5)
      ! for (equal,opposite) spin with general
      ! Jastrow factor exponent -0.5*sum_k u_ik without term k=i
      call JEXP(VJAS,VJDI,V2POT,V2PDI)
      qj = exp(-VJDI(IE))
      !
      ! Calculate single particle wavefunction part
      call DETUPD(DNEW(IES), DOLD(IES))
      qd =  DNEW(IES)/DOLD(IES)
      ! Test on acceptance
      q = (qd*qj)**2
      if( q < 1.0_dp ) then
        call GENRAN(rannumb)
        MCSCHRITT = (dble(rannumb) < q)
      else
        MCSCHRITT = .true.
      endif
      !
      ! Update of observables
      if( MCSCHRITT ) then
        RE(1:3,IE) = RNEU(1:3,IE)
        LOP2 = 0.5d0*V2POT(IE)
        call ERGLOC(LOK, LOP1, LOKOLD)
        call VIRIAL(RE(1:3,IE), VIR)
        DOLD(IES) = DNEW(IES)
        VEL = VELEN(IE)
        MCOUNT = MCOUNT + 1
      else
        RNEU(1:3,IE) = RE(1:3,IE) ! necessary for DETUPD
        LOP2 = 0.5_dp*(V2POT(IE) - V2PDI(IE))
        call ERGLOC(LOK,LOP1,LOKOLD)
        call VIRIAL(RNEU(1:3,IE),VIR)
        LOK = LOKOLD
        VEL = VELENOLD(IE)
      endif
      !
      LOP = LOP1 + LOP2
      ! write(*,*)'LOK= ',LOK,'LOP= ',LOP
      !
      ! Factor 0.5 is correct, LOCPOT=0.5 sum_ik v_ik, sum i appears as
      ! loop over electrons IE with contributions that are summed
      ! and divided by Nelectrons, thus energy per electron is calculated
      LOCEN = LOCEN + LOK + LOP
      LOKIN = LOKIN + LOK
      LOCPOT = LOCPOT + LOP1
      LOCWW = LOCWW + LOP2
      LOCLKD = LOCLKD + LKDETAIL
      LOCVIR = LOCVIR + VIR
      LOCVEL = LOCVEL + VEL
      LOCENVEL = LOCENVEL+VEL+LOP
    !
    enddo lelmai
    !
    ! energy per particle
    LOCEN = LOCEN/DBLE(Nelectrons)
    LOKIN = LOKIN/DBLE(Nelectrons)
    LOCPOT = LOCPOT/DBLE(Nelectrons)
    LOCWW = LOCWW/DBLE(Nelectrons)
    LOCLKD = LOCLKD/DBLE(Nelectrons)
    LOCVIR = LOCVIR/DBLE(Nelectrons)
    LOCVEL = LOCVEL/DBLE(Nelectrons)
    LOCENVEL = LOCENVEL/DBLE(Nelectrons)
    ERWEN = DBLE(IMC-1)/DBLE(IMC)*ERWEN + LOCEN/DBLE(IMC)
    ERWKIN = DBLE(IMC-1)/DBLE(IMC)*ERWKIN + LOKIN/DBLE(IMC)
    ERWPOT = DBLE(IMC-1)/DBLE(IMC)*ERWPOT + LOCPOT/DBLE(IMC)
    ERWWW = DBLE(IMC-1)/DBLE(IMC)*ERWWW + LOCWW/DBLE(IMC)
    ERWLKD = DBLE(IMC-1)/DBLE(IMC)*ERWLKD + LOCLKD/DBLE(IMC)
    ERWVIR = DBLE(IMC-1)/DBLE(IMC)*ERWVIR + LOCVIR/DBLE(IMC)
    ERWVELEL = DBLE(IMC-1)/DBLE(IMC)*ERWVELEL + LOCVEL/DBLE(IMC)
    ERWENVEL = DBLE(IMC-1)/DBLE(IMC)*ERWENVEL + LOCENVEL/DBLE(IMC)
    maxenarr = max(maxenarr,LOCEN)
    minenarr = min(minenarr,LOCEN)
    !
    if( IMC > 1 ) then
      VAREN = DBLE(IMC-1)/DBLE(IMC)*VAREN + 1/DBLE(IMC-1)*(ERWEN-LOCEN)**2
      VARKIN = DBLE(IMC-1)/DBLE(IMC)*VARKIN + 1/DBLE(IMC-1)*(ERWKIN-LOKIN)**2
      VARPOT = DBLE(IMC-1)/DBLE(IMC)*VARPOT + 1/DBLE(IMC-1)*(ERWPOT-LOCPOT)**2
      VARWW = DBLE(IMC-1)/DBLE(IMC)*VARWW + 1/DBLE(IMC-1)*(ERWWW-LOCWW)**2
      VARVELEL = DBLE(IMC-1)/DBLE(IMC)*VARVELEL + 1/DBLE(IMC-1)*(ERWVELEL-LOCVEL)**2
      VARENVEL = DBLE(IMC-1)/DBLE(IMC)*VARENVEL + 1/DBLE(IMC-1)*(ERWENVEL-LOCENVEL)**2
    endif
    LOCEN = 0.D0
    LOKIN = 0.D0
    LOCPOT = 0.D0
    LOCWW = 0.D0
    LOCLKD = 0.D0
    LOCVIR = 0.D0
    LOCVEL = 0.D0
    LOCENVEL = 0.D0
    !
    ! Density
    if( SWIRHO ) then
      call DENSITY(RHO)
      AVRHO(1:NDIV,1:NDIV,1:NDIV) = DBLE(IMC-1)/DBLE(IMC) * AVRHO(1:NDIV,1:NDIV,1:NDIV) + &
         & RHO(1:NDIV,1:NDIV,1:NDIV)/DBLE(IMC)
    endif
    !
    ! Madelung charge counting
    if( SWICHA ) then
      call CHARGE(DKX/2.d0,CHA)
      AVCHA(1:Natoms*Nelectrons+1) = DBLE(IMC-1)/DBLE(IMC)*AVCHA(1:Natoms*Nelectrons+1) + &
        & CHA(1:Natoms*Nelectrons+1)/DBLE(IMC)
    endif
    CHA = 0.d0
  !
  enddo lmainrun
  !
  write(*,*)
  write(*,*) 'End of lmainrun'
  write(*,*)
  !
  FORCE = -ERWVIR/DKX - 0.5d0/DKX**2
  ! end MC loop
  !
  write(*,'(1x,a,2e12.3)') 'minenarr,maxenarr=', minenarr, maxenarr
  write(*,'(1x,a,2e14.5)') 'ERWEN, VAREN = ', ERWEN, VAREN
  write(*,'(1x,a,2e14.5)') 'ERWKIN, VARKIN = ', ERWKIN, VARKIN
  write(*,'(1x,a,2e14.5)') 'ERWPOT, VARPOT = ', ERWPOT, VARPOT
  write(*,'(1x,a,2e14.5)') 'ERWWW, VARWW = ', ERWWW, VARWW
  write(*,'(1x,a,2e14.5)') 'ERWVELEL, VARVELEL = ', ERWVELEL, VARVELEL
  write(*,'(1x,a,2e14.5)') 'ERWENVEL, VARENVEL = ', ERWENVEL, VARENVEL
  !
  write(*,'(1x,A,ES14.5)') '2*ERWEN/(ERWPOT+ERWWW) = ', 2._dp*ERWEN/(ERWPOT+ERWWW)
  write(*,'(1x,A,ES14.5)') 'ERWVIR/2*(ERWPOT+ERWWW) = ', ERWVIR/(ERWPOT+ERWWW)
  write(*,'(1x,A,ES14.5)') 'FORCE = ', FORCE
  !
  write(*,'(1x,A,I9)') 'main run: MCMAX = ', MCMAX
  write(*,'(1x,A,F14.5)') 'main run: Acceptance ratio (%) = ', 100.*DBLE(MCOUNT)/DBLE(Nelectrons*MCMAX)
  !
  ! Output density on file
  if( SWIRHO ) then
    open(unit=36,file=PRONAME//"_DENSITY.dat",position="append", status="unknown")
    do nz=1,NDIV
    do ny=1,NDIV
    do nx=1,NDIV
      write(36,'(t3,3i3,e12.3)') nx,ny,nz,AVRHO(nx,ny,nz)
    enddo
    enddo
    enddo
    close(36)
  endif
        !
        !
        TOTALEN = ( 2.d0*ERWEN + 1.d0/DKX + 1.d0 )*HARTREE
        TOTENVEL = ( 2.d0*ERWENVEL + 1.d0/DKX + 1.d0 )*HARTREE
        !
        write(*,'(1x,A,F18.10)') 'ALPHA = ', ALPHA
        write(*,'(1x,A,F18.10)') 'CJAS  = ', CJAS
        write(*,'(1x,A,F18.10)') 'DKX (A) = ', DKX*BOHR
        !
        write(*,'(1x,A,F18.10)') 'TOTALEN (eV) = ', TOTALEN
        write(*,'(1x,A,F18.10)') 'SIGMA (eV)   = ', 2.d0*sqrt(VAREN/MCMAX)*HARTREE
        !
        write(*,'(1x,A,F18.10)') 'TOTENVEL (eV) = ', TOTENVEL
        write(*,'(1x,A,F18.10)') 'SIGMA (eV)   = ', 2.d0*sqrt(VARENVEL/MCMAX)*HARTREE
        !
        write(*,'(1x,A,5F7.3)') 'AVCHA (1,2),(2,1),(1v2,0),(0,1v2),(0,0) = ', AVCHA(1:Natoms*Nelectrons+1)
        Eminalpha = min(Eminalpha, TOTALEN)
        write(*,*)
      !
      write(*,'(1x,A,F18.10)') 'Eminalpha (eV) = ', Eminalpha

end program
