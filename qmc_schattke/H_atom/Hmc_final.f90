! A full QMC run for the H atom
program HMC
  use m_random
  use m_position

  implicit none
  !integer,parameter    :: NXFA1=59,NXFA2=59
  integer,parameter    :: NXFA1=1, NXFA2=1

  real(8)             :: STEPMAX,ALPHA,DALPHA,ALPHA0,ALPHA1
  real(8)             :: WAVEC,DWAVEC,WAVEC0,WAVEC1
  real(8)             :: RDIF,QUOT,LOCEN
  real(8)             :: ERWEN,VAREN,ERWKIN,VARKIN,ERWPOT,VARPOT
  real(8)             :: RAD1,RAD2,RADNEU,EK,EP
  logical              :: MCSCHRITT
  
  ! Local variables
  integer                     :: ix,ih,n1,n2,n3
  real(8)                    :: rannumb
  real(8),dimension(NXFA2+1) :: work1,work2

  ! Controll output
  open(unit=35,FILE="MC21.OUT",STATUS="UNKNOWN")
  open(unit=38,FILE="Hmcvar_splot21.dat",STATUS="UNKNOWN")
  open(unit=39,FILE="Hmcerw_splot21.dat",STATUS="UNKNOWN")
  write(36,*) 'Variance on (c,alpha) - plane'
  write(36,*) NXFA2+1,'  ',NXFA1+1
  write(37,*) 'Energy expectation value on (c,alpha) - plane'
  write(37,*) NXFA2+1,'  ',NXFA1+1
  
  ! Number MC steps
  MCPRE = 1000000
  MCMAX = 1000000
  DRHO = 0.01d0
  
  !  Start data:
  lnxfa1:do n1 = 1,NXFA1+1
  lnxfa2:do n2 = 1,NXFA2+1
    COUNTU = 0
    COUNTD = 0
    ! Wave function coefficients, 0<ALPHA<=2., -1.<=WAVEC<=2.
    ALPHA0 = 0.1d0
    ALPHA1 = 1.9d0
    WAVEC0 = -0.8d0
    WAVEC1 = +1.0d9
    DALPHA=(ALPHA1-ALPHA0)/DBLE(NXFA1+1)
    DWAVEC=(WAVEC1-WAVEC0)/DBLE(NXFA2+1)
    ALPHA=ALPHA0+(n1-1)*DALPHA
    WAVEC=WAVEC0 + (n2-1)*DWAVEC
    
    call INITRAN()

    write(*,*) 'Pass here 51'
  
  ! Initial electron position
       RE(1)=0.1d0
       RE(2)=0.1d0
       RE(3)=0.1d0
        RNEU = RE
!  Maximum step width, KONTUZ: check with acceptance ratio!
        STEPMAX = 2.D-0/(1.D0+ALPHA)
!  Counts the acceptance number
       MCOUNT = 0
!  Observables
       LOCEN = 0.0d0
       ERWEN = 0.0d0
       VAREN = 0.0d0
       ERWKIN = 0.0d0
       VARKIN = 0.0d0
       ERWPOT = 0.0d0
       VARPOT = 0.0d0
       AVRHORAD = 0.0d0
!
!
!  MC loop: prerun for thermalizing,KONTUZ: does not change the bad
!           sampling of energy!!!
       lrunpre:do IMC=1,MCPRE
         do ix=1,3
!  Shift position at random within +-STEPMAX/2
          call GENRAN(rannumb)
          RDIF = (rannumb-0.5)*STEPMAX
          RNEU(ix) = RE(ix)+RDIF
         end do
!  Calculate wave function ratio psi=(1+c*r)exp(-alpha*r)
         RAD1 = SQRT(RE(1)**2+RE(2)**2+RE(3)**2)
         RAD2 = SQRT(RNEU(1)**2+RNEU(2)**2+RNEU(3)**2)
         QUOT = ((1.0_dp+WAVEC*RAD2)/(1.0_dp+WAVEC*RAD1))**2 * EXP(-2.0d0*ALPHA*(RAD2-RAD1))
!  Test on acceptance
         if (QUOT < 1) THEN
           call GENRAN(rannumb)
           MCSCHRITT = dble(rannumb) < QUOT
         else
           MCSCHRITT = .TRUE.
         end if
         if (MCSCHRITT) THEN
             RE = RNEU
             MCOUNT = MCOUNT + 1
         else
             RNEU = RE
         end if
  end do lrunpre

  write(*,*) 'lrunpre is finished'
  
  write(35,*) 'STEPMAX = ',STEPMAX
  write(35,*) 'prerun: MCPRE= ',MCPRE,' acc. ratio = ', 100.*DBLE(MCOUNT)/DBLE(MCPRE),' %  '
  MCOUNT = 0
  COUNTU = 0
  COUNTD = 0

!  MC loop: main run after thermalizing
  lrun:do IMC=1,MCMAX

    if(mod(imc, 10000) == 0) then
      write(*,*) 'Loop imc = ', imc
    endif

    do ix=1,3
!  Shift position at random within +-STEPMAX/2
          call GENRAN(rannumb)
          RDIF = (rannumb-0.5_dp)*STEPMAX
          RNEU(ix) = RE(ix)+RDIF
    end do
!  Calculate wave function ratio psi=(1+c*r)exp(-alpha*r)
         RAD1 = SQRT(RE(1)**2 + RE(2)**2 + RE(3)**2)
         RAD2 = SQRT(RNEU(1)**2 + RNEU(2)**2 + RNEU(3)**2)
         QUOT=((1.0d0 + WAVEC*RAD2)/(1.0d0 + WAVEC*RAD1))**2 * EXP(-2.0d0*ALPHA*(RAD2-RAD1))
!  Test on acceptance
         if (QUOT < 1) THEN
           call GENRAN(rannumb)
           MCSCHRITT = dble(rannumb) < QUOT
           if (MCSCHRITT) COUNTU = COUNTU +1
         else
           MCSCHRITT = .TRUE.
           COUNTD = COUNTD + 1
         end if
         if (MCSCHRITT) THEN
             RE = RNEU
             MCOUNT = MCOUNT + 1
         else
             RNEU = RE
         end if
         RADNEU = DSQRT(RE(1)**2 + RE(2)**2 + RE(3)**2)
!         write (*,*)'RADNEU = ',RADNEU
!  Update of observables
         if (RADNEU .LT. EMACH) THEN
           LOCEN = -0.5_dp*ALPHA**2 + WAVEC**2 + ALPHA*WAVEC + 3.0_dp*(ALPHA-WAVEC)/2.0_dp/EMACH
           EK = 0.0_dp
           EP = 0.0_dp
         else if (DABS(RADNEU*WAVEC+1) .LT. EMACH) THEN
           EK = -0.5_dp*ALPHA**2
           EK = EK + ALPHA - WAVEC*(1.0_dp+2.0_dp*(ALPHA+WAVEC**2))
           EP = -1.0_dp/RADNEU
           LOCEN = EK + EP
         else
           EK = -0.5_dp*ALPHA**2
           EK = EK + ALPHA/RADNEU-WAVEC*(1.0_dp-ALPHA*RADNEU)/(1.0_dp+WAVEC*RADNEU)/RADNEU
           EP = -1.0_dp/RADNEU
           LOCEN = EK + EP
         end if
!  ERWKIN and ERWPOT miss the correction close to the nucleus
         ERWKIN = dble(IMC-1)/dble(IMC)*ERWKIN +EK/dble(IMC)
         ERWPOT = dble(IMC-1)/dble(IMC)*ERWPOT +EP/dble(IMC)
         call DENSITY1D
         AVRHORAD(1:NRHO) = AVRHORAD(1:NRHO)*dble(IMC-1)/dble(IMC) + RHORAD(1:NRHO)/dble(IMC)
         ERWEN = dble(IMC-1)/dble(IMC)*ERWEN+LOCEN/dble(IMC)
         if (IMC == 1) THEN
          VAREN = 0.0_dp
         else
          VAREN = dble(IMC-1)/dble(IMC)*VAREN + 1/dble(IMC-1)*(ERWEN-LOCEN)**2
         end if
  end do lrun

  write(*,*) 'lrun is finished'

  work1(n2) = VAREN
  work2(n2) = ERWEN
  write(35,35)'main run: MCMAX= ',MCMAX,' acc. ratio = ', 100.*dble(MCOUNT)/dble(MCMAX),' %  '
  write(35,*)'downhill steps, towards higher probability, COUNTS in %= ',100.D0*COUNTD/dble(MCMAX),' % '
  write(35,*)'uphill steps, towards lower probability, COUNTS in %= ',100.D0*COUNTU/dble(MCMAX),' % '
  write(35,*)'ALPHA = ',ALPHA,' WAVEC = ',WAVEC
  write(35,*)'energy+0.5*ALPHA**2 = ',ERWEN+0.5*ALPHA**2
  write(35,*)'ERWKIN = ',ERWKIN,'   ERWPOT = ',ERWPOT
  write(35,*)'ERWEN = ',ERWEN,'  VAREN = ',VAREN
  write(35,*)
  write(35,*)
  
  end do lnxfa2

  write(*,*) 'lnxfa2 is finished' 

  
  do n3=1,NXFA2+1
    !  Cut off variance above 0.01 because data plot
    !  gets too large spread in z-values
        if (work1(n3) .GT. 0.01) work1(n3) = 0.05
        write(38,25) ALPHA,WAVEC0 + (n3-1)*DWAVEC,work1(n3)
        write(39,25) ALPHA,WAVEC0 + (n3-1)*DWAVEC,work2(n3)
       end do
  
  end do lnxfa1
  
  write(*,*) 'lnxfa1 is finished' 



  do n3=1,NRHO
       write(47,*) n3,AVRHORAD(n3)
      end do
      close(35)
      close(38)
      close(39)
25    format(1x,2f7.3,e12.4)
35    format(1x,a,i11,a,f7.3,a)

end program
