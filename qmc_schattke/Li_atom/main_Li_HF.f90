program Li_atom
!  QMC main program of Li atom
  use m_highlevel
  use m_midlevel
  use m_random

  use m_orbital_Li_HF
  use m_determinant
  use m_jastrow, JEXP => JASEXPATOM
  use m_observables
  use m_output

       implicit none
       integer          :: i,ii,hi,n1
       real(dp)         :: q,rd,rannumb
!  Open output data files
       call INITOUT
! Fill an array with a set of random numbers to be available;
! contained in M_random.f.
! RANDNAME = name of the random generator
  RANDNAME="random generator from REC_PJN           "
  MCPRE = 100000
  MCMAX = 20000000
  MAXA  = 100000
  write(*,*)'Prerun and main run number of steps, MCPRE,MCMAX=', MCPRE, MCMAX
  call INITRAN()

  ! Set the initial electron positions;
  ! contained in M_random.f
  call INITRANNUMB()

!  Associate electron with single particle wavefunction, normalize
!  wavefunction, and determine initial wavefunction matrix;
!  contained in M_orbital_Li.f
  SLAP(1) = 0.94_dp
  SLAP(2) = 0.93_dp
  write(*,*)'Slater parameter loop, initial and increment, SLAP(1),SLAP(2)=',SLAP(1),SLAP(2)
    
  ljas: do n1=1,15
       CJAS=0.000000001_dp+n1*0.2_dp
       write(*,*)'CJAS=',CJAS
       call INITORB
!  Set initial values for Jastrow parameters
       call INITJAS
!  Set initial values for inverse matrices;
!  contained in M_determinant.f
       call INITDET
!  Set the initial values of the MC run;
!  contained in M_observables.f
       call INITOBS
       STEPMAX = 0.5_dp
       write(*,*)'Maximum step interval= +-0.5*STEPMAX, STEPMAX=', STEPMAX

!
!  MC loop:
  MCOUNT = 0
!  First, the prerun for thermalizing
  lrun:do IMCR=1,MCMAX

    if(mod(imcr,10000) == 0) then
      write(*,*) 'lrun: imcr = ', imcr
    endif

! Start the main run with statistical sampling of the observables if
! the number of prerun steps is reached. Initialize the observables
! and their statistics
  if( IMCR == MCPRE) then
    write(*,*)'prerun: MCPRE= ',MCPRE,' acc. ratio = ', 100.*DBLE(MCOUNT)/DBLE(NE*(MCPRE)),' %  '
    call INITRANOBS()
  end if
  
  if (IMCR >= MCPRE) IMC=IMCR-MCPRE+1

  ! Inner loop on electron index:
  lelrun:do IE=1,NE

    !write(*,*) 'lelrun: ie = ', ie

    IES=1
    IEES=IE
    if (IE > NES(1)) IES=2
    hi=(IES-1)*NES(1)
    IEES=IE-hi
    IENS=NES(IES)
    do i=1,3
  ! Shift position at random within +-STEPMAX/2
         call GENRAN(rannumb)
         rd = (rannumb-0.5_dp)*STEPMAX
         RENEW(i,IE) = RE(i,IE)+rd
        end do
!  Calculate Jastrow factor quantities
        call JEXP
!  Calculate single particle wavefunction part
        call ORBWAV(RENEW(1:3,IE),PSINEW)
        do ii=1,IENS
         PSIMAT(ii,IEES,IES)=PSINEW(NELORB(hi+ii))
        end do
        call SLAQUOT(AOLD(1:IENS,1:IENS,IES), PSIMAT(1:IENS,1:IENS,IES),QD(IES))
!  Test on acceptance
        q = (QD(IES)*QJ(IES))**2
        if (q < 1) then
         call GENRAN(rannumb)
         MCSTEP = (rannumb < q)
        else
         MCSTEP = .true.
        end if
!  Calculates observables
        if (MCRUN) then
         call OBSERV
        end if
!  Update for acceptance or not
        lmcstep:if (MCSTEP) then
         RE(1:3,IE) = RENEW(1:3,IE)
         call SLASM(QD(IES),AOLD(1:IENS,1:IENS,IES), PSIMAT(1:IENS,1:IENS,IES),ANEW(1:IENS,1:IENS,IES))
         AOLD(1:IENS,1:IENS,IES) = ANEW(1:IENS,1:IENS,IES)
         MCOUNT = MCOUNT + 1
        end if lmcstep
! Calculate statistical results of observables per electron
        if (MCRUN) call OBSERVSTATEL
  end do lelrun

!  Calculate statistical results of observables for all electrons
  if (MCRUN) call OBSERVSTATALL()
  
  end do lrun

!  Do the output
!  Write data on files
  call OUTWRITE
!  Write logfile and close data files
  call OUTLOG
  
  end do ljas
  
end program

