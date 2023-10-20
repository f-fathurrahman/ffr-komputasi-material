program cubic_cluster
!
!  QMC main program of cubic cluster with infinite walls
!
       use m_highlevel
       use m_midlevel
       use m_orbital
       use m_determinant
       use m_random
       use m_jastrow, JEXP=>JASEXPATOM
       use m_observables
       use m_output
!
       implicit none
       integer                :: i,ii
       real(dp)               :: q,rd
       write(*,*)'Length of cubic edge, LCLUSTER = ',LCLUSTER,' a.u.'
       write(*,*)'r_s =',RS
       write(*,*)'Number of electrons, NE = ',NE
!  Open input and output data files
!  Contained in M_output.f
       call INITOUT
!  Read variables if not an input file is present 
!  Fill an array with a set of random numbers to be available;
!  Contained in M_random.f
       call INITRAN
!  Set the initial electron positions;
!  contained in M_random.f
       call INITRANNUMB
!  Associate electron with single particle wavefunction 
!  and determine initial wavefunction matrix;
!  Contained in M_orbital_cubicqd_cluster.f
       call INITORB
!  Set initial values for Jastrow parameters;
!  Contained in M_jastrow.f
       call INITJAS
!  Set initial values for inverse matrices;
!  Contained in M_determinant.f
       call INITDET
!  Set the initial values of the MC run;
!  Contained in M_observables.f
       call INITOBS
       write(*,*)'Maximum step interval= +-0.5*STEPMAX, STEPMAX = ', STEPMAX
       write(*,*)'Prerun and main run number of steps, MCPRE,MCMAX=', MCPRE,MCMAX
!
!  MC loop:
!  First, the prerun for thermalizing
      lrun:do IMCR=1,MCMAX

      if(mod(imcr, 100000) == 0) then
        write(*,*) "imcr ", imcr
      endif
    
    IMC = 0
!  Start the main run with statistical sampling of the observables if
!  the number of prerun steps is reached. Initialize the observables
!  and their statistics
       if ( IMCR == MCPRE) then
        write(*,*)'prerun: MCPRE= ',MCPRE,' acc. ratio = ', 100.*DBLE(MCOUNT)/DBLE(NE*(MCPRE)),' %  '
        call INITRANOBS
       end if
       IMC=IMCR-MCPRE+1
!
!  Inner loop on electron index:
       lelrun:do IE=1,NE
!  Define some spin variables
        IES=1
        if (IE > NES(1)) IES=2
        SPINSEL=(IES-1)*NES(1)
        IEES=IE-SPINSEL
        IENS=NES(IES)
        do i=1,3
!  Shift position at random within +-STEPMAX/2
         rd = (FRAN(IFRAN)-0.5)*STEPMAX
         IFRAN = IFRAN + 1
         if (IFRAN > MRAN) then
          IRAN = IRAN + MRAN
          call FILLRAN
         end if
         RENEW(i,IE) = RE(i,IE)+rd
        end do
         call RDIST(RE,RENEW)
!  Calculate Jastrow factor quantities
        call JEXP
!  Calculate single particle wavefunction part
        call ORBWAV(RENEW(1:3,IE),PSINEW)
!  The matrix PSIMAT of the determinant gets its column IEES
        do ii=1,IENS
         PSIMAT(ii,IEES,IES)=PSINEW(SPINSEL+ii)
        end do
!   Calculate quotient QD, the ratio between new and old determinant,
!   as determinantal acceptance ratio
!   for the offered move of electron IEES of spin IES, calculated by
!   decomposition with respect to the cofactors
        QD(IES) = dot_product (AOLD(IEES,1:IENS,IES),PSIMAT(1:IENS,IEES,IES))
!  Test on acceptance
        q = (QD(IES)*QJ(IES))**2
        if ((dabs(RENEW(1,IE))/LCLUSTER >= 0.5) .or. &
     &      (dabs(RENEW(2,IE))/LCLUSTER >= 0.5) .or. &
     &      (dabs(RENEW(3,IE))/LCLUSTER >= 0.5)) then
         MCSTEP = .false.
        else 
         if (q < 1) then
          MCSTEP = (dble(FRAN(IFRAN)) < q)
          IFRAN = IFRAN + 1
          if (IFRAN .gt. MRAN) then
           IRAN=IRAN+MRAN
           call FILLRAN
          end if
         else
          MCSTEP = .true.
         end if
        end if
!        if (q < 1) then
!         MCSTEP = (dble(FRAN(IFRAN)) < q)
!         IFRAN = IFRAN + 1
!         if (IFRAN .gt. MRAN) then
!          IRAN=IRAN+MRAN
!          call FILLRAN
!         end if
!        else
!         MCSTEP = .true.
!        end if
!        if ((dabs(RENEW(1,IE))/LCLUSTER >= 0.5) .or.
!     &      (dabs(RENEW(2,IE))/LCLUSTER >= 0.5) .or.
!     &      (dabs(RENEW(3,IE))/LCLUSTER >= 0.5)) MCSTEP = .false.
!        if (MCSTEP .and. MCRUN) then
!        write(*,*)'main 2nd:IE,WAV = ',IE,PSINEW
!        write(*,*)'RENEW = ',RENEW(1:3,IE)
!        end if
!  Calculates observables
        if (MCRUN) then
         call OBSERV
        end if
!  Update for acceptance or not
        lmcstep:if (MCSTEP) then
         RE(1:3,IE) = RENEW(1:3,IE)
         call SLASM(QD(IES),AOLD(1:IENS,1:IENS,IES), &
     &         PSIMAT(1:IENS,1:IENS,IES),ANEW(1:IENS,1:IENS,IES))
         AOLD(1:IENS,1:IENS,IES) = ANEW(1:IENS,1:IENS,IES)
         MCOUNT = MCOUNT + 1
        else lmcstep
         RENEW(1:3,IE) = RE(1:3,IE) ! forget the move
        end if lmcstep
!  Calculate statistical results of observables per electron
        if (MCRUN) call OBSERVSTATEL
       end do lelrun
!
!  Calculate statistical results of observables for all electrons
       if (MCRUN) call OBSERVSTATALL
      end do lrun
!
!  Do the output
!  Write density data on file
      call OUTWRITE
!  Write logfile with output data and close data files
      call OUTLOG
      stop

end program
