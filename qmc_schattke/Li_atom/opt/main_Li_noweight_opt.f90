CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCYCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program Li_atom
C  QMC main program of Li atom with parameter optimization
C  Ending OPT hints to the reference wave function used
C  for the random walk. Old parameters refer to test
C  wave functions. Files M_orbital_Li_HF_opt.f, M_jastrow_opt.f,
C  M_output_opt.f.
C  have been introduced for this change. Energies are calculated
C  only for test wave function.
       use highlevel
       use midlevel
       use orbital
       use determinant
       use random
       use jastrow, JEXP=>JASEXPATOM
       use observables
       use output
       implicit none
       integer          :: i,ii,hi,j,jj,islap,n1,nslap1,nslap2
       real(dp)         :: q,rd,rannumb
C  Open output data files
       call INITOUT
C  Fill an array with a set of random numbers to be available;
C  contained in M_random.f.
       MINOPT = 3.0_dp
C                        EGUESS      
C 0.940  0.930  1.000   -2.487579    1.851805 MCMAX=101000
C 0.950  0.935  0.300   -2.512534    1.647485
C 0.945  0.925  1.100   -2.475796    1.538236
C 0.949  0.923  1.300   -2.472657    1.755663
C 0.950  0.924  1.180   -2.488754    1.756638
C 0.949  0.926  1.210   -2.488102    1.632266
C 0.948  0.925  1.225   -2.473639    1.749917
C 0.975  0.915  1.230   -2.487199    1.134109
C 1.000  0.880  1.330   -2.486059    0.659036
C 1.000  0.880  1.150   -2.475601    0.727887
C 1.050  0.835  1.150   -2.472405    0.381200
C 1.150  0.745  1.130   -2.404860    0.108247
C 0.950  1.025  1.050   -2.498521    1.987200 MCMAX=110000
C 0.960  0.835  1.070   -2.491376    2.10175        200000  10h
C 0.950  1.175  0.910   -2.505192    1.976475       110000  20m noweight
C 0.950  1.125  0.910   -2.487895    2.069241          "     "     "
C 0.950  0.725  0.810   -2.513399    2.170238          "     "     "
C 0.950  1.100  0.810   -2.505883    1.841552          "     "     "
C 0.950  0.700  0.810   -2.517777    2.135827          "     "     "
C 0.950  0.950  1.310   -2.495371    1.852254          "     "     "
C 0.950  0.750  0.910   -2.504001    2.056949          "     "     "
C 0.950  0.550  0.610   -2.518231    1.907251          "     "     "
C 0.950  0.500  0.900   -2.542906    1.961285          "     "     "
C 0.950  0.900  0.600   -2.467180    2.664075          "     "     "
C 0.950  0.500  0.550   -2.549482    1.967433          "     "     "
C 0.950  0.800  0.650   -2.431475    2.620268          "     "     "
C 0.950  0.450  1.000   -2.563791    1.924225          "     "     "
C 0.950  0.850  0.650   -2.483594    1.971149          "     "     "
C                         EGUESS
C-------------------------------------------
C 0.950  0.930  1.150   -2.488704   0.024794
C 0.950  0.930  1.300   -2.491726   0.024092
C 0.950  0.980  1.320   -2.492347   0.019977
C 0.950  0.890  1.160   -2.494581   0.017444
C 0.950  0.870  1.300   -2.493565   0.020170
C 0.950  0.910  1.300   -2.492650   0.029421
C 0.950  0.990  1.280   -2.493505   0.023255
C 0.950  0.980  1.320   -2.491471   0.019521
C 0.950  0.890  1.180   -2.492096   0.019686
C 0.950  0.960  1.260   -2.491024   0.027201         2 000 000 12h noweight
C 0.950  0.800  1.120   -2.495546   0.021911           200 000 30m    "
C 0.950  0.900  1.240   -2.491278   0.027951              "     "     "
C 0.950  0.970  1.220   -2.492179   0.021202              "     "     "
C 0.950  0.930  1.280   -2.492673   0.020245              "     "     "
C 0.950  0.840  1.180   -2.491543   0.021130              "     "     "
C 0.950  0.940  1.200   -2.492275   0.020975              "     "     "
C 0.950  0.920  1.260   -2.491518   0.025747              "     "     "
C 0.950  0.970  1.280   -2.491808   0.021371              "     "     "
C 0.950  0.940  1.300   -2.490781   0.022112              "     "     "
C 0.950  0.870  1.100   -2.492614   0.025355         2 000 000  9h    "
C 0.950  0.960  1.260   -2.489977   0.028248           200 000 20m    "
C 0.950  0.870  1.100   -2.491164   0.020883              "     "     "
C 0.950  0.960  1.260   -2.493183   0.021801              "     "     "
C 0.950  0.870  1.100   -2.491164   0.020883              "     "     "
C 0.950  0.960  1.260   -2.493183   0.021801              "     "     "
C A limit cycle with period 2 has established
C 0.950  0.960  1.200   -2.491499   0.026298         2 000 000  9h    "
C 0.950  0.920  1.300   -2.490890   0.024929         2 000 000  9h    "
C 0.950  0.950  1.240   -2.491564   0.025825         2 000 000  9h    "
C 0.950  0.950  1.260   -2.490819   0.024948         2 000 000  9h    "  from reference 0.95,0.95,1.24
C 0.950  0.860  1.180   -2.491980   0.025639         2 000 000  9h    "  from reference 0.95,0.95,1.26
C 0.950  0.970  1.260   -2.490606   0.029184         2 000 000  9h    "  from reference 0.95,0.95,1,26
C 0.950  0.760  1.180   -2.496410   0.027195         2 000 000  9h    "  from reference 0.95,0.95,1,26
C 0.950  0.950  1.320   -2.490632   0.026367         2 000 000  9h    "  from reference 0.95,0.95,1,26 REC_PJN
C from here: EGUESS is input and remains constant during run 
C                         EGUESS     AVTOTALL     sigma**2
C 0.950  0.950  1.260   -2.491604   -2.490632    0.023147 200 000     "  from reference 0.95,0.95,1,26
C 0.950  0.930  1.200   -2.494109   -2.491604    0.018766    "    1h  "  from reference 0.95,0.93,1,20 F90/95
C 0.950  0.950  1.280   -2.490072   -2.494109    0.024280    "    1h  "  from reference 0.95,0.93,1,20 TAO
C 0.950  0.860  1.280   -2.492327   -2.490072    0.023708    "    1h  "  from reference 0.95,0.95,1,28 TAO
C (0.950  0.860  1.280   -2.492327   -2.490072    0.023708    "    1h  "  from reference 0.95,0.95,1,28 TAO?doubled)
C                        AVTOTALL    EGUESS
C 0.950  0.950  1.260   -2.492773   -2.490072    0.022673    "    1h  "  from reference 0.95,0.86,1,28 TAO
C 0.950  0.970  1.240   -2.491770   -2.492773    0.023659    "    1h  "  from reference 0.95,0.95,1,26 TAO
C (0.950  0.880  1.280   -2.491785   -2.491770    0.031466    "    1h  "  from reference 0.95,0.86,1,24 TAO)
C 0.950  0.990  1.280   -2.490125   -2.491770    0.033044    "    1h  "  from reference 0.95,0.97,1,24 TAO
C 0.950  0.930  1.280   -2.491447   -2.490125    0.033143    "    1h  "  from reference 0.95,0.99,1,28 TAO
C 0.950  0.990  1.280   -2.491789   -2.491447    0.027774    "    1h  "  from reference 0.95,0.93,1,28 TAO
C 0.950  0.930  1.280   -2.491447   -2.491789    0.033142    "    1h  "  from reference 0.95,0.99,1,28 TAO
C from here: EGUESS varies during run --------------------------------
C 0.950  0.990  1.280   -2.491789   -2.490285    0.02796    "    1h  "  from reference 0.95,0.93,1,28 TAO
C 0.950  0.930  1.280   -2.491447   -2.492598    0.033176    "    1h  "  from reference 0.95,0.99,1,28 TAO
C 0.950  0.990  1.280   -2.491789   -2.490285    0.027965    "    1h  "  from reference 0.95,0.99,1,28 TAO
C  limit cycle period 2
C from here: 400 000 -------------------------------------------------
C 0.950  0.920  1.280   -2.491528   -2.492051    0.029621    "    2.5h  "  from reference 0.95,0.99,1,28 TAO
C 0.950  1.000  1.240   -2.491554   -2.491293    0.025880    "    2.5h  "  from reference 0.95,0.92,1,28 TAO
C 0.950  0.930  1.280   -2.490272   -2.490395    0.025847    "    2.5h  "  from reference 0.95,1.00,1,24 TAO
C 0.950  0.960  1.260   -2.490688   -2.490900    0.037107    "    2.5h  "  from reference 0.95,0.93,1,28 TAO
C 0.950  0.930  1.280   -2.490585   -2.489941    0.029174    "    2.5h  "  from reference 0.95,0.96,1,26 TAO
C 0.950  0.960  1.260   -2.490688   -2.490900    0.037107    "    2.5h  "  from reference 0.95,0.93,1,28 TAO
C  from here: 50 000 and SLAP1 also floating---------------------------
C 0.990  0.850  1.080   -2.497994   -2.500177    0.020057    "    1.5h  "  from reference 0.95,0.93,1,28 TAO
C 0.980  1.000  1.080   -2.488369   -2.489118    0.017908    "    1.5h  "  from reference 0.99,0.85,1,08 TAO
C 0.980  0.950  1.080   -2.487978   -2.488067    0.016324    "    1.5h  "  from reference 0.98,1.00,1,08 TAO
C 0.980  0.950  1.080   -2.484399   -2.485455    0.017681    "    1.5h  "  from reference 0.98,0.95,1,08 TAO
C converged------------------------------------------------------------
C 0.980  0.950  1.080   -2.489521  0.017425 calculation for this point only at 2 000 000 steps
       EGUESS = -2.491447_dp
C       EGUESS = -2.488_dp
       lslap1: do nslap1=1,10
       lslap2: do nslap2=1,20
C  RANDNAME = name of the random generator
C       RANDNAME="random generator from REC_PJN           "
       RANDNAME="random generator from TAO               "
C       RANDNAME="random generator from G95               "
C       RANDNAME="random generator from F90/95            "
       write(*,*)
       STEPMAX = 0.5_dp
       write(*,*)'Maximum step interval= +-0.5*STEPMAX, STEPMAX=',
     &  STEPMAX
       MCPRE = 10000
       MCMAX = 50000
       MAXA  = 100
       write(*,*)'Prerun and main run number of steps, MCPRE,MCMAX=',
     &                                                 MCPRE,MCMAX
       SLAPOPT(1)=0.980_dp
       SLAPOPT(2)=0.950_dp
       SLAP(1)=0.950_dp + (nslap1-5)*0.01_dp
       SLAP(2)=0.900_dp + (nslap2-10)*0.01_dp
       write(*,*)
       write(*,*)'Slater parameter loop, initial and increment,
     &            SLAP(1),SLAP(2)=',SLAP(1),SLAP(2)
       CJASOPT = 1.08_dp
       ljas: do n1=1,20
       CJAS=0.000000001_dp+1.26_dp+(n1-10)*0.02_dp
       write(*,*)'SLAP(1),SLAP(2),CJAS=',SLAP(1),SLAP(2),CJAS
       write(*,*)
       call INITRAN
C  Set the initial electron positions;
C  contained in M_random.f
       write(*,*) 'check the random numbers to be repeatable'
       call GENRAN(rannumb)
       write(*,*)'after call INITRAN SEED=',SEED,'rannumb=',rannumb
       call INITRANNUMB
       write(*,*)'after call INITRANNUMB: RE(1,1)=',RE(1,1)
C  Associate electron with single particle wavefunction, normalize
C  wavefunction, and determine initial wavefunction matrix;
C  contained in M_orbital_Li.f
       call INITORB
       call INITORB1
C  Set initial values for Jastrow parameters
       call INITJAS
C  Set initial values for inverse matrices;
C  contained in M_determinant.f
       call INITDET
C  Set the initial values of the MC run;
C  contained in M_observables.f
       call INITOBS
       WPSI = 1.0_dp
       SWEIGHT = 0.0_dp
       EGUESS=0.0_dp
C
C  MC loop:
      MCOUNT = 0
C  First, the prerun for thermalizing
      lrun:do IMCR=1,MCMAX
C  Start the main run with statistical sampling of the observables if
C  the number of prerun steps is reached. Initialize the observables
C  and their statistics
       if ( IMCR == MCPRE) then
        write(*,*)'prerun: MCPRE= ',MCPRE,' acc. ratio = ',
     &           100.*DBLE(MCOUNT)/DBLE(NE*(MCPRE)),' %  '
        call INITRANOBS
C         call DETARR(NES(1),PSIMATOPT(1:NES(1),1:NES(1),1),DETOPT(1))
C         call DETARR(NES(1),PSIMAT(1:NES(1),1:NES(1),1),DET(1))
C         DETOPT(2) = PSIMATOPT(1,1,2)
C         DET(2) = PSIMAT(1,1,2)
CC         write(*,*)'JPSI(1),JPSI(2),JPSI(3)=',JPSI(1),JPSI(2),JPSI(3)
CC         write(*,*)'JPSIOPT(1),JPSIOPT(2),JPSIOPT(3)=',
CC     &              JPSIOPT(1),JPSIOPT(2),JPSIOPT(3)
C         WJAS = JPSI(1)*JPSI(2)*JPSI(3)
C         WJASOPT = JPSIOPT(1)*JPSIOPT(2)*JPSIOPT(3)
CC         write(*,*)'WJAS,WJASOPT=',WJAS,WJASOPT
CC         write(*,*)'DET,DETOPT=',DET,DETOPT
       end if
       if (IMCR >= MCPRE) IMC=IMCR-MCPRE+1
C
C  Inner loop on electron index:
       lelrun:do IE=1,NE
        IES=1              ! organize the spin parts
        IEES=IE
        if (IE > NES(1)) IES=2
        hi=(IES-1)*NES(1)
        IEES=IE-hi
        IENS=NES(IES)
        do i=1,3
C  Shift position at random within +-STEPMAX/2
         call GENRAN(rannumb)
         rd = (rannumb-0.5_dp)*STEPMAX
         RENEW(i,IE) = RE(i,IE)+rd
        end do
C        if ((n1 == 5 .and. nslap1 == 5) .and.
C     &      (nslap2 == 4 .or. nslap2 ==6) .and. MCRUN)
C     &       write(40,*)IMC,IE,RENEW(1:3,IE)
C  Calculate Jastrow factor quantities
        call JASEXPATOM1
C  Calculate single particle wavefunction part
        call ORBWAV1(RENEW(1:3,IE),PSINEWOPT)
        do ii=1,IENS
         PSIMATOPT(ii,IEES,IES)=PSINEWOPT(NELORB(hi+ii))
        end do
        call SLAQUOT(AOLDOPT(1:IENS,1:IENS,IES),
     &               PSIMATOPT(1:IENS,1:IENS,IES),QDOPT(IES))
C  Test on acceptance
        q = (QDOPT(IES)*QJOPT(IES))**2
        if (q < 1) then
         call GENRAN(rannumb)
         MCSTEP = (rannumb < q)
        else
         MCSTEP = .true.
        end if
        call JEXP
        call ORBWAV(RENEW(1:3,IE),PSINEW)
        do ii=1,IENS
         PSIMAT(ii,IEES,IES)=PSINEW(NELORB(hi+ii))
        end do
        call SLAQUOT(AOLD(1:IENS,1:IENS,IES),
     &              PSIMAT(1:IENS,1:IENS,IES),QD(IES))
C  Calculates observables
        if (MCRUN) then
         call OBSERV
C         DETOPT(IES) = DETOPT(IES)*QDOPT(IES)
C         DET(IES) = DET(IES)*QD(IES)
C         WJASOPT = WJASOPT*QJOPT(IES)
C         WJAS = WJAS*QJ(IES)
CC         write(*,*)'fracjas=',JPSI(1)*JPSI(2)*JPSI(3)/WJAS
CC         write(*,*)'fracjasopt=',JPSIOPT(1)*JPSIOPT(2)*JPSIOPT(3)
CC     &                           /WJASOPT
C         WJAS = JPSI(1)*JPSI(2)*JPSI(3)
C         WJASOPT = JPSIOPT(1)*JPSIOPT(2)*JPSIOPT(3)
CC         write(*,*)'IES,DET(IES),DETOPT(IES)=',IES,DET(IES),DETOPT(IES)
CC         write(*,*)'IE,WJAS,WJASOPT=',IE,WJAS,WJASOPT
C         if (dabs(DETOPT(IES)*WJASOPT) < EMACH) then
C          WPSI = 1.0_dp
C         else
C          WPSI = (DET(IES)*WJAS/DETOPT(IES)/WJASOPT)**2
CC         WPSI = WPSI*(QJ(IES)*QD(IES)/QJOPT(IES)/QDOPT(IES))**2
C         end if
        end if
C  Update for acceptance or not
        lmcstep:if (MCSTEP) then
         RE(1:3,IE) = RENEW(1:3,IE)
         call SLASM(QDOPT(IES),AOLDOPT(1:IENS,1:IENS,IES),
     &      PSIMATOPT(1:IENS,1:IENS,IES),ANEWOPT(1:IENS,1:IENS,IES))
         AOLDOPT(1:IENS,1:IENS,IES) = ANEWOPT(1:IENS,1:IENS,IES)
         call SLASM(QD(IES),AOLD(1:IENS,1:IENS,IES),
     &         PSIMAT(1:IENS,1:IENS,IES),ANEW(1:IENS,1:IENS,IES))
         AOLD(1:IENS,1:IENS,IES) = ANEW(1:IENS,1:IENS,IES)
         MCOUNT = MCOUNT + 1
        end if lmcstep
C  Calculate statistical results of observables per electron
        if (MCRUN) then
         call OBSERVSTATEL
         IMCOPT = NE*(IMC-1) +IE    ! next lines for optimization
C         write (*,*) 'IMCR,IMC,IE=',IMCR,IMC,IE
C         EGUESS = (IMCOPT-1)*EGUESS/dble(IMCOPT) +
C     &        TOTELEN(IE)/dble(IMCOPT)
C         if ( IMCOPT == 1 ) then
CC          OPTVAR = WPSI*(TOTELEN(1)-EGUESS)**2
C         OPTVAR = 0.0_dp
C         else
C          SWEIGHT = (IMCOPT-1)*SWEIGHT/dble(IMCOPT)
C     &              + WPSI/dble(IMCOPT)
C          OPTVAR = (IMCOPT-1)*OPTVAR/dble(IMCOPT)
C     &             + WPSI*(TOTELEN(IE)-EGUESS)**2/dble(IMCOPT)
C          OPTVAR = (IMCOPT-2)*OPTVAR/dble(IMCOPT-1)
C     &           + (TOTELEN(IE)-EGUESS)**2*IMCOPT/dble(IMCOPT-1)**2
C         end if
C         write (*,*)'IMCOPT,WPSI,SWEIGHT,OPTVAR=',IMCOPT,WPSI,SWEIGHT,
C     &              OPTVAR
        end if
       end do lelrun
C
C  Calculate statistical results of observables for all electrons
       if (MCRUN) then
        call OBSERVSTATALL
         EGUESS = (IMC-1)*EGUESS/dble(IMC) +
     &        AVTOTALL/dble(IMC)
         if ( IMC == 1 ) then
          OPTVAR = 0.0_dp
         else
          OPTVAR = (IMC-2)*OPTVAR/dble(IMC-1)
     &           + (AVALLEL-EGUESS)**2*IMC/dble(IMC-1)**2
         end if
       end if
      end do lrun
C
C  Do the output
C      OPTVAR = OPTVAR/SWEIGHT
      MINOPT = min (OPTVAR,MINOPT)
      write(*,*)'EGUESS,minimum_MSQDeviation=',EGUESS,MINOPT
      write(*,*)'AVTOTALL,VARTOTALL,OPTVAR=',AVTOTALL,VARTOTALL,OPTVAR
C  Write data on files
      call OUTWRITE
C  Write logfile and close data files
C      call OUTLOG
      end do ljas
      end do lslap2
      end do lslap1
      close(31)
      stop
      end program Li_atom
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCYCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
