CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  A MC run for the Hydrogen atom
C  a.u. used
C
C
C-----------------------------------------------------------------------
      module random
C  controls the sequence of random numbers,
C  allows for use of different random generators
       implicit none
       public :: FILLG95,FILLTAO,RANF,AUTOCORR,GENRAN
       integer,parameter,public ::
     &   dp=selected_real_kind(2*precision(1.0))
       integer      :: MCMAX,MCPRE,MCOUNT,COUNTD,COUNTU
       integer,parameter,public :: MRAN=10001,ISEED=1820459103
       integer,parameter,public :: MAKF=10000
       integer,public           :: IRAN,IFRAN,ISEED0,ISEED1,SEED
       real(dp),parameter,public   :: PI=3.1415926535897932_dp
       real,dimension(MRAN),public :: FRAN
       character(40) :: RANDNAME
       DOUBLE PRECISION,private :: B=4.294967291D9,C=0.0D0
       DOUBLE PRECISION,dimension(43) :: ZBQLIX
       data  ZBQLIX /8.001441D7,5.5321801D8,
     &  1.69570999D8,2.88589940D8,2.91581871D8,1.03842493D8,
     &  7.9952507D7,3.81202335D8,3.11575334D8,4.02878631D8,
     &  2.49757109D8,1.15192595D8,2.10629619D8,3.99952890D8,
     &  4.12280521D8,1.33873288D8,7.1345525D7,2.23467704D8,
     &  2.82934796D8,9.9756750D7,1.68564303D8,2.86817366D8,
     &  1.14310713D8,3.47045253D8,9.3762426D7 ,1.09670477D8,
     &  3.20029657D8,3.26369301D8,9.441177D6,3.53244738D8,
     &  2.44771580D8,1.59804337D8,2.07319904D8,3.37342907D8,
     &  3.75423178D8,7.0893571D7 ,4.26059785D8,3.95854390D8,
     &  2.0081010D7,5.9250059D7,1.62176640D8,3.20429173D8,
     &  2.63576576D8/
      contains
C
C-----------------------------------------------------------------------
      subroutine INITRAN
       integer                     :: ia
       real(dp)                    :: aver,vari
       real(dp),dimension(MAKF)    :: corr
       real(dp),dimension(MCMAX)   :: rn
       RANDNAME="random generator from REC_PJN           "
C       RANDNAME="random generator from TAO               "
C       RANDNAME="random generator from G95               "
C       RANDNAME="random generator from F90/95            "
       select case (RANDNAME)
         case("random generator from TAO               ")
           write(35,*)'random generator from TAO               '
           ISEED0 = ISEED
           call FILLTAO
         case("random generator from G95               ")
           write(35,*)'random generator from G95               '
           call srand(ISEED)
           call FILLG95
         case("random generator from REC_PJN           ")
           write(35,*)'random generator from REC_PJN           '
           SEED = 11 ! or some other integer
           call ZBQLINI(SEED)
         case("random generator from F90/95            ")
           write(35,*)'random generator from F90/95            '
C  no further initialization for repeatable runs
         case default
           write(35,*)'No random generator! Stop!'
           stop
        end select
C  Test on autocorrelation
C        do ia = 1,MCMAX
C          call GENRAN(rn(ia))
C        end do
C        call AUTOCORR(MCMAX,MAKF,rn,aver,vari,corr)
C        do ia = 1,MAKF
C          write (30,*)ia,corr(ia)
C        end do
      end subroutine INITRAN
C-----------------------------------------------------------------------
       subroutine GENRAN(rannumb)
        DOUBLE PRECISION,intent(out) :: rannumb
C
        DOUBLE PRECISION :: dum
C
        select case (RANDNAME)
         case("random generator from TAO               ")
           rannumb = FRAN(IFRAN)
           IFRAN = IFRAN + 1
           if (IFRAN .gt. MRAN) THEN
             IRAN=IRAN+MRAN
             call FILLTAO
           end if
         case("random generator from G95               ")
           rannumb = FRAN(IFRAN)
           IFRAN = IFRAN + 1
           if (IFRAN .gt. MRAN) THEN
             IRAN=IRAN+MRAN
             call FILLG95
           end if
         case("random generator from REC_PJN           ")
           rannumb = ZBQLU01(dum)
         case("random generator from F90/95            ")
           call random_number(rannumb)
         case default
           write(*,*)'No random generator!'
           stop
        end select
       end subroutine GENRAN
C-----------------------------------------------------------------------
       subroutine FILLG95()
        integer :: i
        do i=1,MRAN
C  below the g95 random generator which does not seem ok
        FRAN(i) = rand()
        end do
        IFRAN = 1
       end subroutine FILLG95
C-----------------------------------------------------------------------
       subroutine FILLTAO()
        integer :: i
        ISEED1 = ISEED0
        do i=1,MRAN
C  below from the book, "An Introduction to Computational Physics"
C  written by Tao Pang and published and copyrighted
C  by Cambridge University Press in 1997
          FRAN(i)  = RANF()
        end do
        ISEED0 = ISEED1
        IFRAN = 1
       end subroutine FILLTAO
C-----------------------------------------------------------------------
       subroutine AUTOCORR(mc,mx,x,aver,vari,corr)
C   mc = no of MC steps
C   mx = no of points for autocorrelation function
        integer,intent(in)                    :: mc,mx
        real(dp),dimension(mx),intent(in)     :: x
        real(dp),intent(out)                  :: aver,vari
        real(dp),dimension(mx), intent(out) :: corr
        integer :: i,k,n
        aver = x(1)
        do n=2,mc
         vari = ((n-2)*vari)/(n-1) + (x(n)-aver)**2/n
         aver = ((n-1)*aver)/n + x(n)/n
        end do
C   correct variance for small number of events
        vari = (vari*mc)/(mc-1)
        write (35,*)'hello from AUTO'
        corr = 0._dp
         do k=1,mx-1
          do i=1,mc-k
            corr(k) = corr(k) + (x(i)-aver)*(x(i+k)-aver)
          end do
          corr(k) = corr(k)/(vari*(mc-1))
         end do
C        write (50,*) 'no. events mx = ',mx,'  average = ', aver,
C     &      '  variance =',vari
C        do k=1,mx-1
C          write (50,*) k,corr(k)
C        end do
       end subroutine AUTOCORR
C-----------------------------------------------------------------------
       function RANF() result(Z)
        real :: Z
        integer,parameter :: ia=16807,ic=2147483647,iq=127773,ir=2836
        integer :: ih,il,it
C  IC = 2**31-1 exists for 32 bit chips
        ih = ISEED1/iq
        il = MOD(ISEED1,iq)
        it = ia*il-ir*ih
        if (it.GT.0) then
          ISEED1 = it
        else
          ISEED1 = ic+it
        end if
        Z = ISEED1/FLOAT(ic)
       end function RANF
C
C-----------------------------------------------------------------------
C      end module random
C***************below from***********************************************
C********   AUTHORS: Richard Chandler       ***********
C********        (richard@stats.ucl.ac.uk)  ***********
C********        Paul Northrop          ***********
C********        (northrop@stats.ox.ac.uk)  ***********
C********   LAST MODIFIED: 26/8/03          ***********
C*******************************************************************
C      module random_rec_pjn
C
C       Initializes seed array etc. for random number generator.
C       The values below have themselves been generated using the
C       NAG generator. The initialization and
C       a subroutine for uniform random numbers in (0,1] is
C       contained.
C
C      integer :: SEED
C      DOUBLE PRECISION :: B=4.294967291D9,C=0.0D0
C      DOUBLE PRECISION,dimension(43) :: ZBQLIX
C      data  ZBQLIX /8.001441D7,5.5321801D8,
c     &  1.69570999D8,2.88589940D8,2.91581871D8,1.03842493D8,
C     &  7.9952507D7,3.81202335D8,3.11575334D8,4.02878631D8,
C     &  2.49757109D8,1.15192595D8,2.10629619D8,3.99952890D8,
C     &  4.12280521D8,1.33873288D8,7.1345525D7,2.23467704D8,
C     &  2.82934796D8,9.9756750D7,1.68564303D8,2.86817366D8,
C     &  1.14310713D8,3.47045253D8,9.3762426D7 ,1.09670477D8,
C     &  3.20029657D8,3.26369301D8,9.441177D6,3.53244738D8,
C     &  2.44771580D8,1.59804337D8,2.07319904D8,3.37342907D8,
C     &  3.75423178D8,7.0893571D7 ,4.26059785D8,3.95854390D8,
C     &  2.0081010D7,5.9250059D7,1.62176640D8,3.20429173D8,
C     &  2.63576576D8/
C
C
C     contains
C******************************************************************
      SUBROUTINE ZBQLINI(SEED)
C       To initialize the random number generator - either
C       repeatably or nonrepeatably. Need double precision
C       variables because integer storage can't handle the
C       numbers involved
C   SEED    (integer, input). User-input number which generates
C       elements of the array ZBQLIX, which is subsequently used
C       in the random number generation algorithm. If SEED=0,
C       the array is seeded using the system clock if the
C       FORTRAN implementation allows it. This feature for
C       nonrepeatable runs has been cancelled here.
C   LFLNO   (integer). Effect is obsolete here.
C       Number of lowest file handle to try when
C       opening a temporary file to copy the system clock into.
C       Default is 80 to keep out of the way of any existing
C       open files (although the program keeps searching till
C       it finds an available handle). If this causes problems,
C               (which will only happen if handles 80 through 99 are
C               already in use), decrease the default value.
C******************************************************************
      INTEGER LFLNO
      PARAMETER (LFLNO=80)
C******************************************************************
C   ZBQLIX  Seed array for the random number generator. Defined
C       in ZBQLBD01
C   B,C Used in congruential initialisation of ZBQLIX
C   FILNO   File handle used for temporary file
C   INIT    Indicates whether generator has already been initialised
C
      INTEGER SEED,SS,MM,HH,DD,FILNO,I
      DOUBLE PRECISION TMPVAR1,DSS,DMM,DHH,DDD
C
      IF (SEED.EQ.0) THEN
       write(*,*)'SEED=0 not allowed in this package'
       stop
      ELSE
       TMPVAR1 = DMOD(DBLE(SEED),B)
      ENDIF
      ZBQLIX(1) = TMPVAR1
      DO 100 I = 2,43
       TMPVAR1 = ZBQLIX(I-1)*3.0269D4
       TMPVAR1 = DMOD(TMPVAR1,B)
       ZBQLIX(I) = TMPVAR1
 100  CONTINUE

 1    FORMAT(//5X,'****WARNING**** You have called routine ZBQLINI ',
     +'more than',/5X,'once. I''m ignoring any subsequent calls.',//)
 2    FORMAT(//5X,'**** ERROR **** In routine ZBQLINI, I couldn''t',
     +' find an',/5X,
     +'available file number. To rectify the problem, decrease the ',
     +'value of',/5X,
     +'the parameter LFLNO at the start of this routine (in file ',
     +'random_rec_pjn.f)',/5X,
     +'and recompile. Any number less than 100 should work.')
C
      end subroutine ZBQLINI
C*****************************************************************
      FUNCTION ZBQLU01(DUMMY)
C
C       Returns a uniform random number between 0 & 1, using
C       a Marsaglia-Zaman type subtract-with-borrow generator.
C       Uses double precision, rather than integer, arithmetic
C       throughout because MZ's integer constants overflow
C       32-bit integer storage (which goes from -2^31 to 2^31).
C       Ideally, we would explicitly truncate all integer
C       quantities at each stage to ensure that the double
C       precision representations do not accumulate approximation
C       error; however, on some machines the use of DNINT to
C       accomplish this is *seriously* slow (run-time increased
C       by a factor of about 3). This double precision version
C       has been tested against an integer implementation that
C       uses long integers (non-standard and, again, slow) -
C       the output was identical up to the 16th decimal place
C       after 10^10 calls, so we're probably OK ...
C
      DOUBLE PRECISION ZBQLU01,DUMMY,X,B2,BINV
      INTEGER CURPOS,ID22,ID43

      SAVE CURPOS,ID22,ID43
      DATA CURPOS,ID22,ID43 /1,22,43/

      B2 = B
      BINV = 1.0D0/B
 5    X = ZBQLIX(ID22) - ZBQLIX(ID43) - C
      IF (X.LT.0.0D0) THEN
       X = X + B
       C = 1.0D0
      ELSE
       C = 0.0D0
      ENDIF
      ZBQLIX(ID43) = X
C
C     Update array pointers. Do explicit check for bounds of each to
C     avoid expense of modular arithmetic. If one of them is 0 the others
C     won't be
C
      CURPOS = CURPOS - 1
      ID22 = ID22 - 1
      ID43 = ID43 - 1
      IF (CURPOS.EQ.0) THEN
       CURPOS=43
      ELSEIF (ID22.EQ.0) THEN
       ID22 = 43
      ELSEIF (ID43.EQ.0) THEN
       ID43 = 43
      ENDIF
C
C     The integer arithmetic there can yield X=0, which can cause
C     problems in subsequent routines (e.g. ZBQLEXP). The problem
C     is simply that X is discrete whereas U is supposed to
C     be continuous - hence if X is 0, go back and generate another
C     X and return X/B^2 (etc.), which will be uniform on (0,1/B).
C
      IF (X.LT.BINV) THEN
       B2 = B2*B
       GOTO 5
      ENDIF
C
      ZBQLU01 = X/B2
C
      end function ZBQLU01
C-----------------------------------------------------------------------
C      end module random_rec_pjn
       end module random
C-----------------------------------------------------------------------
