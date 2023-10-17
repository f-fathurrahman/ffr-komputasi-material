module m_random

  use m_highlevel
  use m_midlevel
  ! controls the sequence of random numbers,
  ! allows for use of different random generators
       implicit none
       public :: INITRAN,INITRANNUMB,FILLRAN,FILLG95,FILLTAO,RANF, &
     &           AUTOCORR,AVVAR,BLOCKING, GENRAN

  integer, public :: MCPRE,MCMAX,MAXA,MAXZ


       integer,parameter,public :: MRAN=10001,ISEED=1820459103
       integer,public :: IRAN,IFRAN,ISEED0,ISEED1
       real,dimension(MRAN),public :: FRAN
       character(40), public :: RANDNAME

       real(dp),private :: B=4.294967291D9,C=0.0D0
       real(dp),dimension(43),private :: ZBQLIX
       data  ZBQLIX /8.001441D7,5.5321801D8, &
     &  1.69570999D8,2.88589940D8,2.91581871D8,1.03842493D8, &
     &  7.9952507D7,3.81202335D8,3.11575334D8,4.02878631D8, &
     &  2.49757109D8,1.15192595D8,2.10629619D8,3.99952890D8, &
     &  4.12280521D8,1.33873288D8,7.1345525D7,2.23467704D8, &
     &  2.82934796D8,9.9756750D7,1.68564303D8,2.86817366D8, &
     &  1.14310713D8,3.47045253D8,9.3762426D7 ,1.09670477D8, &
     &  3.20029657D8,3.26369301D8,9.441177D6,3.53244738D8, &
     &  2.44771580D8,1.59804337D8,2.07319904D8,3.37342907D8, &
     &  3.75423178D8,7.0893571D7 ,4.26059785D8,3.95854390D8, &
     &  2.0081010D7,5.9250059D7,1.62176640D8,3.20429173D8, &
     &  2.63576576D8/



contains



subroutine INITRAN
 MCPRE=100000
 MCMAX=20000000
 MAXA=100000
 RANDNAME="random generator from TAO               "
 ISEED0 = ISEED
 CALL FILLRAN
end subroutine INITRAN


subroutine INITRANNUMB
! Set the initial electron positions at random
       integer   :: j,jj
       real(dp)  :: rd
       do j=1,NE
        do jj=1,3
         rd = FRAN(IFRAN)-0.5
         RE(jj,j) = rd
         RENEW(jj,j) = rd
         IFRAN = IFRAN + 1
         if (IFRAN > MRAN) then
          IRAN = IRAN + MRAN
          call FILLRAN
         end if
        end do
       end do
end subroutine INITRANNUMB


subroutine FILLRAN()
 select case (RANDNAME)
  case("random generator from g95               ")
    call FILLG95()
  case("random generator from TAO               ")
    call FILLTAO()
  case default
    call FILLG95()
 end select
end subroutine FILLRAN




! ????? Dead code
subroutine FILLG95()
!       integer :: i
!       ISEED1 = ISEED
!       do i=1,MRAN
! below the g95 random generator which seems ok
!         FRAN(i) = rand(IRAN+i)
!       end do
!       IFRAN = 1
end subroutine FILLG95


subroutine FILLTAO()
  integer :: i
  ISEED1 = ISEED0
  do i=1,MRAN
! below from the book, "An Introduction to Computational Physics"
! written by Tao Pang and published and copyrighted
! by Cambridge University Press in 1997
          FRAN(i)  = RANF()
        end do
        ISEED0 = ISEED1
        IFRAN = 1
end subroutine FILLTAO


       subroutine GENRAN(rannumb)
        real(dp),intent(out) :: rannumb

        real(dp) :: dum

        select case (RANDNAME)
         case("random generator from TAO               ")
           rannumb = FRAN(IFRAN)
           IFRAN = IFRAN + 1
           if (IFRAN .gt. MRAN) THEN
             IRAN=IRAN+MRAN
             call FILLTAO()
           end if
         case("random generator from G95               ")
           rannumb = FRAN(IFRAN)
           IFRAN = IFRAN + 1
           if (IFRAN .gt. MRAN) THEN
             IRAN=IRAN+MRAN
             call FILLG95()
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


SUBROUTINE ZBQLINI(SEED)

!   To initialize the random number generator - either
!   repeatably or nonrepeatably. Need double precision
!   variables because integer storage can't handle the
!   numbers involved
! SEED  (integer, input). User-input number which generates
!   elements of the array ZBQLIX, which is subsequently used
!   in the random number generation algorithm. If SEED=0,
!   the array is seeded using the system clock if the
!   FORTRAN implementation allows it. This feature for
!   nonrepeatable runs has been cancelled here.
! LFLNO (integer). Effect is obsolete here.
!   Number of lowest file handle to try when
!   opening a temporary file to copy the system clock into.
!   Default is 80 to keep out of the way of any existing
!   open files (although the program keeps searching till
!   it finds an available handle). If this causes problems,
!           (which will only happen if handles 80 through 99 are
!           already in use), decrease the default value.
      INTEGER LFLNO
      PARAMETER (LFLNO=80)

! ZBQLIX  Seed array for the random number generator. Defined
!   in ZBQLBD01
! B,C Used in congruential initialisation of ZBQLIX
! FILNO File handle used for temporary file
! INIT  Indicates whether generator has already been initialised

      INTEGER SEED,I
      DOUBLE PRECISION TMPVAR1

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

 !1    FORMAT(//5X,'****WARNING**** You have called routine ZBQLINI ',
 !    +'more than',/5X,'once. I''m ignoring any subsequent calls.',//)
 !2    FORMAT(//5X,'**** ERROR **** In routine ZBQLINI, I couldn''t',
 !    +' find an',/5X,
 !    +'available file number. To rectify the problem, decrease the ',
 !    +'value of',/5X,
 !    +'the parameter LFLNO at the start of this routine (in file ',
 !    +'random_rec_pjn.f)',/5X,
 !    +'and recompile. Any number less than 100 should work.')

end subroutine ZBQLINI


FUNCTION ZBQLU01(DUMMY)

!   Returns a uniform random number between 0 & 1, using
!   a Marsaglia-Zaman type subtract-with-borrow generator.
!   Uses double precision, rather than integer, arithmetic
!   throughout because MZ's integer constants overflow
!   32-bit integer storage (which goes from -2^31 to 2^31).
!   Ideally, we would explicitly truncate all integer
!   quantities at each stage to ensure that the double
!   precision representations do not accumulate approximation
!   error; however, on some machines the use of DNINT to
!   accomplish this is *seriously* slow (run-time increased
!   by a factor of about 3). This double precision version
!   has been tested against an integer implementation that
!   uses long integers (non-standard and, again, slow) -
!   the output was identical up to the 16th decimal place
!   after 10^10 calls, so we're probably OK ...

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

! Update array pointers. Do explicit check for bounds of each to
! avoid expense of modular arithmetic. If one of them is 0 the others
! won't be

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

! The integer arithmetic there can yield X=0, which can cause
! problems in subsequent routines (e.g. ZBQLEXP). The problem
! is simply that X is discrete whereas U is supposed to
! be continuous - hence if X is 0, go back and generate another
! X and return X/B^2 (etc.), which will be uniform on (0,1/B).

      IF (X.LT.BINV) THEN
       B2 = B2*B
       GOTO 5
      ENDIF

      ZBQLU01 = X/B2

      end function ZBQLU01




subroutine AUTOCORR()
        real :: AV=0.0,VA=0.0
        real,dimension(MRAN) :: CORR,AVER,VARI
        integer :: i,k,n
        IRAN = 1
        IFRAN =1
        call FILLG95()
        IRAN = IRAN + MRAN
        do k=1,MRAN-1
          CORR(K) = 0.0
          AV = FRAN(1)
          do i=2,MRAN
            AV = AV*real(i-1)/real(i) + FRAN(i)/real(i)
            VA = VA*real(i-1)/real(i) + (AV-FRAN(i))**2/real(i-1)
            n = modulo(i+k,MRAN)
            CORR(k)=CORR(k)+FRAN(i)*FRAN(n)
          end do
          AVER(k) = AV
          VARI(k) = VA
          CORR(k)=(CORR(k)/real(MRAN)-AV**2)/VA
        end do
        do k=1,MRAN-1
          write(50,*) k,CORR(k),'AV= ',AVER(k),'VA= ',VARI(k)
        end do
end subroutine AUTOCORR


function RANF() result(Z)
        real :: Z
        integer,parameter :: ia=16807,ic=2147483647,iq=127773,ir=2836
        integer :: ih,il,it
! IC = 2**31-1 exists for 32 bit chips
        ih = ISEED1/iq
        il = MOD(ISEED1,iq)
        it = ia*il-ir*ih
        if (it.GT.0) then
          ISEED1 = it
        else
          ISEED1 = ic+it
        end if
        Z = real(ISEED1)/real(ic) ! ffr
end function RANF


subroutine AVVAR(i,f,avgf,varf)
       integer,intent(in)       :: i
       real(dp),intent(in)      :: f
       real(dp),intent(inout)   :: avgf,varf
       if (i == 1) then
        avgf=f
        varf=0.0_dp
       else
        avgf=avgf+(f-avgf)/dble(i)
        varf=dble(i-2)/dble(i-1)*varf+ &
     &       dble(i)/dble(i-1)**2*(avgf-f)**2
       end if
      end subroutine AVVAR


subroutine BLOCKING(fa,avgfa,varfa,avgblocka,varblocka)
       real(dp),intent(in),dimension(NE)   :: fa
       real(dp),intent(inout),dimension(NE)  :: avgfa,varfa,avgblocka, varblocka
! Blocks of length MAXA, where index IMCA counts within a block and
! index IBLOCKA counts the MCMAX/MAXA blocks.
! Input is the random variable fa and output is average and variance:
! at each Monte-Carlo step avgfa,varfa and
! after each block avblocka, varblocka.
! The blocking is used for every electron separately.

    if (IBLOCKA > MCMAX/MAXA) write(*,*) 'Error BLOCKING: Step IMC=', &
       & IMC,' goes beyond last block IBLOCKA=',IBLOCKA, &
       & 'which should be MCMAX/MAXA=',MCMAX/MAXA
    

       if (IMCA == 1) then
        avgfa(IE) = fa(IE)
        varfa(IE) = 0.0_dp
       else
        avgfa(IE)=avgfa(IE)+(fa(IE)-avgfa(IE))/dble(IMCA)
        varfa(IE)=dble(IMCA-2)/dble(IMCA-1)*varfa(IE)+ &
     &       dble(IMCA)/dble(IMCA-1)**2*(avgfa(IE)-fa(IE))**2
       end if
       if (IE == NE) then
        IMCA=IMCA+1
        if (IMCA == MAXA+1) then
! Collect the blocking intermediate results
         if (IBLOCKA == 1) then
          avgblocka=avgfa
          varblocka=0.0_dp
         else
  
          avgblocka=avgblocka+ &
     &               (avgfa-avgblocka)/dble(IBLOCKA)
  
  varblocka=dble(IBLOCKA-2)/dble(IBLOCKA-1)*varblocka+ &
     &    dble(IBLOCKA)/dble(IBLOCKA-1)**2*(avgblocka-avgfa)**2
  
         end if
         IBLOCKA = IBLOCKA + 1
         IMCA = 1
        end if
       end if
      end subroutine BLOCKING

end module