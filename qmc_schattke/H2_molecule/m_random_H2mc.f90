!--------------
module m_random
!--------------

use m_highlevel
use m_midlevel

implicit none

public :: FILLG95, FILLTAO, RANF, AUTOCORR
public :: GENRAN, INITRAN, ZBQLINI, ZBQLU01

integer, parameter, public :: MRAN=10001,ISEED=1720459103
integer, public :: MCPRE,MCMAX,MAXA,MAXZ
integer, parameter, public :: MAKF=1000
integer, public :: IRAN,IFRAN,ISEED0,ISEED1,SEED
real(dp), dimension(MRAN), public :: FRAN
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


!-------------------
subroutine INITRAN()
!-------------------
  implicit none
!  RANDNAME = name of the random generator, see main program
!   RANDNAME="random generator from REC_PJN           "
!   RANDNAME="random generator from TAO               "
!   RANDNAME="random generator from G95               "
!   RANDNAME="random generator from F90/95            "
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
!  no further initialization for repeatable runs
         case default
           write(35,*)'No random generator! Stop!'
           stop
        end select
!  Test on autocorrelation
!    do ia = 1,MCMAX
!      call GENRAN(rn(ia))
!    end do
!    call AUTOCORR(MCMAX,MAKF,rn,aver,vari,corr)
!    do ia = 1,MAKF
!      write (30,*)ia,corr(ia)
!    end do
end subroutine



!-------------------------
subroutine GENRAN(rannumb)
!-------------------------
  implicit none
  real(dp), intent(out) :: rannumb
  !
  select case (RANDNAME)
  !
  case("random generator from TAO               ")
    rannumb = FRAN(IFRAN)
    IFRAN = IFRAN + 1
    if (IFRAN .gt. MRAN) THEN
      IRAN=IRAN+MRAN
      call FILLTAO()
    end if
  !
  case("random generator from G95               ")
    rannumb = FRAN(IFRAN)
    IFRAN = IFRAN + 1
    if( IFRAN .gt. MRAN ) THEN
      IRAN = IRAN + MRAN
      call FILLG95()
    end if
  !
  case("random generator from REC_PJN           ")
    rannumb = ZBQLU01()
  !
  case("random generator from F90/95            ")
    call random_number(rannumb)
  !
  case default
    write(*,*)'No random generator!'
    stop
  end select

end subroutine



!-------------------
subroutine FILLG95()
!------------------
  implicit none
  integer :: i
  do i=1,MRAN
    FRAN(i) = rand()
  end do
  IFRAN = 1
end subroutine


!-------------------
subroutine FILLTAO()
!-------------------
  implicit none
  integer :: i
  ISEED1 = ISEED0
  do i=1,MRAN
!  below from the book, "An Introduction to Computational Physics"
!  written by Tao Pang and published and copyrighted
!  by Cambridge University Press in 1997
    FRAN(i)  = RANF()
  end do
  ISEED0 = ISEED1
  IFRAN = 1
end subroutine FILLTAO



!------------------------------------------       
subroutine AUTOCORR(mc,mx,x,aver,vari,corr)
!------------------------------------------
  implicit none
! mc = no of MC steps
! mx = no of points for autocorrelation function
  integer,intent(in)                    :: mc,mx
  real(dp),dimension(mx),intent(in)     :: x
  real(dp),intent(out)                  :: aver,vari
  real(dp),dimension(mx), intent(out)   :: corr
  integer :: i,k,n
  aver = x(1)
  do n=2,mc
    vari = ((n-2)*vari)/(n-1) + (x(n)-aver)**2/n
    aver = ((n-1)*aver)/n + x(n)/n
  end do
  ! correct variance for small number of events
  vari = (vari*mc)/(mc-1)
  write (35,*)'hello from AUTO'
  corr = 0._dp
  do k=1,mx-1
    do i=1,mc-k
      corr(k) = corr(k) + (x(i)-aver)*(x(i+k)-aver)
    end do
    corr(k) = corr(k)/(vari*(mc-k-1))
  end do
  ! write (50,*) 'no. events mx = ',mx,'  average = ', aver, '  variance =',vari
  ! do k=1,mx-1
  !   write (50,*) k,corr(k)
  ! end do
end subroutine


!------------------------
function RANF() result(Z)
!------------------------
  implicit none
  real(dp) :: Z
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
  Z = real(ISEED1, kind=dp)/real(ic, kind=dp)
end function RANF


!-----------------------
SUBROUTINE ZBQLINI(SEED)
!-----------------------
  implicit none
!   To initialize the random number generator - either
!   repeatably or nonrepeatably. Need double precision
!   variables because integer storage can't handle the
!   numbers involved
!	SEED	(integer, input). User-input number which generates
!		elements of the array ZBQLIX, which is subsequently used
!		in the random number generation algorithm. If SEED=0,
!		the array is seeded using the system clock if the
!		FORTRAN implementation allows it. This feature for
!   nonrepeatable runs has been cancelled here.
!	LFLNO	(integer). Effect is obsolete here.
!   Number of lowest file handle to try when
!		opening a temporary file to copy the system clock into.
!		Default is 80 to keep out of the way of any existing
!		open files (although the program keeps searching till
!		it finds an available handle). If this causes problems,
!           (which will only happen if handles 80 through 99 are
!           already in use), decrease the default value.
  INTEGER, parameter :: LFLNO=80

!	ZBQLIX	Seed array for the random number generator. Defined
!		in ZBQLBD01
!	B,C	Used in congruential initialisation of ZBQLIX
!	FILNO	File handle used for temporary file
!	INIT	Indicates whether generator has already been initialised

  INTEGER :: SEED, I
  DOUBLE PRECISION :: TMPVAR1

  IF( SEED .EQ. 0) THEN
    write(*,*)'SEED=0 not allowed in this package'
    stop
  ELSE
    TMPVAR1 = DMOD(DBLE(SEED),B)
  ENDIF
  
  ZBQLIX(1) = TMPVAR1
  DO I = 2,43
    TMPVAR1 = ZBQLIX(I-1)*3.0269D4
    TMPVAR1 = DMOD(TMPVAR1,B)
    ZBQLIX(I) = TMPVAR1
  enddo

end subroutine ZBQLINI


!----------------------
FUNCTION ZBQLU01()
!----------------------
  implicit none
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

  DOUBLE PRECISION ZBQLU01,X,B2,BINV
  INTEGER CURPOS,ID22,ID43

  SAVE CURPOS,ID22,ID43
  DATA CURPOS,ID22,ID43 /1,22,43/

  B2 = B
  BINV = 1.0D0/B

5 X = ZBQLIX(ID22) - ZBQLIX(ID43) - C
  IF(X .LT. 0.0D0) THEN
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




end module
