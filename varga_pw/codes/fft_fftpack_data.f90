!
! $Id: fft_fftpack_data.f90,v 1.6 1999/03/22 11:52:49 noe Exp $
!
module fft_data

  !
  ! fft array storage parameters, fftpack version
  ! 
  ! The fft array will be stored as :
  !   c_ccf( n1 + fftinc1, n2 + fftinc2, n3 + fftinc3 )
  !
  public

  INTEGER             ::  fftinc1 = 1
  INTEGER             ::  fftinc2 = 0
  INTEGER             ::  fftinc3 = 0

  ! fft configuration :
  !   0 ->  no grid adaption necessary
  !   1 ->  T3E
  !
  INTEGER        ::  fft_conf = 0

  !
  ! fft allowed mesh sizes, FFTPACK version
  ! In principle the mesh sizes should be factors of small primes
  ! but we allow all sizes here, i.e. 1st element of the mesh array is 0.
  ! 
  INTEGER , pointer      :: fftmesh(:)

  INTEGER , target       :: fft_defaultmesh(42) =           &
       (/6,8,10,12,14,16,18,20,24,28,30,32,36,40,42,48,    &
         56,60,64,70,72,80,84,90,96,112,120,126,128,140,   &
         144,160,168,180,192,210,224,240,252,256,280,288/)


  interface fft_getmesh
    module procedure getmesh
  end interface

contains

include 'fft_getmesh.f90'

end module
