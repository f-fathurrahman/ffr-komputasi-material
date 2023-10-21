!----------------
module m_position
!----------------

use m_constants, only: pi, dp

implicit none

public :: DENSITY1D

integer, parameter,public          :: NRHO=1000
real(8), parameter,public    :: EMACH=1.0d-8
integer, public :: IMC
real(8), public :: LENGTH,DRHO
real(8), dimension(3), public :: RE,RNEU
real(8), public, dimension(NRHO) :: RHORAD,AVRHORAD

  contains

!---------------------
subroutine DENSITY1D()
!---------------------
  ! Radial density RHORAD is a function of distance s from nucleus and
  ! normalized to 1=sum_s 4*PI*s**2 ds RHORAD(s). It is discretized
  ! in units of DRHO with NRHO sections. Values below DRHO are added
  ! to first unit and those above NRHO*DRHO added to last unit.
  integer :: j
  real(8) :: s,h
  RHORAD = 0.d0
  h = 4.d0*PI*DRHO
  s = max(sqrt(sum(RNEU(1:3)**2)), DRHO+EMACH)
  s = min(s,NRHO*DRHO)
  j = int(s/DRHO)
  RHORAD(j) = 1/h/s**2
end subroutine DENSITY1D

end module

