!----------------
module m_midlevel
!----------------

use m_highlevel, only: DP, NK, NE
implicit none

SAVE

! IE=1,NE index of electron
! IES=1,2 index of spin: 1=up, 2=down
! IEES=1,NES(IES) index of electrons of same spin
! IENS=NES(IES) no. of electrons with specified spin
! IK=1,NK index of nuclei
! SPINSEL= spin select
integer, public :: IE, IES, IK, IEES, IENS, IMCR, IBLOCKA, SPINSEL
integer, public :: IMC, IMCZ, IMCA, IA, IZ, MCOUNT

! RK position array of nuclei
real(dp), dimension(3,NK), public :: RK

! RE actual electron position
! RNEU updated new electron position
real(dp), dimension(3,NE), public :: RE,RNEU

! QJC acceptance ratio for James and Coolidge series
real :: QJC

! MCSCHRITT = .true. if step is accepted
! SWIRHO = .true. if density is calculated
! SWICHA = .true. if Madelung charge is calculated
logical, public :: MCSTEP, MCRUN, MCSCHRITT, SWIRHO, SWICHA


end module

