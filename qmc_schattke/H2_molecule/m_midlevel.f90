!----------------
module m_midlevel
!----------------

use m_highlevel, only: DP, Natoms, Nelectrons
implicit none

SAVE

! IE=1,Nelectrons index of electron
! IES=1,2 index of spin: 1=up, 2=down
! IEES=1,NES(IES) index of electrons of same spin
! IENS=NES(IES) no. of electrons with specified spin
! IK=1,Natoms index of nuclei
! SPINSEL= spin select
integer, public :: IE, IES, IK, IEES, IENS, IMCR, IBLOCKA, SPINSEL
integer, public :: IMC, IMCZ, IMCA, IA, IZ, MCOUNT

! RK position array of nuclei
real(dp), dimension(3,Natoms), public :: RK

! RE actual electron position
! RNEU updated new electron position
real(dp), dimension(3,Nelectrons), public :: RE, RNEU

! QJC acceptance ratio for James and Coolidge series
real(dp) :: QJC ! ffr: change to real(dp)

! MCSCHRITT = .true. if step is accepted
! SWIRHO = .true. if density is calculated
! SWICHA = .true. if Madelung charge is calculated
logical, public :: MCSTEP, MCRUN, MCSCHRITT, SWIRHO, SWICHA


end module

