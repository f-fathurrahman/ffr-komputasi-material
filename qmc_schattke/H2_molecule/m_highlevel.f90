!-----------------
module m_highlevel
!-----------------

implicit none

! Double Precision
integer, parameter, public :: dp=selected_real_kind(2*precision(1.0))

! Constants    
real(dp), parameter, public :: EMACH=1.0e-8_dp, PI=3.1415926535897932_dp

! Physical units
real(dp), parameter, public :: HARTREE=27.21168_dp, BOHR=0.52917706_dp

! Number of nuclei
integer, parameter, public :: Natoms = 2

! Number of electrons
integer, parameter, public :: Nelectrons = 2

! Array maxima for nuclei und electrons
integer, parameter, public :: NelectronsMax=2, NatomsMax=2

! NES number electrons per spin
integer, parameter, dimension(2), public :: NES=(/1, 1/)
integer, parameter, public :: NelectronsPerSpin=1 ! no of electron perspin



end module