!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE MESH
! ===========
!
! This simple module holds the information that describe the KS states:
! the number of orbitals, and the variables that hold them, used by most
! of the programs and subroutines.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module states

  integer :: N_occ    ! Number of occupied orbitals (occupation number is alwasy two)
  integer :: N_empty  ! Number of unoccupied orbitals.
  integer :: N_wf     ! Total (it should always be N_occ + N_wf)

  ! The next two variables hold the wavefunctions; wfs is real, and is used 
  ! by the "gs" program, whereas zwfs is complex and used by the "td" program.
  real(8), allocatable    :: wfs(:, :, :)
  complex(8), allocatable :: zwfs(:, :, :)

end module states
